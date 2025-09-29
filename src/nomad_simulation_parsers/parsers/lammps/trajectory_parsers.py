from typing import Any

import numpy as np
from nomad.parsing.file_parser import FileParser, Quantity, TextParser

from nomad_simulation_parsers.parsers.utils.constants_definitions import (
    CHEMICAL_SYMBOLS,
    REFERENCE_MASSES,
)


class TrajParser(TextParser):
    def __init__(self) -> None:
        self._masses = None
        self._chemical_symbols = None
        super().__init__(None)

    def get_pbc_cell(self, val: str) -> tuple[list, np.ndarray]:
        # TODO: extend logic to handle all LAMMPS-supported pbc styles
        # TODO: collect example outputs!
        # https://docs.lammps.org/boundary.html
        # https://docs.lammps.org/Howto_triclinic.html
        val = val.split()
        cell = np.zeros((3, 3))
        # ! LB edit, untested!
        if 'xy' == val[0]:
            pbc = [v == 'pp' for v in val[3:6]]
            tilt_factors = np.zeros(3)
            for i in range(3):
                tilt_factors[i] = float(val[i * 3 + 8])
                cell[i][i] = float(val[i * 3 + 7]) - float(val[i * 3 + 6])
            xy, yz, xz = tilt_factors
            cell[1][0] = xy
            cell[2][0] = xz
            cell[2][1] = yz
        else:  # TODO: orthogonal can be ff or ss (or mm?)
            pbc = [v == 'pp' for v in val[:3]]
            for i in range(3):
                cell[i][i] = float(val[i * 2 + 4]) - float(val[i * 2 + 3])

        return pbc, cell

    def init_quantities(self) -> None:
        def get_atoms_info(val: str) -> dict[str, float]:
            val = val.split('\n')
            keys = val[0].split()
            values = np.array([v.split() for v in val[1:] if v], dtype=float)
            values = values[values[:, 0].argsort()].T
            return {keys[i]: values[i] for i in range(len(keys))}

        self._quantities = [
            Quantity(
                'time_step',
                r'\s*ITEM:\s*TIMESTEP\s*\n\s*(\d+)\s*\n',
                comment='#',
                repeats=True,
            ),
            Quantity(
                'n_atoms',
                r'\s*ITEM:\s*NUMBER OF ATOMS\s*\n\s*(\d+)\s*\n',
                comment='#',
                repeats=True,
            ),
            Quantity(
                'pbc_cell',
                # TODO: LB - Check why pbc none for atom_run
                r'\s*ITEM: BOX BOUNDS\s*([\s\w]+)\n([\+\-\d\.eE\s]+)\n',
                str_operation=self.get_pbc_cell,
                comment='#',
                repeats=True,
            ),
            Quantity(
                'atoms_info',
                r'\s*ITEM:\s*ATOMS\s*([ \w]+\n)*?([\+\-eE\d\.\n ]+)',
                str_operation=get_atoms_info,
                comment='#',
                repeats=True,
            ),
        ]

    @property
    def with_trajectory(self) -> bool:
        return self.get('atoms_info') is not None

    @property
    def n_frames(self) -> int:
        return len(self.get('atoms_info', []))

    @property
    def masses(self) -> np.ndarray:
        return self._masses

    # TODO: handle non-atomistic representations
    @masses.setter
    def masses(self, val: Any) -> None:
        if not val:
            return
        if not isinstance(val, np.ndarray):
            try:
                val = np.asarray(val)
            except (ValueError, TypeError) as e:
                raise ValueError(f'Cannot convert masses to array: {e}')
        min_dim = 2
        if val.ndim < min_dim:
            raise ValueError(f'masses must be at least 2D, got {val.ndim}D array')

        self._masses = val

    @property
    def chemical_symbols(self) -> dict | None:
        """Chemical symbols derived from particle masses."""
        if self._chemical_symbols is None and self._masses is not None:
            self._derive_chemical_symbols()
        return self._chemical_symbols

    def _derive_chemical_symbols(self) -> None:
        """Derive chemical symbols from mass data."""
        masses = self._masses[0][1]
        self._chemical_symbols = {}
        for i in range(len(masses)):
            symbol_idx = np.argmin(abs(REFERENCE_MASSES - masses[i][1]))
            self._chemical_symbols[masses[i][0]] = CHEMICAL_SYMBOLS[symbol_idx]

    def get_atom_labels(self, idx: int) -> list[str] | None:
        atoms_info = self.get('atoms_info')
        if atoms_info is None:
            return
        atoms_info = atoms_info.get(idx)

        atoms_id = atoms_info.get('id')
        default = ['CGX' for _ in atoms_id] if atoms_id is not None else None
        atoms_type = atoms_info.get('type')
        if atoms_type is None:
            return default
        if self._chemical_symbols is None:
            return default

        atom_labels = [self._chemical_symbols[atype] for atype in atoms_type]

        return [label for label in atom_labels if label is not None]

    def _get_frame_atoms_info(self, idx: int) -> dict | None:
        """Helper to get atoms_info for a specific frame."""
        atoms_info = self.get('atoms_info')
        # ? Is atoms_info really a list, or is it a dict with frame indices as keys?
        if atoms_info is None or idx < 0 or idx >= len(atoms_info):
            return None
        return atoms_info[idx]

    def _extract_vector_components(
        self, atoms_info: dict, *keys: str
    ) -> np.ndarray | None:
        """Extract vector components from atoms_info."""
        print(keys)
        if all(k in atoms_info for k in keys):
            return np.transpose([atoms_info.get(key) for key in keys])
        return None

    def get_positions(self, idx: int) -> np.ndarray | None:
        frame_atoms_info = self._get_frame_atoms_info(idx)
        if frame_atoms_info is None:
            return None

        positions = None

        def has_coords(*keys):
            return all(k in frame_atoms_info for k in keys)

        cell = self.get('pbc_cell')
        cell = None if cell is None else cell[idx][1]

        # Cell required
        if cell is not None:
            if has_coords('xs', 'ys', 'zs'):
                positions = self._extract_vector_components(
                    frame_atoms_info, 'xs', 'ys', 'zs'
                )
                positions = positions * np.linalg.norm(cell, axis=1) + np.amin(
                    cell, axis=1
                )
            elif has_coords('xsu', 'ysu', 'zsu'):
                positions = self._extract_vector_components(
                    frame_atoms_info, 'xsu', 'ysu', 'zsu'
                )
                positions = positions * np.linalg.norm(cell, axis=1) + np.amin(
                    cell, axis=1
                )

        # Unwrapped
        if positions is None and has_coords('xu', 'yu', 'zu'):
            positions = self._extract_vector_components(
                frame_atoms_info, 'xu', 'yu', 'zu'
            )

        # Absolute positions with optional image correction
        if positions is None and has_coords('x', 'y', 'z'):
            positions = self._extract_vector_components(frame_atoms_info, 'x', 'y', 'z')
            if cell is not None and has_coords('ix', 'iy', 'iz'):
                positions_img = self._extract_vector_components(
                    frame_atoms_info, 'ix', 'iy', 'iz'
                )
                positions += positions_img * np.linalg.norm(cell, axis=1)

        return positions

    def get_velocities(self, idx: int) -> np.ndarray | None:
        frame_atoms_info = self._get_frame_atoms_info(idx)
        if frame_atoms_info is None:
            return None
        return self._extract_vector_components(frame_atoms_info, 'vx', 'vy', 'vz')

    def get_forces(self, idx: int) -> np.ndarray | None:
        frame_atoms_info = self._get_frame_atoms_info(idx)
        if frame_atoms_info is None:
            return None
        return self._extract_vector_components(frame_atoms_info, 'fx', 'fy', 'fz')

    def get_lattice_vectors(self, idx: int) -> np.ndarray | None:
        pbc_cell = self.get('pbc_cell')
        if pbc_cell is None:
            return
        return pbc_cell[idx][1]

    def get_pbc(self, idx: int) -> list[bool] | None:
        pbc_cell = self.get('pbc_cell')
        if pbc_cell is None:
            return
        return pbc_cell[idx][0]

    def get_n_atoms(self, idx: int) -> int | None:
        n_atoms = self.get('n_atoms')
        if n_atoms is None:
            return len(self.get_positions(idx))
        return n_atoms[idx]

    def get_step(self, idx: int) -> int | None:
        step = self.get('time_step')
        if step is None:
            return
        return step[idx]


class XYZTrajParser(TrajParser):
    def __init__(self) -> None:
        super().__init__()

    def init_quantities(self) -> None:
        def get_atoms_info(val_in: str) -> dict[str, int | float]:
            val = [v.split('#')[0].split() for v in val_in.strip().split('\n')]
            symbols = []
            for v in val:
                if v[0].isalpha():
                    if v[0] not in symbols:
                        symbols.append(v[0])
                    v[0] = symbols.index(v[0]) + 1
            N_VALS = 4
            val = np.transpose(
                np.array([v for v in val if len(v) == N_VALS], dtype=float)
            )
            # val[0] is the atomic number
            val[0] = [list(set(val[0])).index(v) + 1 for v in val[0]]
            return dict(type=val[0], x=val[1], y=val[2], z=val[3])

        self.quantities = [
            Quantity(
                'atoms_info',
                r'((?:\d+|[A-Z][a-z]?) [\s\S]+?)(?:\s\d+\n|\Z)',
                str_operation=get_atoms_info,
                comment='#',
                repeats=True,
            )
        ]


class TrajParsers:
    def __init__(self, parsers: list[TextParser | FileParser]) -> None:
        self._parsers = parsers
        for parser in parsers:
            parser.parse()

    # ? Is this function used anywhere?
    # ? Should the return types for the else case be handled better?
    def __getitem__(self, index) -> TrajParser | None:
        if self._parsers:
            return self._parsers[index]

    # ? Also here, should we make the negative return type more explicit?
    def eval(self, key: str, *args, **kwargs) -> Any | None:
        for parser in self._parsers:
            parser_method = getattr(parser, key)
            if parser_method is not None:
                val = (
                    parser_method(*args, **kwargs) if args or kwargs else parser_method
                )
                if val is not None:
                    return val
