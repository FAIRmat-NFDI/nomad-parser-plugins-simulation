from typing import Any

import numpy as np
from nomad_file_parser import Quantity, TextParser

from nomad_simulation_parsers.parsers.utils.constants import (
    CHEMICAL_SYMBOLS,
    REFERENCE_MASSES,
)
from nomad_simulation_parsers.parsers.utils.mdanalysisparser import MDAnalysisParser


class TrajParser(TextParser):
    """
    Parser for LAMMPS trajectory dump files.

    Extracts positions, velocities, forces, and cell information from file
    specified by LAMMPS 'dump' command.
    """

    def __init__(self) -> None:
        self._masses = None
        self._chemical_symbols = None
        # TODO: add self._atom_type_overrides: dict[int, int] = {} here
        # once virtual-type support is implemented.  Maps atom_id to a
        # synthetic type_id whose chemical_symbol holds the overridden label.
        super().__init__(None)

    def reset(self):
        super().reset()
        self._masses = None
        self._chemical_symbols = None
        # TODO: reset self._atom_type_overrides = {} here as well.

    def get_pbc_cell(self, val: str) -> tuple[list, np.ndarray]:
        # TODO: extend logic to handle all LAMMPS-supported pbc styles
        # TODO: collect example outputs!
        # https://docs.lammps.org/boundary.html
        # https://docs.lammps.org/Howto_triclinic.html
        val = val.split()
        cell = np.zeros((3, 3))
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
        if val is None:
            return
        if not isinstance(val, np.ndarray):
            try:
                val = np.asarray(val)
            except (ValueError, TypeError):
                return None
        min_dim = 2
        if val.ndim < min_dim:
            return None

        self._masses = val

    @property
    def chemical_symbols(self) -> dict | None:
        """Chemical symbols derived from particle masses."""
        if self._chemical_symbols is None and self._masses is not None:
            self._chemical_symbols = {}
            for i in range(len(self._masses)):
                symbol_idx = np.argmin(abs(REFERENCE_MASSES - self._masses[i][1]))
                self._chemical_symbols[self._masses[i][0]] = CHEMICAL_SYMBOLS[
                    symbol_idx
                ]
        return self._chemical_symbols

    # TODO: Inspect other parsers for consistent default return type
    def get_atom_labels(self, idx: int) -> list[str] | list:
        atoms_info = self.get('atoms_info')
        try:
            atoms_info = atoms_info[idx]
        except (TypeError, IndexError, AttributeError):
            return []

        atoms_id = atoms_info.get('id')
        default = ['CGX' for _ in atoms_id] if atoms_id is not None else []
        atoms_type = atoms_info.get('type')

        if atoms_type is None:
            return default
        # Access property to trigger lazy computation from masses
        if self.chemical_symbols is None:
            return default

        # TODO: when virtual-type support is implemented, resolve per-atom
        # type overrides here before the chemical_symbols lookup:
        #   effective_types = [
        #       self._atom_type_overrides.get(int(aid), int(atype))
        #       for aid, atype in zip(atoms_id, atoms_type)
        #   ]
        # then replace atoms_type with effective_types below.
        atom_labels = [self._chemical_symbols[atype] for atype in atoms_type]

        return [
            'CGX' if label == 'X' else label
            for label in atom_labels
            if label is not None
        ]

    def _get_frame_atoms_info(self, idx: int) -> dict:
        """Helper to get atoms_info for a specific frame."""
        atoms_info = self.get('atoms_info')
        if atoms_info is None or idx >= len(atoms_info):
            return {}
        return atoms_info[idx]

    def _extract_vector_components(
        self, atoms_info: dict, *keys: str
    ) -> np.ndarray | None:
        """Extract vector components from atoms_info."""
        if all(k in atoms_info for k in keys):
            return np.transpose([atoms_info.get(key) for key in keys])
        return None

    def get_positions(self, idx: int) -> np.ndarray | None:
        frame_atoms_info = self._get_frame_atoms_info(idx)
        positions = None

        def has_coords(*keys):
            return all(k in frame_atoms_info for k in keys)

        cell = self.get('pbc_cell')
        cell = None if cell is None else cell[idx][1]

        # Cell required
        if cell is not None:
            for coord_set in [('xs', 'ys', 'zs'), ('xsu', 'ysu', 'zsu')]:
                if has_coords(*coord_set):
                    positions = self._extract_vector_components(
                        frame_atoms_info, *coord_set
                    )
                    if positions is not None:
                        # Apply cell scaling transformation
                        positions = positions * np.linalg.norm(cell, axis=1) + np.amin(
                            cell, axis=1
                        )
                        break

        # Unwrapped
        if positions is None and has_coords('xu', 'yu', 'zu'):
            positions = self._extract_vector_components(
                frame_atoms_info, 'xu', 'yu', 'zu'
            )

        # Absolute positions with optional image correction
        if positions is None and has_coords('x', 'y', 'z'):
            positions = self._extract_vector_components(frame_atoms_info, 'x', 'y', 'z')
            if (
                cell is not None
                and positions is not None
                and has_coords('ix', 'iy', 'iz')
            ):
                positions_img = self._extract_vector_components(
                    frame_atoms_info, 'ix', 'iy', 'iz'
                )
                if positions_img is not None:
                    positions += positions_img * np.linalg.norm(cell, axis=1)

        return positions

    def get_velocities(self, idx: int) -> np.ndarray | None:
        frame_atoms_info = self._get_frame_atoms_info(idx)
        return self._extract_vector_components(frame_atoms_info, 'vx', 'vy', 'vz')

    def get_forces(self, idx: int) -> np.ndarray | None:
        frame_atoms_info = self._get_frame_atoms_info(idx)
        return self._extract_vector_components(frame_atoms_info, 'fx', 'fy', 'fz')

    def _get_cell_component(
        self, idx: int, component: int
    ) -> np.ndarray | list[bool] | None:
        """
        Helper to extract components from cell data.
        """
        pbc_cell = self.get('pbc_cell', [])
        try:
            return pbc_cell[idx][component]
        except (IndexError, TypeError):
            return None

    def get_lattice_vectors(self, idx: int) -> np.ndarray | None:
        _LATTICE_VECTORS_IDX = 1
        return self._get_cell_component(idx, _LATTICE_VECTORS_IDX)

    def get_pbc(self, idx: int) -> list[bool] | None:
        _PBC_FLAGS_IDX = 0
        return self._get_cell_component(idx, _PBC_FLAGS_IDX)

    def _get_trajectory_info(self, idx: int, key: str) -> int | None:
        """Helper to get specific trajectory info."""
        info = self.get(key)
        if info is None or idx >= len(info):
            return None
        return info[idx]

    def get_n_atoms(self, idx: int) -> int | None:
        n_atoms = self._get_trajectory_info(idx, 'n_atoms')
        if n_atoms is None:
            positions = self.get_positions(idx)
            return len(positions) if positions is not None else None
        return n_atoms

    def get_step(self, idx: int) -> int | None:
        return self._get_trajectory_info(idx, 'time_step')


class XYZTrajParser(TrajParser):
    """
    Parser for XYZ trajectory files.
    """

    def init_quantities(self) -> None:
        def get_atoms_info(val_in: str) -> dict[str, int | float]:
            val = [v.split('#')[0].split() for v in val_in.strip().splitlines()]
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
            return {key: val[n] for n, key in enumerate(['type', 'x', 'y', 'z'])}

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
    """
    Container for multiple trajectory parsers.

    Provides a unified interface to access trajectory data from multiple parser
    instances, with automatic fallback between parsers.
    """

    def __init__(
        self, parsers: list[TrajParser | XYZTrajParser | MDAnalysisParser]
    ) -> None:
        self._parsers = parsers
        self.logger = parsers[0].logger if parsers else get_logger(__name__)
        for parser in parsers:
            parser.parse()

    def __getitem__(
        self, index: int
    ) -> TrajParser | XYZTrajParser | MDAnalysisParser | None:
        try:
            return self._parsers[index]
        except (IndexError, TypeError):
            self.logger.warning('No parsers available or invalid index: %d', index)
            return None

    def eval(self, key: str, *args, **kwargs) -> Any | None:
        """
        Evaluate a method or property across all parsers.

        Returns:
            First non-None result from the parsers, or None if all return None
        """
        found_attribute = False
        val = None

        for parser in self._parsers:
            try:
                parser_attr = getattr(parser, key, None)
                if parser_attr is None:
                    continue

                found_attribute = True

                if callable(parser_attr):
                    # Always call if it's callable (even with no args)
                    val = parser_attr(*args, **kwargs)
                else:
                    # It's a property/attribute - only return if no args expected
                    if args or kwargs:
                        self.logger.warning(
                            'Arguments provided for non-callable attribute'
                        )
                        continue
                    val = parser_attr

            except Exception as e:
                self.logger.debug(
                    'Error evaluating input.',
                    exc_info=e,
                )
                continue

        if not found_attribute:
            self.logger.warning("Attribute '%s' not found in parsers", key)

        return val
