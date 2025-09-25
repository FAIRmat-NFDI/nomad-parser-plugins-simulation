from typing import Any

import numpy as np
from ase import data as asedata
from nomad.parsing.file_parser import Quantity, TextParser


class TrajParser(TextParser):
    def __init__(self) -> None:
        self._masses = None
        self._reference_masses = dict(
            masses=np.array(asedata.atomic_masses), symbols=asedata.chemical_symbols
        )
        self._chemical_symbols = None
        super().__init__(None)

    def get_pbc_cell(self, val) -> tuple[list, np.ndarray]:
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
        def get_atoms_info(val) -> dict[str, float]:
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
    def masses(self, val) -> None:
        self._masses = val
        if self._masses is None:
            return

        if self._chemical_symbols is None:
            masses = self._masses[0][1]
            self._chemical_symbols = {}
            for i in range(len(masses)):
                symbol_idx = np.argmin(
                    abs(self._reference_masses['masses'] - masses[i][1])
                )
                self._chemical_symbols[masses[i][0]] = self._reference_masses[
                    'symbols'
                ][symbol_idx]

    def get_atom_labels(self, idx) -> list[str] | None:
        atoms_info = self.get('atoms_info')
        if atoms_info is None:
            return

        atoms_id = atoms_info[idx].get('id')
        default = ['CGX' for _ in atoms_id] if atoms_id is not None else None
        atoms_type = atoms_info[idx].get('type')
        if atoms_type is None:
            return default
        if self._chemical_symbols is None:
            return default

        try:
            atom_labels = [self._chemical_symbols[atype] for atype in atoms_type]
        except Exception:
            self.logger.error('Error resolving atom labels.')
            return

        return atom_labels

    def get_positions(self, idx) -> np.ndarray | None:
        atoms_info = self.get('atoms_info')
        if atoms_info is None:
            return

        atoms_info = atoms_info[idx]

        cell = self.get('pbc_cell')
        cell = None if cell is None else cell[idx][1]
        if 'xs' in atoms_info and 'ys' in atoms_info and 'zs' in atoms_info:
            if cell is None:
                return
            positions = np.array(
                [atoms_info['xs'], atoms_info['ys'], atoms_info['zs']]
            ).T
            positions = positions * np.linalg.norm(cell, axis=1) + np.amin(cell, axis=1)

        elif 'xu' in atoms_info and 'yu' in atoms_info and 'zu' in atoms_info:
            positions = np.array(
                [atoms_info['xu'], atoms_info['yu'], atoms_info['zu']]
            ).T

        elif 'xsu' in atoms_info and 'ysu' in atoms_info and 'zsu' in atoms_info:
            if cell is None:
                return
            positions = np.array(
                [atoms_info['xsu'], atoms_info['ysu'], atoms_info['zsu']]
            ).T
            positions = positions * np.linalg.norm(cell, axis=1) + np.amin(cell, axis=1)

        elif 'x' in atoms_info and 'y' in atoms_info and 'z' in atoms_info:
            positions = np.array([atoms_info['x'], atoms_info['y'], atoms_info['z']]).T
            if 'ix' in atoms_info and 'iy' in atoms_info and 'iz' in atoms_info:
                if cell is None:
                    return
                positions_img = np.array(
                    [atoms_info['ix'], atoms_info['iy'], atoms_info['iz']]
                ).T

                positions += positions_img * np.linalg.norm(cell, axis=1)
        else:
            positions = None

        return positions

    def get_velocities(self, idx) -> np.ndarray | None:
        atoms_info = self.get('atoms_info')
        if atoms_info is None:
            return
        atoms_info = atoms_info[idx]
        if 'vx' not in atoms_info or 'vy' not in atoms_info or 'vz' not in atoms_info:
            return

        return np.array([atoms_info['vx'], atoms_info['vy'], atoms_info['vz']]).T

    def get_forces(self, idx) -> np.ndarray | None:
        atoms_info = self.get('atoms_info')
        if atoms_info is None:
            return
        atoms_info = atoms_info[idx]
        if 'fx' not in atoms_info or 'fy' not in atoms_info or 'fz' not in atoms_info:
            return
        return np.array([atoms_info['fx'], atoms_info['fy'], atoms_info['fz']]).T

    def get_lattice_vectors(self, idx) -> np.ndarray | None:
        pbc_cell = self.get('pbc_cell')
        if pbc_cell is None:
            return
        return pbc_cell[idx][1]

    def get_pbc(self, idx) -> list[bool] | None:
        pbc_cell = self.get('pbc_cell')
        if pbc_cell is None:
            return
        return pbc_cell[idx][0]

    def get_n_atoms(self, idx) -> int | None:
        n_atoms = self.get('n_atoms')
        if n_atoms is None:
            return len(self.get_positions(idx))
        return n_atoms[idx]

    def get_step(self, idx) -> int | None:
        step = self.get('time_step')
        if step is None:
            return
        return step[idx]


class XYZTrajParser(TrajParser):
    def __init__(self) -> None:
        super().__init__()

    def init_quantities(self) -> None:
        def get_atoms_info(val_in) -> dict[str, int | float]:
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
    def __init__(self, parsers) -> None:
        self._parsers = parsers
        for parser in parsers:
            parser.parse()

    # ? Is this function used anywhere?
    # ? Should the return types for the else case be handled better?
    def __getitem__(self, index) -> TrajParser | None:
        if self._parsers:
            return self._parsers[index]

    # ? Also here, should we make the negative return type more explicit?
    def eval(self, key, *args, **kwargs) -> Any | None:
        for parser in self._parsers:
            parser_method = getattr(parser, key)
            if parser_method is not None:
                val = (
                    parser_method(*args, **kwargs) if args or kwargs else parser_method
                )
                if val is not None:
                    return val
