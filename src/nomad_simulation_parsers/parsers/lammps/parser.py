import os
import re
import sys
from collections.abc import Iterable

import numpy as np
from ase import data as asedata
from nomad.datamodel import EntryArchive
from nomad.parsing.file_parser import ArchiveWriter, Quantity, TextParser
from nomad.parsing.parser import MatchingParser

# ? Do we want to migrate the MDAnalysisParser to nomad-parser-plugins-simulation, or is there a better solution?
from atomisticparsers.utils import MDAnalysisParser, MDParser
from nomad.units import ureg
from nomad_simulations.schema_packages.general import Program, Simulation
from structlog.stdlib import BoundLogger

re_float = r'[-+]?\d+\.*\d*(?:[Ee][-+]\d+)?'
re_n = r'[\n\r]'


def get_unit(units_type, property_type=None, dimension=3):
    mole = 6.022140857e23

    units_type = units_type.lower()
    if units_type == 'real':
        units = dict(
            mass=ureg.g / mole,
            distance=ureg.angstrom,
            time=ureg.fs,
            energy=ureg.J * 4184.0 / mole,
            velocity=ureg.angstrom / ureg.fs,
            force=ureg.J * 4184.0 / ureg.angstrom / mole,
            torque=ureg.J * 4184.0 / mole,
            temperature=ureg.K,
            pressure=ureg.atm,
            dynamic_viscosity=ureg.poise,
            charge=ureg.elementary_charge,
            dipole=ureg.elementary_charge * ureg.angstrom,
            electric_field=ureg.V / ureg.angstrom,
            density=ureg.g / ureg.cm**dimension,
        )

    elif units_type == 'metal':
        units = dict(
            mass=ureg.g / mole,
            distance=ureg.angstrom,
            time=ureg.ps,
            energy=ureg.eV,
            velocity=ureg.angstrom / ureg.ps,
            force=ureg.eV / ureg.angstrom,
            torque=ureg.eV,
            temperature=ureg.K,
            pressure=ureg.bar,
            dynamic_viscosity=ureg.poise,
            charge=ureg.elementary_charge,
            dipole=ureg.elementary_charge * ureg.angstrom,
            electric_field=ureg.V / ureg.angstrom,
            density=ureg.g / ureg.cm**dimension,
        )

    elif units_type == 'si':
        units = dict(
            mass=ureg.kg,
            distance=ureg.m,
            time=ureg.s,
            energy=ureg.J,
            velocity=ureg.m / ureg.s,
            force=ureg.N,
            torque=ureg.N * ureg.m,
            temperature=ureg.K,
            pressure=ureg.Pa,
            dynamic_viscosity=ureg.Pa * ureg.s,
            charge=ureg.C,
            dipole=ureg.C * ureg.m,
            electric_field=ureg.V / ureg.m,
            density=ureg.kg / ureg.m**dimension,
        )

    elif units_type == 'cgs':
        units = dict(
            mass=ureg.g,
            distance=ureg.cm,
            time=ureg.s,
            energy=ureg.erg,
            velocity=ureg.cm / ureg.s,
            force=ureg.dyne,
            torque=ureg.dyne * ureg.cm,
            temperature=ureg.K,
            pressure=ureg.dyne / ureg.cm**2,
            dynamic_viscosity=ureg.poise,
            charge=ureg.esu,
            dipole=ureg.esu * ureg.cm,
            electric_field=ureg.dyne / ureg.esu,
            density=ureg.g / ureg.cm**dimension,
        )

    elif units_type == 'electron':
        units = dict(
            mass=ureg.amu,
            distance=ureg.bohr,
            time=ureg.fs,
            energy=ureg.hartree,
            velocity=ureg.bohr / ureg.atomic_unit_of_time,
            force=ureg.hartree / ureg.bohr,
            temperature=ureg.K,
            pressure=ureg.Pa,
            charge=ureg.elementary_charge,
            dipole=ureg.debye,
            electric_field=ureg.V / ureg.cm,
        )

    elif units_type == 'micro':
        units = dict(
            mass=ureg.pg,
            distance=ureg.microm,
            time=ureg.micros,
            energy=ureg.pg * ureg.microm**2 / ureg.micros**2,
            velocity=ureg.microm / ureg.micros,
            force=ureg.pg * ureg.microm / ureg.micros**2,
            torque=ureg.pg * ureg.microm**2 / ureg.micros**2,
            temperature=ureg.K,
            pressure=ureg.pg / (ureg.microm * ureg.micros**2),
            dynamic_viscosity=ureg.pg / (ureg.microm * ureg.micros),
            charge=ureg.pC,
            dipole=ureg.pC * ureg.microm,
            electric_field=ureg.V / ureg.microm,
            density=ureg.pg / ureg.microm**dimension,
        )

    elif units_type == 'nano':
        units = dict(
            mass=ureg.ag,
            distance=ureg.nm,
            time=ureg.ns,
            energy=ureg.ag * ureg.nm**2 / ureg.ns**2,
            velocity=ureg.nm / ureg.ns,
            force=ureg.ag * ureg.nm / ureg.ns**2,
            torque=ureg.ag * ureg.nm**2 / ureg.ns**2,
            temperature=ureg.K,
            pressure=ureg.ag / (ureg.nm * ureg.ns**2),
            dynamic_viscosity=ureg.ag / (ureg.nm * ureg.ns),
            charge=ureg.elementary_charge,
            dipole=ureg.elementary_charge * ureg.nm,
            electric_field=ureg.V / ureg.nm,
            density=ureg.ag / ureg.nm**dimension,
        )

    elif units_type == 'lj':
        units = dict(
            mass=1,
            distance=1,
            time=1,
            energy=1,
            velocity=1,
            force=1,
            torque=1,
            temperature=1,
            pressure=1,
            dynamic_viscosity=1,
            charge=1,
            dipole=1,
            electric_field=1,
            density=1,
        )
    else:
        units = dict()

    if property_type:
        return units.get(property_type, None)
    else:
        return units


class TrajParser(TextParser):
    def __init__(self):
        self._masses = None
        self._reference_masses = dict(
            masses=np.array(asedata.atomic_masses), symbols=asedata.chemical_symbols
        )
        self._chemical_symbols = None
        super().__init__(None)

    def init_quantities(self):
        def get_pbc_cell(val):
            val = val.split()

            pbc = [v == 'pp' for v in val[:3]]

            cell = np.zeros((3, 3))
            for i in range(3):
                cell[i][i] = float(val[i * 2 + 4]) - float(val[i * 2 + 3])

            return pbc, cell

        def get_atoms_info(val):
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
                r'\s*ITEM: BOX BOUNDS\s*([\s\w]+)\n([\+\-\d\.eE\s]+)\n',
                str_operation=get_pbc_cell,
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
    def with_trajectory(self):
        return self.get('atoms_info') is not None

    @property
    def n_frames(self):
        return len(self.get('atoms_info', []))

    @property
    def masses(self):
        return self._masses

    @masses.setter
    def masses(self, val):
        self._masses = val
        if self._masses is None:
            return

        self._masses = val
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

    def get_atom_labels(self, idx):
        atoms_info = self.get('atoms_info')
        if atoms_info is None:
            return

        atoms_id = atoms_info[idx].get('id')
        default = ['X' for _ in atoms_id] if atoms_id is not None else None
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

    def get_positions(self, idx):
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

    def get_velocities(self, idx):
        atoms_info = self.get('atoms_info')
        if atoms_info is None:
            return
        atoms_info = atoms_info[idx]
        if 'vx' not in atoms_info or 'vy' not in atoms_info or 'vz' not in atoms_info:
            return

        return np.array([atoms_info['vx'], atoms_info['vy'], atoms_info['vz']]).T

    def get_forces(self, idx):
        atoms_info = self.get('atoms_info')
        if atoms_info is None:
            return
        atoms_info = atoms_info[idx]
        if 'fx' not in atoms_info or 'fy' not in atoms_info or 'fz' not in atoms_info:
            return
        return np.array([atoms_info['fx'], atoms_info['fy'], atoms_info['fz']]).T

    def get_lattice_vectors(self, idx):
        pbc_cell = self.get('pbc_cell')
        if pbc_cell is None:
            return
        return pbc_cell[idx][1]

    def get_pbc(self, idx):
        pbc_cell = self.get('pbc_cell')
        if pbc_cell is None:
            return
        return pbc_cell[idx][0]

    def get_n_atoms(self, idx):
        n_atoms = self.get('n_atoms')
        if n_atoms is None:
            return len(self.get_positions(idx))
        return n_atoms[idx]

    def get_step(self, idx):
        step = self.get('time_step')
        if step is None:
            return
        return step[idx]


class XYZTrajParser(TrajParser):
    def __init__(self):
        super().__init__()

    def init_quantities(self):
        def get_atoms_info(val_in):
            val = [v.split('#')[0].split() for v in val_in.strip().split('\n')]
            symbols = []
            for v in val:
                if v[0].isalpha():
                    if v[0] not in symbols:
                        symbols.append(v[0])
                    v[0] = symbols.index(v[0]) + 1
            val = np.transpose(np.array([v for v in val if len(v) == 4], dtype=float))
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
    def __init__(self, parsers):
        self._parsers = parsers
        for parser in parsers:
            parser.parse()

    def __getitem__(self, index):
        if self._parsers:
            return self._parsers[index]

    def eval(self, key, *args, **kwargs):
        for parser in self._parsers:
            parser_method = getattr(parser, key)
            if parser_method is not None:
                val = (
                    parser_method(*args, **kwargs) if args or kwargs else parser_method
                )
                if val is not None:
                    return val


class DataParser(TextParser):
    def __init__(self):
        self._headers = [
            'atoms',
            'bonds',
            'angles',
            'dihedrals',
            'impropers',
            'atom types',
            'bond types',
            'angle types',
            'dihedral types',
            'improper types',
            'extra bond per atom',
            'extra/bond/per/atom',
            'extra angle per atom',
            'extra/angle/per/atom',
            'extra dihedral per atom',
            'extra/dihedral/per/atom',
            'extra improper per atom',
            'extra/improper/per/atom',
            'extra special per atom',
            'extra/special/per/atom',
            'ellipsoids',
            'lines',
            'triangles',
            'bodies',
        ]
        self._sections = [
            'Atoms',
            'Velocities',
            'Masses',
            'Ellipsoids',
            'Lines',
            'Triangles',
            'Bodies',
            'Bonds',
            'Angles',
            'Dihedrals',
            'Impropers',
            'Pair Coeffs',
            'PairIJ Coeffs',
            'Bond Coeffs',
            'Angle Coeffs',
            'Dihedral Coeffs',
            'Improper Coeffs',
            'BondBond Coeffs',
            'BondAngle Coeffs',
            'MiddleBondTorsion Coeffs',
            'EndBondTorsion Coeffs',
            'AngleTorsion Coeffs',
            'AngleAngleTorsion Coeffs',
            'BondBond13 Coeffs',
            'AngleAngle Coeffs',
        ]
        self._interactions = [
            section for section in self._sections if section.endswith('Coeffs')
        ]
        super().__init__(None)

    def init_quantities(self):
        self._quantities = [
            Quantity(header, rf'{re_n} *(\d+) +{header}', repeats=True, dtype=np.int32)
            for header in self._headers
        ]

        def get_section_value(val):
            val = val.strip().splitlines()
            name = None

            if val[0][0] == '#':
                name = val[0][1:].strip()
                val = val[1:]

            value = []
            for i in range(len(val)):
                v = val[i].split('#')[0].split()
                if not v:
                    continue

                try:
                    value.append(np.array(v, dtype=float))
                except Exception:
                    break

            return name, np.array(value)

        self._quantities.extend(
            [
                Quantity(
                    section,
                    rf'{section} *(#*.*{re_n}\s+(?:[\d ]+{re_float}.+\s+)+)',
                    str_operation=get_section_value,
                    repeats=True,
                )
                for section in self._sections
            ]
        )

    def get_interactions(self):
        styles_coeffs = []
        for interaction in self._interactions:
            coeffs = self.get(interaction, None)
            if coeffs is None:
                continue
            if isinstance(coeffs, tuple):
                coeffs = list(coeffs)

            styles_coeffs += coeffs

        return styles_coeffs


class LogParser(TextParser):
    def __init__(self):
        self._commands = [
            'angle_coeff',
            'angle_style',
            'atom_modify',
            'atom_style',
            'balance',
            'bond_coeff',
            'bond_style',
            'bond_write',
            'boundary',
            'change_box',
            'clear',
            'comm_modify',
            'comm_style',
            'compute',
            'compute_modify',
            'create_atoms',
            'create_bonds',
            'create_box',
            'delete_bonds',
            'dielectric',
            'dihedral_coeff',
            'dihedral_style',
            'dimension',
            'displace_atoms',
            'dump',
            'dump_modify',
            'dynamical_matrix',
            'echo',
            'fix',
            'fix_modify',
            'group',
            'group2ndx',
            'ndx2group',
            'hyper',
            'if',
            'improper_coeff',
            'improper_style',
            'include',
            'info',
            'jump',
            'kim_init',
            'kim_interactions',
            'kim_query',
            'kim_param',
            'kim_property',
            'kspace_modify',
            'kspace_style',
            'label',
            'lattice',
            'log',
            'mass',
            'message',
            'min_modify',
            'min_style',
            'minimize',
            'minimize/kk',
            'molecule',
            'neb',
            'neb/spin',
            'neigh_modify',
            'neighbor',
            'newton',
            'next',
            'package',
            'pair_coeff',
            'pair_modify',
            'pair_style',
            'pair_write',
            'partition',
            'prd',
            'print',
            'processors',
            'quit',
            'read_data',
            'read_dump',
            'read_restart',
            'region',
            'replicate',
            'rerun',
            'reset_atom_ids',
            'reset_mol_ids',
            'reset_timestep',
            'restart',
            'run',
            'run_style',
            'server',
            'set',
            'shell',
            'special_bonds',
            'suffix',
            'tad',
            'temper/grem',
            'temper/npt',
            'thermo',
            'thermo_modify',
            'thermo_style',
            'third_order',
            'timer',
            'timestep',
            'uncompute',
            'undump',
            'unfix',
            'units',
            'variable',
            'velocity',
            'write_coeff',
            'write_data',
            'write_dump',
            'write_restart',
        ]
        self._interactions = [
            'atom',
            'pair',
            'bond',
            'angle',
            'dihedral',
            'improper',
            'kspace',
        ]
        self._units = None
        super().__init__(None)

    def init_quantities(self):
        def str_op(val):
            val = val.split('#')[0]
            val = val.replace('&\n', ' ').split()
            val = val if len(val) > 1 else val[0]
            return val

        self._quantities = [
            Quantity(
                name,
                r'\n\s*%s\s+(?!.*\$\{)([${}\w\. \/\#\-]+)(\&\n[\w\. \/\#\-]*)*' % name,
                str_operation=str_op,
                comment='#',
                repeats=True,
            )
            for name in self._commands
        ]

        self._quantities.append(
            Quantity(
                'program_version',
                r'\s*LAMMPS\s*\(([\w ]+)\)\n',
                dtype=str,
                repeats=False,
                flatten=False,
            )
        )

        self._quantities.append(
            Quantity('finished', r'\s*Dangerous builds\s*=\s*(\d+)', repeats=False)
        )

        self._quantities.append(
            Quantity(
                'minimization_stats',
                r'\s*Minimization stats:\s*([\s\S]+?)\n\n',
                flatten=False,
            )
        )

        def str_to_thermo(val):
            res = {}
            if val.count('Step') > 1:
                val = (
                    val.replace('--', '').replace('=', '').replace('(sec)', '').split()
                )
                val = [v.strip() for v in val]

                for i in range(len(val)):
                    if val[i][0].isalpha():
                        res.setdefault(val[i], [])
                        res[val[i]].append(float(val[i + 1]))

            else:
                val = val.split('\n')
                keys = [v.strip() for v in val[0].split()]
                val = np.array([v.split() for v in val[1:] if v], dtype=float).T

                res = {key: [] for key in keys}
                for i in range(len(keys)):
                    res[keys[i]] = val[i]

            return res

        self._quantities.append(
            Quantity(
                'thermo_data',
                r'\s*\-*(\s*Step\s*[\-\s\w\.\=\(\)]*[ \-\.\d\n]+)Loop',
                str_operation=str_to_thermo,
                repeats=False,
                convert=False,
            )
        )

    @property
    def units(self):
        if self._units is None:
            units_type = self.get('units', ['lj'])[0]
            self._units = get_unit(units_type)
        return self._units

    def get_thermodynamic_data(self):
        thermo_data = self.get('thermo_data')

        if thermo_data is None:
            return

        data = {}
        for key, val in thermo_data.items():
            low_key = key.lower()
            if low_key.startswith('e_') or low_key.endswith('eng'):
                data[key] = val * self.units.get('energy', 1)
            elif low_key == 'press':
                data[key] = val * self.units.get('pressure', 1)
            elif low_key == 'temp':
                data[key] = val * self.units.get('temperature', 1)
            else:
                data[key] = val
        return data

    def get_traj_files(self):
        dump = self.get('dump')
        if dump is None:
            self.logger.warning('Trajectory not specified in directory, will scan.')
            # TODO improve matching of traj file
            traj_files = os.listdir(self.maindir)
            traj_files = [
                f for f in traj_files if f.endswith('trj') or f.endswith('xyz')
            ]
            # further eliminate
            if len(traj_files) > 1:
                prefix = os.path.basename(self.mainfile).rsplit('.', 1)[0]
                traj_files = [f for f in traj_files if prefix in f]
        else:
            traj_files = []
            if type(dump[0]) in [str, int]:
                dump = [dump]
            traj_files = [d[4] for d in dump]
        traj_files = [
            i for n, i in enumerate(traj_files) if i not in traj_files[:n]
        ]  # remove duplicates

        return [os.path.join(self.maindir, f) for f in traj_files]

    def get_data_files(self):
        def check_file_header(file_path, regex_pattern):
            print('regex_pattern:', regex_pattern)
            header_size = 1024
            file_path = f'{self.maindir}/{file_path}'
            try:
                with open(file_path, 'rb') as file:
                    file_header = file.read(header_size)
                    file_header_str = file_header.decode(errors='ignore')
            except Exception:
                file_header_str = ''

            return re.search(regex_pattern, file_header_str)

        read_data = self.get('read_data')
        if read_data is None or 'CPU' in read_data:
            self.logger.warning('Data file not specified in directory, will scan.')
            data_files = os.listdir(self.maindir)
            data_files = [
                f for f in data_files if f.endswith('data') or f.startswith('data')
            ]
            if not data_files:
                data_files = os.listdir(self.maindir)
                data_files = [
                    f for f in data_files if check_file_header(f, 'LAMMPS data file')
                ]  # TODO: Should this be the default?
            if len(data_files) > 1:
                prefix = os.path.basename(self.mainfile).rsplit('.', 1)
                prefix = (
                    prefix[1] if len(prefix) > 1 and prefix[1] != 'log' else prefix[0]
                )
                data_files = [f for f in data_files if prefix in f]
        else:
            data_files = read_data

        return [os.path.join(self.maindir, f) for f in data_files]

    def get_pbc(self):
        pbc = self.get('boundary', ['p', 'p', 'p'])
        return [v == 'p' for v in pbc]

    def get_sampling_method(self):
        fix_style = self.get('fix', [[''] * 3])[0][2]

        sampling_method = (
            'langevin_dynamics' if 'langevin' in fix_style else 'molecular_dynamics'
        )
        return sampling_method, fix_style

    def get_thermostat_settings(self):
        fix = self.get('fix', [None])[0]
        if fix is None:
            return {}

        try:
            fix_style = fix[2]
        except IndexError:
            return {}

        temp_unit = self.units.get('temperature', 1)
        press_unit = self.units.get('pressure', 1)
        time_unit = self.units.get('time', 1)

        res = dict()
        if fix_style.lower() == 'nvt':
            try:
                res['target_T'] = float(fix[5]) * temp_unit
                res['thermostat_tau'] = float(fix[6]) * time_unit
            except Exception:
                pass

        elif fix_style.lower() == 'npt':
            try:
                res['target_T'] = float(fix[5]) * temp_unit
                res['thermostat_tau'] = float(fix[6]) * time_unit
                res['target_P'] = float(fix[9]) * press_unit
                res['barostat_tau'] = float(fix[10]) * time_unit
            except Exception:
                pass

        elif fix_style.lower() == 'nph':
            try:
                res['target_P'] = float(fix[5]) * press_unit
                res['barostat_tau'] = float(fix[6]) * time_unit
            except Exception:
                pass

        elif fix_style.lower() == 'langevin':
            try:
                res['target_T'] = float(fix[4]) * temp_unit
                res['langevin_gamma'] = float(fix[5]) * time_unit
            except Exception:
                pass

        else:
            self.logger.warning('Fix style not supported', data=dict(style=fix_style))

        return res

    def get_interactions(self):
        styles_coeffs = []
        for interaction in self._interactions:
            styles = self.get('%s_style' % interaction, None)
            if styles is None:
                continue

            if isinstance(styles[0], str):
                styles = [styles]

            for i in range(len(styles)):
                if interaction == 'kspace':
                    coeff = [[float(c) for c in styles[i][1:]]]
                    style = styles[i][0]

                else:
                    coeff = self.get('%s_coeff' % interaction)
                    style = ' '.join([str(si) for si in styles[i]])

                styles_coeffs.append((style.strip(), coeff))

        return styles_coeffs


class LammpsArchiveWriter(ArchiveWriter):
    """
    TODO: Docstring
    """

    def write_to_archive(self) -> None:
        self.archive.data = Simulation(program=Program(name='LAMMPS'))
        # LAMMPS mainfile is the main log file
        self.basename = os.path.basename(self.mainfile)
        self.basedir = os.path.dirname(self.mainfile)
        self.log_parser = LogParser()
        self.log_parser.mainfile = self.mainfile
        self.log_parser.logger = self.logger
        self.log_parser._units = None

        # parse data from auxiliary log file
        self.aux_log_parser = LogParser()
        if self.log_parser.get('log') is not None:
            self.aux_log_parser.mainfile = os.path.join(
                self.log_parser.maindir,
                self.log_parser.get('log')[0],
            )
            # we assign units here which is read from log parser
            self.aux_log_parser._units = self.log_parser.units
            self.aux_log_parser.logger = self.logger

        self._traj_parser = TrajParser()
        self._traj_parser.logger = self.logger
        self._traj_parser._chemical_symbols = None
        self._xyztraj_parser = XYZTrajParser()
        self._xyztraj_parser.logger = self.logger
        self._mdanalysistraj_parser = MDAnalysisParser(
            topology_format='DATA', format='LAMMPSDUMP'
        )
        self._mdanalysistraj_parser.logger = self.logger
        self.data_parser = DataParser()
        self.data_parser.logger = self.logger

        # parse data file associated with calculation
        data_files = self.log_parser.get_data_files()
        print(data_files)
        if len(data_files) > 1:
            self.logger.warning('Multiple data files are specified')
        if data_files:
            self.data_parser.mainfile = data_files[0]

        sys.exit()

        # TODO extend


class LammpsParser(MatchingParser):
    """
    Main parser interface to NOMAD.
    """

    archive_writer = LammpsArchiveWriter()

    # ? Really needed for the LAMMPS parser?
    # ? Would it make sense to handle a potential auxillary log file here,
    # ? since there seems to be an issue with its handling in the old parser?
    def is_mainfile(
        self,
        filename: str,
        mime: str,
        buffer: bytes,
        decoded_buffer: str,
        compression: str = None,
    ) -> bool | Iterable[str]:
        """
        TODO: Documentation
        """
        is_mainfile = super().is_mainfile(
            filename, mime, buffer, decoded_buffer, compression
        )

        if is_mainfile:
            # ? Handle check for auxiliary log file here?
            return is_mainfile

    def parse(
        self,
        mainfile: str,
        archive: EntryArchive,
        logger: BoundLogger = None,
        child_archives: dict[str, EntryArchive] = {},
    ) -> None:
        self.archive_writer.write(mainfile, archive, logger, child_archives)
