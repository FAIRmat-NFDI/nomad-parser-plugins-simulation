import os
import re
from difflib import SequenceMatcher
from typing import Any

import numpy as np
import pint
from nomad.units import ureg
from nomad_file_parser import Quantity, TextParser

from nomad_simulation_parsers.parsers.utils.constants import (
    MOLE,
    RE_FLOAT,
    RE_N,
)
from nomad_simulation_parsers.parsers.utils.general import search_files


def get_unit(
    units_type: str, property_type: Any | None = None, dimension: int = 3
) -> dict[str, Any]:
    units_type = units_type.lower()
    if units_type == 'real':
        units = dict(
            mass=ureg.g / MOLE,
            distance=ureg.angstrom,
            time=ureg.fs,
            energy=ureg.J * 4184.0 / MOLE,
            velocity=ureg.angstrom / ureg.fs,
            force=ureg.J * 4184.0 / ureg.angstrom / MOLE,
            torque=ureg.J * 4184.0 / MOLE,
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
            mass=ureg.g / MOLE,
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

    # TODO: Add support for lj systems
    # elif units_type == 'lj':
    #    units = dict(
    #        mass=ureg('dimensionless'),
    #        distance=ureg('dimensionless'),
    #        time=ureg('dimensionless'),
    #        energy=ureg('dimensionless'),
    #        velocity=ureg('dimensionless'),
    #        force=ureg('dimensionless'),
    #        torque=ureg('dimensionless'),
    #        temperature=ureg('dimensionless'),
    #        pressure=ureg('dimensionless'),
    #        dynamic_viscosity=ureg('dimensionless'),
    #        charge=ureg('dimensionless'),
    #        dipole=ureg('dimensionless'),
    #        electric_field=ureg('dimensionless'),
    #        density=ureg('dimensionless'),
    #    )
    #
    # ! Temporary, untested fix
    # LJ units according to the LAMMPS documentation
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


class DataParser(TextParser):
    def __init__(self) -> None:
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

    def init_quantities(self) -> None:
        self._quantities = [
            Quantity(header, rf'{RE_N} *(\d+) +{header}', repeats=True, dtype=np.int32)
            for header in self._headers
        ]

        def get_section_value(val: str) -> tuple[str | None, np.ndarray]:
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
                    return name, None

            return name, np.array(value)

        self._quantities.extend(
            [
                Quantity(
                    section,
                    rf'{section} *(#*.*{RE_N}\s+(?:[\d ]+{RE_FLOAT}.+\s+)+)',
                    str_operation=get_section_value,
                    repeats=True,
                )
                for section in self._sections
            ]
        )
        self._quantities.append(
            Quantity(
                'box_bounds',
                (
                    rf'({RE_FLOAT})\s+({RE_FLOAT})\s+xlo\s+xhi\s+'
                    rf'({RE_FLOAT})\s+({RE_FLOAT})\s+ylo\s+yhi\s+'
                    rf'({RE_FLOAT})\s+({RE_FLOAT})\s+zlo\s+zhi'
                ),
                dtype=float,
            )
        )
        self._quantities.append(
            Quantity(
                'tilt_factors',
                rf'({RE_FLOAT})\s+({RE_FLOAT})\s+({RE_FLOAT})\s+xy\s+xz\s+yz',
                dtype=float,
            )
        )

    def get_lattice_vectors(self) -> np.ndarray | None:
        box = self.get('box_bounds')
        if box is None:
            return None
        xlo, xhi, ylo, yhi, zlo, zhi = box
        lx, ly, lz = xhi - xlo, yhi - ylo, zhi - zlo
        xy, xz, yz = 0.0, 0.0, 0.0
        tilt = self.get('tilt_factors')
        if tilt is not None:
            xy, xz, yz = tilt
        matrix = np.array([[lx, 0.0, 0.0], [xy, ly, 0.0], [xz, yz, lz]])
        return matrix * ureg.angstrom

    def get_interactions(self) -> list[list]:
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
    # TODO: add parsers for other LAMMPS-supported trajectory formats
    _supported_traj_extensions = ['trj', 'xyz', 'dcd', 'lammpstrj']  # , 'h5', 'nc'
    _data_file_patterns = ['data', 'dat']

    def __init__(self) -> None:
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
        self._file_scan_warning = '{} file not specified in directory, will scan.'
        super().__init__(None)

    def init_quantities(self) -> None:
        def str_op(val: str) -> str | list[str]:
            val = val.split('#', maxsplit=1)[0]
            val = re.sub(f'&{RE_N}+', ' ', val)
            val = val.split()
            val = val if len(val) > 1 else val[0]
            return val

        self._quantities = [
            Quantity(
                name,
                (
                    rf'\n\s*{name}\s+'  # Name with whitespace
                    r'(?!.*\$\{)'  # No variable substitution
                    r'([${}\w\. \/\#\-]+)'  # Command arguments
                    r'(\&\n[\w\. \/\#\-]*)*'  # Line continuation with &
                ),
                # TODO: test!
                # LB - Edited regex (added \b \b) - word boundaries
                # (
                #     rf'\n\s*\b{name}\b\s+'  # Name with whitespace and word boundaries
                #     r'(?!.*\$\{)'  # No variable substitution
                #     r'([${}\w\. \/\#\-]+)'  # Command arguments
                #     r'(\&\n[\w\. \/\#\-]*)*'  # Line continuation with &
                # ),
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
                # TODO: test!
                # LB edit - Edited regex for '(2 Aug 2023 - Update 1)' searches for any
                # character except ')' now, not just word chars
                # r'\s*LAMMPS\s*\(([^)]+)\)\n',
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

        def str_to_thermo(val: str) -> dict[str, float]:
            res = {}
            if val.count('Step') > 1:
                # TODO: Test to make sure the regex substitution works as intended
                val = re.sub(r'--|=|\(sec\)', '', val).split()
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
                # r'([ \-]*Step\s*[\-/\s\w\.\=\(\)]*[ \-\.\d\n]+)Loop',  # TODO: LB edit
                str_operation=str_to_thermo,
                repeats=False,
                convert=False,
            )
        )

    @property
    def units(self) -> dict[str, pint.Quantity]:
        if self._units is None:
            units_type = self.get('units', ['lj'])[0]
            self._units = get_unit(units_type)
        return self._units

    def reset(self):
        super().reset()
        self._units = None

    def get_thermodynamic_data(self) -> dict[str, float] | None:
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

    # TODO: move to utils
    def find_best_matching_file(
        self, traj_files: list[str], mainfile_basename: str
    ) -> list[str]:
        """Find the best matching file based on the main input file name.
        Fall back to the first file if no good match is found."""

        prefix = os.path.basename(mainfile_basename).rsplit('.', 1)[1]
        matching_files = None

        # Strategy 1: Exact prefix match
        matching_files = [f for f in traj_files if prefix in f]

        if not matching_files:
            traj_data = [
                (os.path.basename(f).rsplit('.', 1)[0].lower(), f) for f in traj_files
            ]

            def get_tokens(filename: str) -> set[str]:
                """Extract tokens from filename by splitting on delimiters."""
                basename = os.path.basename(filename).rsplit('.', 1)[0].lower()
                tokens = re.split(r'[_\-\.\d]+', basename)
                MIN_CHARS = 2
                return {t for t in tokens if len(t) > MIN_CHARS}

            # TODO: check difflib or bisect to get best match
            def get_best_from_scores(
                scores: list[tuple[float, str]], min_score: float
            ) -> list[str] | None:
                """Return best match if score exceeds threshold, else None."""
                if not scores:
                    return None
                scores.sort(reverse=True, key=lambda x: x[0])
                best_score, best_file = scores[0]
                if best_score > min_score:
                    return [best_file]
                return None

            # Strategy 2: Token-based similarity
            mainfile_tokens = get_tokens(prefix)
            token_scores = []
            for traj_basename, traj_file in traj_data:
                traj_tokens = get_tokens(traj_basename)
                if mainfile_tokens or traj_tokens:
                    intersection = len(mainfile_tokens & traj_tokens)
                    union = len(mainfile_tokens | traj_tokens)
                    score = intersection / union if union > 0 else 0
                    token_scores.append((score, traj_file))

            matching_files = get_best_from_scores(token_scores, 0.0)

            # Strategy 3: String sequence similarity
            if not matching_files:
                sequence_scores = [
                    (
                        SequenceMatcher(None, prefix.lower(), traj_basename).ratio(),
                        traj_file,
                    )
                    for traj_basename, traj_file in traj_data
                ]
                matching_files = get_best_from_scores(sequence_scores, 0.3)

            # Fallback: Use first file
            if not matching_files:
                self.logger.warning(
                    (
                        'No match found for "%(mainfile)s". '
                        'Using first file: %(first_file)s'
                    ),
                    extra={
                        'mainfile': mainfile_basename,
                        'first_file': os.path.basename(traj_files[0]),
                    },
                )
                matching_files = [traj_files[0]]

        return matching_files

    def get_traj_files(self) -> list[str]:
        dump = self.get('dump', None)
        if dump is None:
            self.logger.warning(self._file_scan_warning.format('Trajectory'))
            traj_files = []
            for ext in self._supported_traj_extensions:
                found_files = search_files(
                    pattern=f'*.{ext}',
                    basedir=self.maindir,
                    deep=False,  # Only search current directory
                )
                traj_files.extend(found_files)
            # further eliminate
            if len(traj_files) > 1:
                traj_files = self.find_best_matching_file(traj_files, self.mainfile)
        else:
            traj_files = []
            if type(dump[0]) in [str, int]:
                dump = [dump]
            traj_files = [d[4] for d in dump]
        traj_files = list(dict.fromkeys(traj_files))  # remove duplicates

        return [os.path.join(self.maindir, f) for f in traj_files]

    def get_data_files(self) -> list[str]:
        """Get the data files either from the input script or by searching the main
        directory for files with LAMMPS data file extensions or header patterns."""

        def check_file_header(file_path: str, regex_pattern: str) -> None:
            header_size = 1024
            if not os.path.isabs(file_path):
                file_path = os.path.join(self.maindir, file_path)
            try:
                with open(file_path, 'rb') as file:
                    file_header = file.read(header_size)
                    file_header_str = file_header.decode(errors='ignore')
            except Exception:
                file_header_str = ''
            return re.search(regex_pattern, file_header_str)

        read_data = self.get('read_data')
        # TODO: Check in with Leo on the reasoning behind this change
        # # Chop out 'CPU' before, then just check none
        # if read_data is not None:
        #     try:
        #         read_data.remove('CPU')
        #     except Exception:
        #         pass
        if read_data is None or 'CPU' in read_data:
            self.logger.warning(self._file_scan_warning.format('Data'))
            for ext in self._data_file_patterns:
                data_files = search_files(
                    pattern=f'*{ext}*',
                    basedir=self.maindir,
                    deep=False,  # Only search current directory
                )
            if not data_files:
                # Search any file for the LAMMPS data file header.
                # Fallback to the LAMMPS input structure, if no run data file is found.
                data_files = os.listdir(self.maindir)
                patterns = ['LAMMPS data file', 'LAMMPS Description']
                for pattern in patterns:
                    data_files = [
                        f for f in data_files if check_file_header(f, pattern)
                    ]
                    if data_files:
                        break
            if len(data_files) > 1:
                # Search data files for the LAMMPS data file header.
                data_files = [
                    f for f in data_files if check_file_header(f, 'LAMMPS data file')
                ]
        else:
            data_files = read_data

        if not data_files:
            self.logger.warning('No data_files found to match the log file.')

        return [os.path.join(self.maindir, f) for f in data_files]

    def get_pbc(self) -> list[str]:
        pbc = self.get('boundary', ['p', 'p', 'p'])
        return [v == 'p' for v in pbc]

    def get_sampling_method(self) -> tuple[str, str]:
        fix_style = self.get('fix', [[''] * 3])[0][2]

        sampling_method = (
            'langevin_dynamics' if 'langevin' in fix_style else 'molecular_dynamics'
        )
        return sampling_method, fix_style

    def get_thermostat_settings(self) -> dict:
        fix = self.get('fix', [[]])[0]
        # TODO: check which items in <fix> need conversion to float
        fix = [float(x) if re.fullmatch(RE_FLOAT, str(x)) else x for x in fix]

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
                res['target_T'] = fix[5] * temp_unit
                res['thermostat_tau'] = fix[6] * time_unit
            except Exception:
                pass

        elif fix_style.lower() == 'npt':
            try:
                res['target_T'] = fix[5] * temp_unit
                res['thermostat_tau'] = fix[6] * time_unit
                res['target_P'] = fix[9] * press_unit
                res['barostat_tau'] = fix[10] * time_unit
            except Exception:
                pass

        elif fix_style.lower() == 'nph':
            try:
                res['target_P'] = fix[5] * press_unit
                res['barostat_tau'] = fix[6] * time_unit
            except Exception:
                pass

        elif fix_style.lower() == 'langevin':
            try:
                res['target_T'] = fix[4] * temp_unit
                res['langevin_gamma'] = fix[5] * time_unit
            except Exception:
                pass

        else:
            self.logger.warning('Fix style not supported', data=dict(style=fix_style))

        return res

    def get_interactions(self) -> list[tuple[str, list]]:
        styles_coeffs = []
        for interaction in self._interactions:
            styles = self.get(f'{interaction}_style', None)
            if styles is None:
                continue

            if isinstance(styles[0], str):
                styles = [styles]

            for i in range(len(styles)):
                if interaction == 'kspace':
                    coeff = [[float(c) for c in styles[i][1:]]]
                    style = styles[i][0]

                else:
                    coeff = self.get(f'{interaction}_coeff')
                    style = ' '.join([str(si) for si in styles[i]])

                styles_coeffs.append((style.strip(), coeff))

        return styles_coeffs
