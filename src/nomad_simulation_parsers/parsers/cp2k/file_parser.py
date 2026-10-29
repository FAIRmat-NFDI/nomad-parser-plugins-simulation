#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD.
# See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# ruff: noqa: E501, PLR0911, PLR0912, PLR2004, UP031
import re

import numpy as np
from nomad.units import ureg
from nomad_file_parser import Quantity, TextParser

units_map = {
    'hbar': ureg.hbar,
    'hartree': ureg.hartree,
    'angstrom': ureg.angstrom,
    'au_t': ureg.hbar / ureg.hartree,
    'fs': ureg.femtosecond,
    'K': ureg.kelvin,
}


def resolve_unit(unit_str, parts=[]):
    unit_str = unit_str.lower().replace(' ', '')
    parts = list(parts)

    if unit_str in units_map:
        return units_map[unit_str]

    try:
        return float(unit_str)
    except Exception:
        pass

    if unit_str == '':
        return 1

    open_p = unit_str.rfind('(')
    if open_p > -1:
        n_groups = unit_str.count('(')
        if n_groups != unit_str.count(')'):
            return unit_str
        for n in range(n_groups):
            part = unit_str[open_p + 1 :]
            part = part[: part.find(')')]
            parts.append(resolve_unit(part, parts))
            unit_str = unit_str.replace('(%s)' % part, '[%d]' % n)
            open_p = unit_str.rfind('(')
        return resolve_unit(unit_str, parts)

    vals = unit_str.split('/')
    if len(vals) > 1:
        vals = [resolve_unit(v, parts) for v in vals]
        val = vals[0]
        for v in vals[1:]:
            val /= v
        return val

    vals = unit_str.split('**')
    if len(vals) > 1:
        vals = [resolve_unit(v, parts) for v in vals]
        val = vals[0]
        for v in reversed(vals[1:]):
            val = val**v
        return val

    vals = unit_str.split('^')
    if len(vals) > 1:
        vals = [resolve_unit(v, parts) for v in vals]
        val = vals[0]
        for v in reversed(vals[1:]):
            val = val**v
        return val

    vals = unit_str.split('*')
    if len(vals) > 1:
        vals = [resolve_unit(v, parts) for v in vals]
        unit = 1
        for v in vals:
            unit *= v
        return unit

    vals = unit_str.split('-1')
    if len(vals) == 2:
        return 1 / resolve_unit(vals[0], parts)

    vals = re.match(r'\[(\d+)\]', unit_str)
    if vals:
        return parts[int(vals.group(1))]


class Property:
    def __init__(self, **kwargs):
        self._data = kwargs

    def __getattr__(self, key):
        return self._data.get(key, None)


class Trajectory(Property):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)


re_float = r'[-+]?\d+\.?\d*(?:[Ee][-+]\d+)?'


class CP2KOutParser(TextParser):
    def __init__(self):
        super().__init__()

    def init_quantities(self):
        def str_to_header(val_in):
            val = val_in.split('  ', 1)
            return [val[0].strip().replace(' ', '_').lower(), val[-1].strip()]

        def md_extract(val_in):
            result = re.search(
                r' MD\| (?P<key>.+?)(?: \[(?P<unit>.+)\])? {2,}(?P<value>.+)', val_in
            )
            value = result.group('value')
            unit = units_map.get(result.group('unit'))
            key = result.group('key').strip().replace(' ', '_').lower()
            if unit:
                value = float(value) * unit
            return [key, value]

        def str_to_program(val_in):
            val = val_in.split(' ', 2)
            return [
                '_'.join(val[:2]).lower(),
                ''.join([v.strip() for v in val[2].split('\n')]),
            ]

        def str_to_atomic_coordinates(val_in):
            val = [v.split() for v in val_in.split('\n')]
            length_unit = val[0][0].lower()
            val = np.transpose(np.array([v for v in val if len(v) == 9]))
            labels = val[2]
            positions = np.transpose(np.array(val[4:7], dtype=float)) * resolve_unit(
                length_unit
            )
            atomic_numbers = {
                element: int(val[3][n]) for n, element in enumerate(val[2])
            }
            return Trajectory(
                labels=labels, positions=positions, atomic_numbers=atomic_numbers
            )

        def str_to_stress_eigenvalues(val_in):
            val = [v.split() for v in val_in.strip().split('\n')]
            val = np.array([v for v in val if v], dtype=float)
            return val[0] * ureg.GPa, val[1:]

        def str_to_iteration(val_in):
            val = val_in.strip().split()
            return {
                'time_calculation': float(val[-4]) * ureg.s,
                'energy_total': float(val[-2]) * ureg.hartree,
                'energy_change': float(val[-1]) * ureg.hartree,
            }

        def str_to_information(val_in):
            val = [v.split('=') for v in val_in.strip().split('\n')]
            return {
                v[0].strip().lower().replace(' ', '_').replace('.', ''): v[1]
                for v in val
                if len(v) == 2
            }

        n_orbital_basis_quantities = [
            Quantity(
                'basis_set_number_of_%s' % key.lower().replace(' ', '_'),
                r'Number of %s:\s+(\d+)' % key,
                dtype=int,
            )
            for key in [
                'orbital shell sets',
                'orbital shells',
                'primitive Cartesian functions',
                'Cartesian basis functions',
                'spherical basis functions',
            ]
        ]

        energy_quantities = [
            Quantity(
                '%s' % key.lower().replace(' ', '_').replace('-', '_'),
                rf'%s:\s*({re_float})' % key,
                dtype=float,
                unit='hartree',
                repeats=True,
            )
            for key in [
                'Hartree energy',
                'Exchange-correlation energy',
                'Electronic kinetic energy',
                'Total energy',
                'Fermi energy',
            ]
        ]
        # what is the difference between Total energy and ENERGY| Total

        scf_wavefunction_optimization_quantities = [
            Quantity(
                'iteration',
                rf'(?m:^[ \t]*)(\d+[ \t]+\S+(?:[ \t]+\S+)?'
                rf'[ \t]+{re_float}[ \t]+{re_float}[ \t]+{re_float}'
                rf'[ \t]+{re_float}[ \t]+{re_float})[ \t]*\r?\n',
                repeats=True,
                convert=False,
                str_operation=str_to_iteration,
            ),
            # TODO add minimizer info
            Quantity('converged', r'SCF run converged in\s*(\d+) steps', dtype=int),
            # find example with cueb file
            Quantity(
                'cube_file',
                r' The electron density is written in cube file format to the file:\s*(.+?\.cube)',
            ),
            # TODO add mulliken, hirschfield
            Quantity(
                'energy_total',
                rf'ENERGY\| Total FORCE_EVAL \( \w+ \) energy \(a\.u\.\):\s*({re_float})',
                dtype=float,
                unit='hartree',
            ),
            Quantity(
                'atom_forces',
                r'ATOMIC FORCES in \[a\.u\.\]\s*.+([\s\S]+?)SUM',
                convert=False,
                str_operation=lambda x: (
                    np.array(
                        np.transpose([v.split() for v in x.strip().split('\n')])[3:6],
                        dtype=float,
                    ).T
                ),
            ),
            # TODO test stress cannot find example
            Quantity(
                'stress_tensor',
                r' (?:NUMERICAL )?STRESS TENSOR \[GPa\]\s+X\s+Y\s+Z\s+([\d\.\-\s]+)',
                str_operation=lambda x: np.array(
                    [v.split() for v in x.strip().split('\n')], dtype=float
                ),
                convert=False,
                unit='GPa',
            ),
            Quantity(
                'stress_tensor_one_third_of_trace',
                rf'  1/3 Trace\(stress tensor\)\s*:\s*({re_float})',
                dtype=float,
                unit='GPa',
            ),
            Quantity(
                'stress_tensor_determinant',
                rf'Det\(stress tensor\)\s*:\s*({re_float})',
                dtype=float,
                unit='GPa**3',
            ),
            Quantity(
                'stress_eigenvalues_eigenvectors',
                r' EIGENVECTORS AND EIGENVALUES OF THE STRESS TENSOR\s*([\d\.\-\s]+)',
                str_operation=str_to_stress_eigenvalues,
                convert=False,
            ),
        ] + energy_quantities

        geometry_optimization_quantities = [
            Quantity(
                'method',
                r'\*{3}\s*((?:CONJUGATE GRADIENTS|L\-BFGS|BFGS))\s*\*{3}',
                flatten=False,
            ),
            Quantity(
                'self_consistent',
                r'SCF WAVEFUNCTION OPTIMIZATION([\s\S]+?)OPTIMIZ',
                repeats=False,
                sub_parser=TextParser(
                    quantities=scf_wavefunction_optimization_quantities
                ),
            ),
            Quantity(
                'optimization_step',
                r'(ATION STEP:\s*\d+[\s\S]+?)(?:\-\s+OPTIMIZ|\Z)',
                repeats=True,
                sub_parser=TextParser(
                    quantities=[
                        # TODO parse atomic positions
                        Quantity('step', r'ATION STEP:\s*(\d+)'),
                        # I do not quite get why there can be multiple scfs in a step
                        Quantity(
                            'information',
                            r'Informations at step\s*=\s*\d+\s*\-+([\s\S]+?)\-{5}',
                            str_operation=str_to_information,
                        ),
                        Quantity(
                            'self_consistent',
                            r'FUNCTION OPTIMIZATION([\s\S]+?)(?: SCF WAVE|\Z)',
                            repeats=True,
                            sub_parser=TextParser(
                                quantities=scf_wavefunction_optimization_quantities
                            ),
                        ),
                    ]
                ),
            ),
        ]

        molecular_dynamics_quantities = [
            Quantity(
                'initial',
                r' INITIAL\| (.+? {2})=\s+(.+)',
                str_operation=str_to_header,
                repeats=True,
            ),
            Quantity(
                'md_step',
                r'(SCF WAVEFUNCTION OPTIMIZATION[\s\S]+?ENSEMBLE TYPE[\s\S]+?\*{50})',
                repeats=True,
                sub_parser=TextParser(
                    quantities=[
                        Quantity('ensemble_type', r'ENSEMBLE TYPE\s*=\s*(.+)'),
                        Quantity('step', r'STEP NUMBER\s*=\s*(\d+)', dtype=int),
                        Quantity(
                            'time', rf'TIME \[fs\]\s*=\s*({re_float})', dtype=float
                        ),
                        Quantity(
                            'conserved_quantity',
                            rf'CONSERVED QUANTITY \[hartree\]\s*=\s*({re_float})',
                            dtype=float,
                            unit='hartree',
                        ),
                        Quantity(
                            'cpu_time',
                            rf'CPU TIME \[s\]\s*=\s*({re_float})\s*(re_float)',
                            dtype=float,
                        ),
                        Quantity(
                            'energy_drift',
                            rf'ENERGY DRIFT PER ATOM \[K\]\s*=\s*({re_float})\s*({re_float})',
                            dtype=float,
                            unit='hartree',
                        ),
                        Quantity(
                            'potential_energy',
                            rf'POTENTIAL ENERGY\[hartree\]\s*=\s*({re_float})\s*({re_float})',
                            dtype=float,
                            unit='hartree',
                        ),
                        Quantity(
                            'kinetic_energy',
                            rf'KINETIC ENERGY\[hartree\]\s*=\s*({re_float})\s*({re_float})',
                            dtype=float,
                            unit='hartree',
                        ),
                        Quantity(
                            'temperature',
                            rf'TEMPERATURE \[K\]\s*=\s*({re_float})\s*({re_float})',
                            dtype=float,
                        ),
                        Quantity(
                            'pressure',
                            rf'PRESSURE \[bar\]\s*=\s*({re_float})\s*({re_float})',
                            dtype=float,
                            unit='bar',
                        ),
                        Quantity(
                            'barostat_temperature',
                            rf'BAROSTAT TEMP\[K\]\s*=\s*({re_float})\s*({re_float})',
                            dtype=float,
                        ),
                        Quantity(
                            'volume',
                            rf'VOLUME\[bohr\^3\]\s*=\s*({re_float})\s*({re_float})',
                            dtype=float,
                            unit='bohr**3',
                        ),
                        Quantity(
                            'cell_length_instantaneous',
                            rf'CELL LNTHS\[bohr\]\s*=\s*({re_float})\s*({re_float})\s*({re_float})',
                            dtype=float,
                        ),
                        Quantity(
                            'cell_length_average',
                            rf'AVE\. CELL LNTHS\[bohr\]\s*=\s*({re_float})\s*({re_float})\s*({re_float})',
                            dtype=float,
                        ),
                        Quantity(
                            'cell_angle_instantaneous',
                            rf'CELL ANGLS\[deg\]\s*=\s*({re_float})\s*({re_float})\s*({re_float})',
                            dtype=float,
                        ),
                        Quantity(
                            'cell_angle_average',
                            rf'AVE\. CELL ANGLS\[deg\]\s*=\s*({re_float})\s*({re_float})\s*({re_float})',
                            dtype=float,
                        ),
                        Quantity(
                            'self_consistent',
                            r'(SCF WAVEFUNCTION OPTIMIZATION[\s\S]+?)\*{50}',
                            repeats=True,
                            sub_parser=TextParser(
                                quantities=scf_wavefunction_optimization_quantities
                            ),
                        ),
                    ]
                ),
            ),
        ]

        quickstep_quantities = [
            Quantity(
                'dft',
                r' DFT\| (.+? {2}) +(.+)',
                str_operation=str_to_header,
                repeats=True,
            ),
            Quantity('dft_u', r'(DFT\+U\|)'),
            Quantity('mp2', r'(MP2\|)'),
            Quantity('rpa', r'(RI-RPA\|)'),
            Quantity('functional', r' FUNCTIONAL\| (\S+):', repeats=True),
            Quantity(
                'vdw',
                r' vdW POTENTIAL\| .+?([A-Z]\. [A-Z].+? \(\d+\))',
                flatten=False,
                repeats=True,
            ),
            Quantity(
                'qs',
                r' QS\| ((?:Method|Density cutoff)).*?:( {2}) +(.+)',
                str_operation=str_to_header,
                repeats=True,
            ),
            Quantity(
                'atomic_kind_information',
                r' ATOMIC KIND INFORMATION([\s\S]+?)\n\n\n',
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'atom',
                            r'(ic kind: [A-Z][a-z]?[\s\S]+?)(?:\d+\. Atom|\Z)',
                            repeats=True,
                            sub_parser=TextParser(
                                quantities=[
                                    Quantity('kind_label', r'ic kind:\s*(\w+)'),
                                    Quantity(
                                        'kind_number_of_atoms',
                                        r'Number of atoms:\s*(\d+)',
                                        dtype=int,
                                    ),
                                    Quantity(
                                        'kind_basis_set_name',
                                        r'Orbital Basis Set\s*(.+)',
                                    ),
                                    Quantity(
                                        'basis_set_norm_type',
                                        r'Norm type:\s*(\d+)',
                                        dtype=int,
                                    ),
                                ]
                                + n_orbital_basis_quantities
                            ),
                        )
                    ]
                ),
            ),
            Quantity(
                'total_maximum_numbers',
                r' TOTAL NUMBERS AND MAXIMUM NUMBERS([\s\S]+?)\n\n\n',
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            '%s'
                            % key.lower()
                            .replace('the ', '')
                            .replace(' ', '_')
                            .replace('-', '_'),
                            r'\- %s:\s*(\d+)' % key,
                            dtype=int,
                        )
                        for key in [
                            'Atomic kinds',
                            'Atoms',
                            'Shell sets',
                            'Shells',
                            'Primitive Cartesian functions',
                            'Cartesian basis functions',
                            'Spherical basis functions',
                            'Orbital basis functions',
                            'Local part of the GTH pseudopotential',
                            'Non-local part of the GTH pseudopotential',
                        ]
                    ]
                ),
            ),
            Quantity(
                'atomic_coordinates',
                r'(?i: atomic coordinates(?: in) (angstrom))([\s\S]+?)\n\n\n',
                convert=False,
                str_operation=str_to_atomic_coordinates,
            ),  # TODO: if we always capture angstrom, then no need to extract the units...
            Quantity(
                'scf_parameters',
                r' SCF PARAMETERS([\s\S]+?)\*{79}',
                sub_parser=TextParser(
                    quantities=[
                        Quantity('n_max_iteration', r'max_scf:\s*(\d+)', dtype=int),
                        Quantity(
                            'threshold_energy_change',
                            rf'eps_scf:\s*({re_float})',
                            dtype=float,
                            unit='hartree',
                        ),
                        Quantity(
                            'md',
                            r'( MD\| .+? {2} +.+)',
                            str_operation=md_extract,
                            convert=False,
                            repeats=True,
                        ),
                    ]
                ),
            ),
            # TODO add mp2, rpa, gw
            Quantity(
                'single_point',
                r'SCF WAVEFUNCTION OPTIMIZATION([\s\S]+?)(?:\-{50}\n\s*\-|MD_ENERGIES|\Z)',
                repeats=False,
                sub_parser=TextParser(
                    quantities=scf_wavefunction_optimization_quantities
                ),
            ),
            Quantity(
                'geometry_optimization',
                r'STARTING.+?OPTIMIZATION([\s\S]+?(?:OPTIMIZATION COMPLETED|\Z))',
                sub_parser=TextParser(quantities=geometry_optimization_quantities),
            ),
            Quantity(
                'molecular_dynamics',
                r'(MD_ENERGIES\| Initialization proceeding[\s\S]+?\-{50}\n\s*\-)',
                sub_parser=TextParser(quantities=molecular_dynamics_quantities),
            ),
        ]

        self._quantities = [
            Quantity(
                'dbcsr',
                r' (DBCSR)\| (.+? {2}) +(.+)',
                str_operation=str_to_header,
                repeats=True,
            ),
            Quantity(
                'program',
                r'\*\*\s*PROGRAM ([\s\S]+?)(?:\*\*|\n\n|\Z)',
                str_operation=str_to_program,
                repeats=True,
            ),
            Quantity(
                'cp2k',
                r' CP2K\| (.+? {2}) +(.+)',
                str_operation=str_to_header,
                repeats=True,
            ),
            Quantity(
                'global',
                r' GLOBAL\| (.+? {2}) +(.+)',
                str_operation=str_to_header,
                repeats=True,
            ),
            Quantity(
                'restart',
                r'RESTART INFORMATION\s*\*+\s*\*+([\s\S]+?)\*{79}',
                sub_parser=TextParser(
                    quantities=[
                        Quantity('filename', r'RESTART FILE NAME: (\S+)'),
                        Quantity(
                            'quantities',
                            r'RESTARTED QUANTITIES:\s*\*\s*([\s\S]+?)\Z',
                            str_operation=lambda x: [
                                v.split('*')[1].strip() for v in x.strip().split('\n')
                            ],
                        ),
                    ]
                ),
            ),
            Quantity(
                'lattice_vectors',
                rf' CELL\| Vector [abc] \[angstrom\]:\s*({re_float})\s*({re_float})\s*({re_float})',
                repeats=True,
            ),
            # TODO add restart find example
            Quantity(
                'quickstep',
                r'\.\.\. make the atoms dance([\s\S]+?(?:\-{79}\s*\- +\-|\Z))',
                sub_parser=TextParser(quantities=quickstep_quantities),
            ),
            Quantity(
                'spin_polarized',
                r'\| Spin unrestricted \(spin\-polarized\) Kohn\-Sham calculation *([a-zA-Z]+)',
                repeats=False,
            ),
            Quantity(
                'qs_dftb',
                r'  #####   #####        # ######  ####### ####### ######\s*'
                r' #     # #     #      #  #     # #          #    #     #\s*'
                r' #     # #           #   #     # #          #    #     #\s*'
                r' #     #  #####     #    #     # #####      #    ######\s*'
                r' #   # #       #   #     #     # #          #    #     #\s*'
                r' #    #  #     #  #      #     # #          #    #     #\s*'
                r'  #### #  #####  #       ######  #          #    ######\s*'
                r'([\s\S]+?(?:\-{79}\s*\-|\Z))',
                sub_parser=TextParser(quantities=quickstep_quantities),
            ),
            # TODO add other calculation types
        ]
