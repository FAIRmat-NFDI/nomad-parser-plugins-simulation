import re
from typing import Any

import numpy as np
import pint
from nomad.parsing.file_parser import Quantity, TextParser
from nomad.units import ureg

RE_FLOAT = r'[-+]?\d+\.\d*(?:[Ee][-+]\d+)?'
RE_SYMBOL = re.compile(r'([A-Z][a-z]?)')


def str_to_array(val_in: str) -> np.ndarray:
    """
    Converts a string block to a numpy array of floats.
    """
    val = [v.split(':')[-1].split() for v in val_in.strip().split('\n')]
    val = val[0] if len(val) == 1 else val
    return np.array(val, dtype=float)


def str_to_atom_properties_dict(val_in: str) -> dict[str, Any]:
    """
    Reads the atom properties from a string block.
    """
    unit = None
    if 'charge' in val_in:
        unit = ureg.elementary_charge
    elif 'moment' in val_in:
        unit = ureg.elementary_charge * ureg.bohr

    val = val_in.strip().split('\n')

    properties = dict()
    atom_resolved = []
    species = None
    n_min = 2
    for val_n in val:
        v = val_n.strip().split(':')
        if len(v) < n_min:
            continue

        elif v[0].startswith('species'):
            species = re.search(RE_SYMBOL, v[-1]).group(1)

        elif v[0].startswith('atom'):
            v[0] = v[0].split()
            v[1] = [float(vi) for vi in v[1].split()]
            v[1] = v[1][0] if len(v[1]) == 1 else v[1]
            if species is None:
                species = v[0][2]
            atom_resolved.append((species, v[1] * unit))

        else:
            vi = [float(vii) for vii in v[1].split()]
            # vi = vi[0] if len(vi) == 1 else vi
            properties[v[0].strip()] = vi * unit

    properties['atom_resolved'] = atom_resolved
    return properties


def strip_parentheses(val_in: str) -> str:
    """
    Strips parentheses from string.
    """
    return val_in.strip().replace('(', '').replace(')', '').split()


def str_to_energy_dict(val_in: str) -> dict[str, pint.Quantity]:
    val = val_in.strip().split('\n')
    energies = dict()
    min_n = 2
    for val_n in val:
        v = val_n.split(':')
        if len(v) < min_n:
            continue
        energies[v[0].strip()] = float(v[1]) * ureg.hartree
    return energies


class InfoFileParser(TextParser):
    def init_quantities(self):
        self._quantities = [
            Quantity(
                'program_version',
                r'\s*EXCITING\s*([\w\-\(\)\. ]+)\s*started',
                repeats=False,
                dtype=str,
                flatten=False,
                str_operation=lambda x: x.strip(),
            ),
            Quantity('hash_id', r'version hash id: +(\S+)', dtype=str),
        ]

        initialization_quantities = [
            Quantity(
                'lattice_vectors',
                r'Lattice vectors\s*[\(cartesian\)]*\s*:\s*([\-0-9\.\s]+)\n',
                str_operation=str_to_array,
                unit=ureg.bohr,
                repeats=False,
                convert=False,
            ),
            Quantity(
                'lattice_vectors_reciprocal',
                r'Reciprocal lattice vectors\s*[\(cartesian\)]*\s*:\s*([\-0-9\.\s]+)\n',
                str_operation=str_to_array,
                unit=1 / ureg.bohr,
                repeats=False,
                convert=False,
            ),
        ]

        self._system_keys_mapping = {
            'x_exciting_unit_cell_volume': ('Unit cell volume', ureg.bohr**3),
            'x_exciting_brillouin_zone_volume': (
                'Brillouin zone volume',
                1 / ureg.bohr**3,
            ),
            'x_exciting_number_of_atoms': ('Total number of atoms per unit cell', None),
            'x_exciting_spin_treatment': ('Spin treatment', None),
            'x_exciting_number_of_bravais_lattice_symmetries': (
                'Number of Bravais lattice symmetries',
                None,
            ),
            'x_exciting_number_of_crystal_symmetries': (
                'Number of crystal symmetries',
                None,
            ),
            'kpoint_grid': (r'k\-point grid', None),
            'kpoint_offset': (r'k\-point offset', None),
            'x_exciting_number_kpoints': (r'Total number of k\-points', None),
            'x_exciting_rgkmax': (r'R\^MT\_min \* \|G\+k\|\_max \(rgkmax\)', None),
            'x_exciting_species_rtmin': (r'Species with R\^MT\_min', None),
            'x_exciting_gkmax': (r'Maximum \|G\+k\| for APW functions', 1 / ureg.bohr),
            'x_exciting_gmaxvr': (
                r'Maximum \|G\| for potential and density',
                1 / ureg.bohr,
            ),
            'x_exciting_gvector_size': (r'G\-vector grid sizes', None),
            'x_exciting_gvector_total': (r'Total number of G\-vectors', None),
            'x_exciting_lmaxapw': (r'   APW functions', None),
            'x_exciting_nuclear_charge': (
                'Total nuclear charge',
                ureg.elementary_charge,
            ),
            'x_exciting_electronic_charge': (
                'Total electronic charge',
                ureg.elementary_charge,
            ),
            'x_exciting_core_charge_initial': (
                'Total core charge',
                ureg.elementary_charge,
            ),
            'x_exciting_valence_charge_initial': (
                'Total valence charge',
                ureg.elementary_charge,
            ),
            'x_exciting_wigner_radius': (r'Effective Wigner radius, r\_s', ureg.bohr),
            'x_exciting_empty_states': ('Number of empty states', None),
            'x_exciting_valence_states': ('Total number of valence states', None),
            'x_exciting_hamiltonian_size': ('Maximum Hamiltonian size', None),
            'x_exciting_pw': (r'Maximum number of plane\-waves', None),
            'x_exciting_lo': (r'Total number of local\-orbitals', None),
        }

        self._method_keys_mapping = {
            'smearing_kind': ('Smearing scheme', None),
            'smearing_width': ('Smearing width', None),
        }

        for name, key_unit in self._system_keys_mapping.items():
            initialization_quantities.append(
                Quantity(
                    name,
                    rf'{key_unit[0]}\s*:\s*([\s\S]*?)\n',
                    unit=key_unit[1],
                    repeats=False,
                )
            )

        for name, key_unit in self._method_keys_mapping.items():
            initialization_quantities.append(
                Quantity(
                    name,
                    rf'{key_unit[1]}\s*:\s*([\s\S]*?)\n',
                    unit=key_unit[1],
                    repeats=False,
                )
            )

        initialization_quantities.append(
            Quantity(
                'species',
                rf'(Species : *\d+ *\(\w+\)[\s\S]+?{RE_FLOAT} *{RE_FLOAT} *{RE_FLOAT}'
                rf'\n\s*\n)',
                repeats=True,
                sub_parser=TextParser(
                    quantities=[
                        Quantity('number', r'Species : *(\d+)', dtype=np.int32),
                        Quantity('symbol', r'\((\w+)\)'),
                        Quantity('file', r'parameters loaded from *: *(.+)'),
                        Quantity('name', r'name *: *(.+)'),
                        Quantity(
                            'nuclear_charge',
                            rf'nuclear charge *: *({RE_FLOAT})',
                            dtype=np.float64,
                            unit=ureg.elementary_charge,
                        ),
                        Quantity(
                            'electronic_charge',
                            rf'electronic charge *: *({RE_FLOAT})',
                            dtype=np.float64,
                            unit=ureg.elementary_charge,
                        ),
                        Quantity(
                            'atomic_mass',
                            rf'atomic mass *: *({RE_FLOAT})',
                            dtype=np.float64,
                            unit=ureg.electron_mass,
                        ),
                        Quantity(
                            'muffin_tin_radius',
                            rf'muffin-tin radius *: *({RE_FLOAT})',
                            dtype=np.float64,
                            unit=ureg.bohr,
                        ),
                        Quantity(
                            'radial_points',
                            rf'radial points in muffin-tin *: *({RE_FLOAT})',
                            dtype=np.int32,
                        ),
                        Quantity(
                            'positions_format',
                            r'atomic positions \((.+?)\)',
                            flatten=False,
                        ),
                        Quantity(
                            'positions',
                            rf'\d+ : *({RE_FLOAT}) *({RE_FLOAT}) *({RE_FLOAT})',
                            repeats=True,
                            dtype=np.dtype(np.float64),
                        ),
                    ]
                ),
            )
        )

        initialization_quantities.append(
            Quantity(
                'potential_mixing',
                r'Using ([\w ]+) potential mixing',
                repeats=False,
                flatten=False,
            )
        )

        initialization_quantities.append(
            Quantity(
                'xc_functional',
                r'(Exchange-correlation type[\s\S]+?\n *\n)',
                sub_parser=TextParser(
                    quantities=[
                        Quantity('type', r'Exchange-correlation type +: +(\S+)'),
                        Quantity(
                            'name_reference',
                            r'\n *(.+?,.+)',
                            str_operation=lambda x: [v.strip() for v in x.split(':')],
                        ),
                        Quantity(
                            'parameters',
                            r'\n *(.+?:.+)',
                            repeats=True,
                            str_operation=lambda x: [v.strip() for v in x.split(':')],
                        ),
                    ]
                ),
            )
        )

        self._quantities.append(
            Quantity(
                'initialization',
                r'(?:All units are atomic|Starting initialization)([\s\S]+?)'
                r'(?:Using|Ending initialization)',
                repeats=False,
                sub_parser=TextParser(quantities=initialization_quantities),
            )
        )

        scf_quantities = [
            Quantity(
                'energy_total',
                r'[Tt]*otal energy\s*:\s*([\-\d\.Ee]+)',
                repeats=False,
                dtype=float,
                unit=ureg.hartree,
            ),
            Quantity(
                'energy_contributions',
                r'(?:Energies|_)([\+\-\s\w\.\:]+?)\n *(?:DOS|Density)',
                str_operation=str_to_energy_dict,
                repeats=False,
                convert=False,
            ),
            Quantity(
                'x_exciting_dos_fermi',
                r'DOS at Fermi energy \(states\/Ha\/cell\)\s*:\s*([\-\d\.Ee]+)',
                repeats=False,
                dtype=float,
                unit=1 / ureg.hartree,
            ),
            Quantity(
                'charge_contributions',
                r'(?:Charges|Electron charges\s*\:*\s*)([\-\s\w\.\:\(\)]+?)\n *[A-Z\+]',
                str_operation=str_to_atom_properties_dict,
                repeats=False,
                convert=False,
            ),
            Quantity(
                'moment_contributions',
                r'(?:Moments\s*\:*\s*)([\-\s\w\.\:\(\)]+?)\n *[A-Z\+]',
                str_operation=str_to_atom_properties_dict,
                repeats=False,
                convert=False,
            ),
        ]

        self._miscellaneous_keys_mapping = {
            'x_exciting_gap': (r'Estimated fundamental gap', ureg.hartree),
            'time_physical': (r'Wall time \(seconds\)', ureg.s),
        }

        for name, key_unit in self._miscellaneous_keys_mapping.items():
            scf_quantities.append(
                Quantity(
                    name,
                    rf'{key_unit[0]}\s*\:*\s*([\-\d\.Ee]+)',
                    repeats=False,
                    unit=key_unit[1],
                )
            )

        self._convergence_keys_mapping = {
            'x_exciting_effective_potential_convergence': (
                r'RMS change in effective potential \(target\)',
                ureg.hartree,
            ),
            'x_exciting_energy_convergence': (
                r'Absolute change in total energy\s*\(target\)',
                ureg.hartree,
            ),
            'x_exciting_charge_convergence': (
                r'Charge distance\s*\(target\)',
                ureg.elementary_charge,
            ),
            'x_exciting_IBS_force_convergence': (
                r'Abs\. change in max\-nonIBS\-force\s*\(target\)',
                ureg.hartree / ureg.bohr,
            ),
        }

        for name, key_unit in self._convergence_keys_mapping.items():
            scf_quantities.append(
                Quantity(
                    name,
                    rf'{key_unit[0]}\s*\:*\s*([\(\)\d\.\-\+Ee ]+)',
                    str_operation=strip_parentheses,
                    unit=key_unit[1],
                    repeats=False,
                )
            )

        module_quantities = [
            Quantity(
                'scf_iteration',
                r'(?:I| i)teration number :([\s\S]+?)(?:\n *\n\+{10}|\+\-{10})',
                sub_parser=TextParser(quantities=scf_quantities),
                repeats=True,
            ),
            Quantity(
                'final',
                r'(?:Convergence targets achieved\. Performing final SCF iteration'
                r'|Reached self-consistent loops maximum)([\s\S]+?)(\n *\n\+{10})',
                sub_parser=TextParser(quantities=scf_quantities),
                repeats=False,
            ),
            Quantity(
                'atomic_positions',
                r'(Atomic positions\s*\([\s\S]+?)\n\n',
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'positions_format', r'Atomic positions\s*\(([a-z]+)\)'
                        ),
                        Quantity(
                            'symbols', r'atom\s*\d+\s*(\w+)', repeats=True, dtype=str
                        ),
                        Quantity(
                            'positions',
                            r'\s*:\s*([\d\.\-]+\s*[\d\.\-]+\s*[\d\.\-]+)',
                            repeats=True,
                            dtype=float,
                        ),
                    ]
                ),
            ),
            Quantity(
                'forces',
                r'Total atomic forces including IBS \(\w+\)\s*\:'
                r'(\s*atom[\-\s\w\.\:]*?)\n *Atomic',
                repeats=False,
                str_operation=str_to_array,
                dtype=float,
                unit=ureg.hartree / ureg.bohr,
            ),
        ]

        self._quantities.append(
            Quantity(
                'groundstate',
                r'(?:Self\-consistent loop started|Groundstate module started)'
                r'([\s\S]+?)Groundstate module stopped',
                sub_parser=TextParser(quantities=module_quantities),
                repeats=False,
            )
        )

        optimization_quantities = [
            Quantity(
                'atomic_positions',
                r'(Atomic positions at this step\s*\([\s\S]+?)\n\n',
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'positions_format',
                            r'Atomic positions at this step\s*\(([a-z]+)\)',
                        ),
                        Quantity(
                            'symbols', r'atom\s*\d+\s*(\w+)', repeats=True, dtype=str
                        ),
                        Quantity(
                            'positions',
                            r'\s*:\s*([\d\.\-]+\s*[\d\.\-]+\s*[\d\.\-]+)',
                            repeats=True,
                            dtype=float,
                        ),
                    ]
                ),
            ),
            Quantity(
                'forces',
                r'Total atomic forces including IBS \(\w+\)\s*\:'
                r'(\s*atom[\-\s\w\.\:]*?)\n *Time',
                repeats=False,
                str_operation=str_to_array,
                convert=False,
                unit=ureg.hartree / ureg.bohr,
            ),
            Quantity('step', r'Optimization step\s*(\d+)', repeats=False, dtype=int),
            Quantity('method', r'method\s*=\s*(\w+)', repeats=False, dtype=str),
            Quantity(
                'n_scf_iterations',
                r'Number of (?:total)* scf iterations\s*\:\s*(\d+)',
                repeats=False,
                dtype=int,
            ),
            Quantity(
                'force_convergence',
                r'Maximum force magnitude\s*\(target\)\s*\:(\s*[\(\)\d\.\-\+Ee ]+)',
                str_operation=strip_parentheses,
                unit=ureg.hartree / ureg.bohr,
                repeats=False,
                dtype=float,
            ),
            Quantity(
                'energy_total',
                r'Total energy at this optimization step\s*\:\s*([\-\d\.Ee]+)',
                unit=ureg.hartree,
                repeats=False,
                dtype=float,
            ),
            Quantity(
                'time_calculation',
                r'Time spent in this optimization step\s*\:\s*([\-\d\.Ee]+)\s*seconds',
                unit=ureg.s,
                repeats=False,
                dtype=float,
            ),
        ]

        self._quantities.append(
            Quantity(
                'structure_optimization',
                r'Structure\-optimization module started([\s\S]+?)Structure'
                r'\-optimization module stopped',
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'optimization_step',
                            r'(Optimization step\s*\d+[\s\S]+?(?:\n *\n\-{10}|Time '
                            r'spent in this optimization step\s*:\s*[\d\.]+ seconds))',
                            sub_parser=TextParser(quantities=optimization_quantities),
                            repeats=True,
                        ),
                        Quantity(
                            'final',
                            r'Force convergence target achieved([\s\S]+?Opt)',
                            sub_parser=TextParser(quantities=scf_quantities),
                            repeats=False,
                        ),
                        Quantity(
                            'atomic_positions',
                            r'(imized atomic positions\s*\([\s\S]+?)\n\n',
                            sub_parser=TextParser(
                                quantities=[
                                    Quantity(
                                        'positions_format',
                                        r'imized atomic positions\s*\(([a-z]+)\)',
                                    ),
                                    Quantity(
                                        'symbols',
                                        r'atom\s*\d+\s*(\w+)',
                                        repeats=True,
                                        dtype=str,
                                    ),
                                    Quantity(
                                        'positions',
                                        r'\s*:\s*([\d\.\-]+\s*[\d\.\-]+\s*[\d\.\-]+)',
                                        repeats=True,
                                        dtype=float,
                                    ),
                                ]
                            ),
                        ),
                        Quantity(
                            'forces',
                            r'Total atomic forces including IBS \(\w+\)\s*\:'
                            r'(\s*atom[\-\s\w\.\:]*?)\n *Atomic',
                            repeats=False,
                            str_operation=str_to_array,
                            dtype=float,
                            unit=ureg.hartree / ureg.bohr,
                        ),
                    ]
                ),
                repeats=False,
            )
        )

        self._quantities.append(
            Quantity(
                'hybrids',
                r'Hybrids module started([\s\S]+?)Hybrids module stopped',
                sub_parser=TextParser(quantities=module_quantities),
            )
        )

        self._quantities.append(
            Quantity(
                'total_time',
                r' Total time spent \(seconds\) +: +([\d\.]+)',
                unit=ureg.s,
                repeats=False,
                dtype=float,
            )
        )

    def get_atom_labels(self, section):
        labels = section.get('symbols')

        if labels is None:
            # we get it by concatenating species symbols
            species = self.get('initialization', {}).get('species', [])
            labels = []
            for specie in species:
                labels += [specie.get('symbol')] * len(specie.get('positions'))
        return labels

    def get_positions_format(self, section):
        positions_format = section.get('positions_format')

        if positions_format is None:
            species = self.get_initialization_parameter('species', [])
            for specie in species:
                positions_format = specie.get('positions_format', None)
                if positions_format is not None:
                    break

        return positions_format

    def get_atom_positions(self, section={}, positions=None, positions_format=None):
        positions = positions if positions is not None else section.get('positions')

        if positions is None:
            species = self.get_initialization_parameter('species', [])
            if species:
                positions = np.vstack([s.get('positions') for s in species])

        if positions is None:
            return

        positions = np.array(positions)
        positions_format = (
            positions_format
            if positions_format is not None
            else self.get_positions_format(section)
        )

        if positions_format == 'lattice':
            cell = self.get_initialization_parameter('lattice_vectors')
            if cell is None:
                return
            positions = np.dot(positions, cell.magnitude)

        return positions * ureg.bohr

    def get_scf_threshold(self, name):
        reference = self.get('groundstate', self.get('hybrids', {}))
        return reference.get('scf_iteration', [{}])[-1].get(name, [None, None])[-1]

    def get_scf_quantity(self, name):
        n_scf = len(self.get('energy_total_scf_iteration', []))
        quantity = self.get(f'{name}_scf_iteration')
        if quantity is None:
            return

        # this is really problematic if some scf steps dont have the quantity
        # the only thing that we can do is to assume that the first steps are the
        # ones with the missing quantity
        if len(quantity) < n_scf:
            quantity = [None] * (n_scf - len(quantity)) + quantity

        return quantity

    def get_xc_functional_name(self):
        # TODO expand list to include other xcf
        xc_functional_map = {
            2: ['LDA_C_PZ', 'LDA_X_PZ'],
            3: ['LDA_C_PW', 'LDA_X_PZ'],
            4: ['LDA_C_XALPHA'],
            5: ['LDA_C_VBH'],
            20: ['GGA_C_PBE', 'GGA_X_PBE'],
            21: ['GGA_C_PBE', 'GGA_X_PBE_R'],
            22: ['GGA_C_PBE_SOL', 'GGA_X_PBE_SOL'],
            26: ['GGA_C_PBE', 'GGA_X_WC'],
            30: ['GGA_C_AM05', 'GGA_C_AM05'],
            300: ['GGA_C_BGCP', 'GGA_X_PBE'],
            406: ['HYB_GGA_XC_PBEH'],
            408: ['HYB_GGA_XC_HSE03'],
        }

        xc_functional = self.get('initialization', {}).get('xc_functional', None)
        if xc_functional is None:
            return []

        name = xc_functional_map.get(xc_functional.type, [])

        return name

    @property
    def n_optimization_steps(self):
        return len(self.get('structure_optimization', {}).get('optimization_step', []))

    def get_number_of_spin_channels(self):
        spin_treatment = self.get('initialization', {}).get(
            'x_exciting_spin_treatment', 'spin-unpolarised'
        )
        n_spin = 1 if spin_treatment.lower() == 'spin-unpolarised' else 2
        return n_spin

    def get_unit_cell_volume(self):
        return self.get('initialization', {}).get(
            'x_exciting_unit_cell_volume', 1.0 * ureg.bohr**3
        )

    def get_initialization_parameter(self, key, default=None):
        return self.get('initialization', {}).get(key, default)
