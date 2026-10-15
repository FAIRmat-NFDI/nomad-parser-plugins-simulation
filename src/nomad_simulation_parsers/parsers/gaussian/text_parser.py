import re

import numpy as np
from nomad.parsing.file_parser import Quantity, TextParser
from nomad.units import ureg

MIN_NORMAL_MODE_TOKENS = 5  # minimum number of tokens expected per normal-mode row


class GaussianOutReader(TextParser):
    def __init__(self):
        super().__init__()

    def init_quantities(self):
        re_float = r'[\d\.\-\+Ee]+'
        re_float_dexp = r'[\d\.\-\+EeDd]+'
        re_eigs = re.compile(r'(\-?\d*\.\d+\s*)')
        # re_force = re.compile(r'(\d+\s*\-?\d*\..*)')

        def str_to_exp(val_in):
            val = np.array(
                val_in.rstrip('.').upper().replace('D', 'E').split(), dtype=float
            )
            return val[0] if len(val) == 1 else val

        def str_to_orbital_symmetries(val_in):
            val = re.findall(r'(?:Occupied|Virtual)\s*((?:\(.+\)\s*)*)', val_in)
            return [v.replace('(', '').replace(')', '').split() for v in val]

        def str_to_eigenvalues(val_in):
            val = val_in.split()
            spin_index = 0 if val[0] == 'Alpha' else 1
            occ = 0.0 if val[1] == 'virt' else 1.0
            eigs, occs = [[], []], [[], []]
            # gaussian sometimes prints values without spaces so we need a re matching
            eigs[spin_index] = re_eigs.findall(val_in)
            occs[spin_index] = [occ] * len(eigs[spin_index])
            return eigs, occs

        def str_to_normal_modes(val_in):
            val = [row.split() for row in val_in.split('\n')]
            return np.array(
                [row[2:] for row in val if len(row) >= MIN_NORMAL_MODE_TOKENS],
                dtype=float,
            )

        def str_to_force_constants(val_in):
            val = val_in.split('\n')
            fc = []
            for line in val:
                parsed = np.array(line.upper().replace('D', 'E').split(), dtype=float)
                index = int(parsed[0])
                if len(fc) < index:
                    fc.append(parsed[1:])
                else:
                    fc[index - 1] += parsed[1:]
            return np.array(fc, dtype=float)

        def str_to_units(unit: str) -> ureg.Unit:
            """Map native Gaussian units to pint units.
            Assumes lower case string input."""
            conv = (
                (r'AMU', 'amu'),
                (r'Dyne', 'dyne'),
                (r'KM', 'km'),
                (r'Mole', 'mole'),
            )
            for u_gauss, u_pint in conv:
                unit = re.sub(u_gauss, u_pint, unit)
            return ureg(unit)

        def normalize_calc_type(val_in: str) -> str:
            return ' '.join(val_in.split())

        orientation_quantities = [
            Quantity(
                'standard_orientation',
                r'd orientation[\s\S]+?X\s*Y\s*Z\s*\-+\s*([\d\.\s\-]+?)\-{2}',
                convert=False,
                str_operation=lambda x: np.array(
                    [v.split() for v in x.strip().split('\n')], dtype=float
                ),
            ),
            Quantity(
                'input_orientation',
                r't orientation[\s\S]+?X\s*Y\s*Z\s*\-+\s*([\d\.\s\-]+?)\-{2}',
                convert=False,
                str_operation=lambda x: np.array(
                    [v.split() for v in x.strip().split('\n')], dtype=float
                ),
            ),
            Quantity(
                'z_matrix_orientation',
                r'x orientation[\s\S]+?X\s*Y\s*Z\s*\-+\s*([\d\.\s\-]+?)\-{2}',
                convert=False,
                str_operation=lambda x: np.array(
                    [v.split() for v in x.strip().split('\n')], dtype=float
                ),
            ),
        ]

        calculation_quantities = [
            Quantity(
                'energy_total',
                rf'\n *(?:Energy=|Electronic Energy)\s*({re_float_dexp})',
                str_operation=str_to_exp,
                convert=False,
                unit='hartree',
                repeats=True,
            ),
            Quantity('hybrid_xc_coeff1', rf'ScaHFX=\s*({re_float})', dtype=float),
            Quantity(
                'hybrid_xc_coeff2',
                rf'ScaDFX=\s*({re_float}\s*{re_float}\s*{re_float}\s*{re_float})',
                flatten=False,
            ),
            Quantity('mp', r'(E2) ='),
            Quantity(
                'energy_total_mp',
                rf'(?:EUMP2|EUMP3|UMP4\(DQ\)|UMP4\(SDQ\)|UMP4\(SDTQ\)|MP5)\s*=\s*({re_float_dexp})',
                str_operation=str_to_exp,
                convert=False,
                unit='hartree',
                repeats=True,
            ),
            Quantity(
                'mp2_correction_energy',
                rf'E2 =\s*({re_float_dexp})',
                str_operation=str_to_exp,
                convert=False,
                unit='hartree',
                repeats=True,
            ),
            Quantity(
                'mp3_correction_energy',
                rf'E3\s*=\s*({re_float_dexp})',
                str_operation=str_to_exp,
                convert=False,
                unit='hartree',
                repeats=True,
            ),
            Quantity(
                'mp4dq_correction_energy',
                rf'E4\(DQ\)=\s*({re_float_dexp})',
                str_operation=str_to_exp,
                convert=False,
                unit='hartree',
                repeats=True,
            ),
            Quantity(
                'mp4sdq_correction_energy',
                rf'E4\(SDQ\)=\s*({re_float_dexp})',
                str_operation=str_to_exp,
                convert=False,
                unit='hartree',
                repeats=True,
            ),
            Quantity(
                'mp4sdtq_correction_energy',
                rf'E4\(SDTQ\)=\s*({re_float_dexp})',
                str_operation=str_to_exp,
                convert=False,
                unit='hartree',
                repeats=True,
            ),
            Quantity(
                'mp5_correction_energy',
                rf'DEMP5 =\s*({re_float_dexp})',
                str_operation=str_to_exp,
                convert=False,
                unit='hartree',
                repeats=True,
            ),
            Quantity('cc', r'((?:CCSD\(T\)|E\(CORR\)))'),
            Quantity(
                'energy_total_cc',
                rf'(?:\n *CCSD\(T\)|E\(CORR\))\s*=\s*({re_float_dexp})',
                str_operation=str_to_exp,
                convert=False,
                unit='hartree',
                repeats=True,
            ),
            Quantity(
                'ccsd_correction_energy',
                rf'DE\(Corr\)=\s*({re_float_dexp})',
                str_operation=str_to_exp,
                convert=False,
                unit='hartree',
                repeats=True,
            ),
            Quantity('qci', r'(Quadratic Configuration Interaction)'),
            Quantity(
                'energy_total_qci',
                rf'(?:QCISD\(T\)|E\(Z\)|QCISD\(TQ\))=\s*({re_float_dexp})',
                str_operation=str_to_exp,
                convert=False,
                unit='hartree',
                repeats=True,
            ),
            Quantity(
                'qcisd_correction_energy',
                rf'DE\(Z\)=\s*({re_float_dexp})',
                str_operation=str_to_exp,
                convert=False,
                unit='hartree',
                repeats=True,
            ),
            Quantity(
                'qcisdtq_correction_energy',
                rf'DE5\s*=\s*({re_float_dexp})',
                str_operation=str_to_exp,
                convert=False,
                unit='hartree',
                repeats=True,
            ),
            Quantity('ci', r'(\s{2}Configuration Interaction)'),
            Quantity(
                'energy_total_ci',
                rf'E\(CI\)=\s*({re_float_dexp})',
                str_operation=str_to_exp,
                convert=False,
                unit='hartree',
                repeats=True,
            ),
            Quantity(
                'ci_correction_energy',
                rf'DE\(CI\)=\s*({re_float_dexp})',
                str_operation=str_to_exp,
                convert=False,
                unit='hartree',
                repeats=True,
            ),
            Quantity(
                'semiempirical_method',
                r'([-A-Z0-9]+\s*calculation of energy[a-zA-Z,. ]+)',
            ),
            Quantity(
                'semiempirical_energy',
                rf'It=\s*\d+\s*PL=\s*[-+0-9EeDd.]+\s*DiagD=[A-Z]\s*ESCF=\s*({re_float})',
                repeats=True,
                unit='hartree',
                dtype=float,
            ),
            Quantity(
                'molmech_method', r'([a-zA-Z0-9]+\s*calculation of energy[a-z,. ]+)'
            ),
            Quantity(
                'excited_state',
                (
                    rf'Excited State\s*(\d+):\s*\S+\s*({re_float})\s*eV\s*'
                    rf'{re_float}\s*nm\s*'
                    rf'f=({re_float})\s*<[\w*]+>=({re_float})\s*(.*)'
                ),
                repeats=True,
            ),
            Quantity(
                'casscf_energy',
                rf'\(\s*[0-9]+\)\s*EIGENVALUE\s*({re_float})',
                repeats=True,
                dtype=float,
                unit='hartree',
            ),
            Quantity(
                'optimization_completed',
                r'(Optimization (?:completed|stopped))',
                flatten=False,
                convert=False,
            ),
            Quantity(
                'population_analysis',
                r'(Population analysis using the SCF density[\s\S]+?Condensed)',
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'orbital_symmetries',
                            r'Orbital symmetries:([\s\S]*?)electro',
                            repeats=False,
                            str_operation=str_to_orbital_symmetries,
                        ),
                        Quantity(
                            'x_gaussian_elstate_symmetry',
                            r'nic state is\s*(.+)\.',
                            flatten=False,
                        ),
                        Quantity(
                            'eigenvalues',
                            r'(Alpha|Beta)\s*(occ|virt)\. eigenvalues \-\-\s*(.+)',
                            repeats=True,
                            str_operation=str_to_eigenvalues,
                            convert=False,
                        ),
                    ]
                ),
            ),
            Quantity(
                'charge',
                rf'\n *Charge=\s*({re_float})\s*electrons',
                dtype=float,
                unit='elementary_charge',
            ),
            Quantity(
                'dipole',
                ''.join([rf'{c}=\s*([\-\d\.]+)\s*' for c in ['X', 'Y', 'Z']]),
                dtype=float,
                unit='debye',
            ),
            Quantity(
                'quadrupole',
                ''.join(
                    [
                        rf'{c}=\s*([\-\d\.]+)\s*'
                        for c in ['XX', 'YY', 'ZZ', 'XY', 'XZ', 'YZ']
                    ]
                ),
                dtype=float,
            ),
            Quantity(
                'octapole',
                r''.join(
                    [
                        rf'{c}=\s*([\-\d\.]+)\s*'
                        for c in [
                            'XXX',
                            'YYY',
                            'ZZZ',
                            'XYY',
                            'XXY',
                            'XXZ',
                            'XZZ',
                            'YZZ',
                            'YYZ',
                            'XYZ',
                        ]
                    ]
                ),
                dtype=float,
                unit='debye * angstrom**2',
            ),
            Quantity(
                'hexadecapole',
                r''.join(
                    [
                        rf'{c}=\s*([\-\d\.]+)\s*'
                        for c in [
                            'XXXX',
                            'YYYY',
                            'ZZZZ',
                            'XXXY',
                            'XXXZ',
                            'YYYX',
                            'YYYZ',
                            'ZZZX',
                            'ZZZY',
                            'XXYY',
                            'XXZZ',
                            'YYZZ',
                            'XXYZ',
                            'YYXZ',
                            'ZZXY',
                        ]
                    ]
                ),
                dtype=float,
                unit='debye * angstrom**3',
            ),
            Quantity(
                'frequency_unit',
                r'[Hh]armonic frequencies \((\S+)\)',
                str_operation=str_to_units,
            ),
            Quantity(
                'reduced_mass_unit',
                r'reduced masses \((\S+)\)',
                str_operation=str_to_units,
            ),
            Quantity(
                'harmonic_force_constant_unit',
                r'force constants \((\S+)\)',
                str_operation=str_to_units,
            ),
            Quantity(
                'ir_intensity_unit',
                r'IR intensities \((\S+)\)',
                str_operation=str_to_units,
            ),
            Quantity(
                'frequencies',
                r'Frequencies\s[\-]{2}([\s\d\.]+)\n',
                dtype=np.float64,
                repeats=True,
            ),  # note the mandatory space after the '--'.
            # Use nested strategy if space is optional
            Quantity(
                'reduced_masses',
                r'Red\. masses\s[\-]{2}([\s\d\.]+)\n',
                dtype=np.float64,
                repeats=True,
            ),  # note the mandatory space after the '--'.
            # Use nested strategy if space is optional
            Quantity(
                'harmonic_force_constants',
                r'Frc consts[\s]{2}[\-]{2}([\s\d\.]+)\n',
                dtype=np.float64,
                repeats=True,
            ),  # note the mandatory space after the '--'.
            # Use nested strategy if space is optional
            Quantity(
                'ir_intensities',
                r'IR Inten[\s]{4}[\-]{2}([\s\d\.]+)\n',
                dtype=np.float64,
                repeats=True,
            ),  # note the mandatory space after the '--'.
            # Use nested strategy if space is optional
            Quantity(
                'normal_modes',
                r'Atom\s*AN.*\s*([\-\d\s\.]+)',
                str_operation=str_to_normal_modes,
                convert=False,
                repeats=True,
            ),
            Quantity(
                'temperature_pressure',
                rf'Temperature\s*({re_float})\s*Kelvin\.\s*Pressure\s*({re_float})\s*Atm\.',
            ),
            Quantity(
                'moments',
                (
                    r'(?:Eigenvalues|EIGENVALUES) \-\-\s*'
                    r'(\d+\.\d{5})\s*(\d+\.\d{5})\s*(\d+\.\d{5})'
                ),
                dtype=float,
                unit='amu*angstrom**2',
            ),
            Quantity(
                'zero_point_energy',
                rf'Zero\-point correction=\s*({re_float})',
                dtype=float,
                unit='hartree',
            ),
            Quantity(
                'thermal_correction_energy',
                rf'Thermal correction to Energy=\s*({re_float})',
                dtype=float,
                unit='hartree',
            ),
            Quantity(
                'thermal_correction_enthalpy',
                rf'Thermal correction to Enthalpy=\s*({re_float})',
                dtype=float,
                unit='hartree',
            ),
            Quantity(
                'thermal_correction_free_energy',
                rf'Thermal correction to Gibbs Free Energy=\s*({re_float})',
                dtype=float,
                unit='hartree',
            ),
            Quantity(
                'forces',
                (
                    r'Forces \(Hartrees/Bohr\)\s*Number\s*Number\s*X\s*Y\s*Z\s*\-+\s*'
                    r'([\d\s\-\.]+?)\s*\-\-'
                ),
                str_operation=lambda x: np.array(
                    [xi.split()[2:5] for xi in x.split('\n')], dtype=float
                ),
            ),
            Quantity(
                'force_constants',
                r'Force constants in Cartesian coordinates:\s*([\s\S]+?)Force',
                str_operation=str_to_force_constants,
                convert=False,
            ),
            Quantity(
                'scf_iteration',
                r'(cle\s*\d+[\s\S]+?(?:Cy|\n\n|Leave))',
                sub_parser=TextParser(
                    quantities=[
                        Quantity('number', r'cle\s*(\d+)', dtype=int),
                        Quantity(
                            'energy_total_scf_iteration',
                            rf' E=\s*({re_float})',
                            dtype=float,
                            unit='hartree',
                        ),
                        Quantity(
                            'x_gaussian_delta_energy_total_scf_iteration',
                            rf'Delta\-E=\s*({re_float})',
                            dtype=float,
                            unit='hartree',
                        ),
                    ]
                ),
                repeats=True,
            ),
            Quantity(
                'scf_iteration_final',
                r'(SCF Done:[\s\S]+?(?:Leave|\n\n|End))',
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'x_gaussian_single_configuration_calculation_converged',
                            r'(SCF Done)',
                            convert=False,
                            flatten=False,
                        ),
                        Quantity('x_gaussian_hf_detect', r'E\((.+?)\)'),
                        Quantity(
                            'energy_total',
                            rf'({re_float})\s*A\.U\.',
                            dtype=float,
                            unit='hartree',
                        ),
                        Quantity(
                            'x_gaussian_energy_error',
                            rf'Conv\s*=\s*({re_float_dexp})',
                            unit='hartree',
                            convert=False,
                            str_operation=str_to_exp,
                        ),
                        Quantity(
                            ' energy_kinetic_electronic',
                            rf'KE\s*=\s*({re_float_dexp})',
                            unit='hartree',
                            convert=False,
                            str_operation=str_to_exp,
                        ),
                        Quantity(
                            'spin_S2',
                            rf'before annihilation\s*({re_float})',
                            dtype=float,
                        ),
                        Quantity(
                            'x_gaussian_after_annihilation_spin_S2',
                            rf',\s*after\s*({re_float})',
                            dtype=float,
                        ),
                        Quantity(
                            'x_gaussian_perturbation_energy',
                            r'[()A-Z0-9]+\s*=\s*[-+0-9D.]+\s*[()A-Z0-9]+\s*=\s*([-+0-9D.]+)',
                            convert=False,
                            str_operation=str_to_exp,
                            unit='hartree',
                        ),
                    ]
                ),
                repeats=False,
            ),
        ]

        run_quantities = [
            Quantity(
                'x_gaussian_settings_corrected',
                r'\-{10}\s*(#[\s\S]+?)\-{10}',
                convert=False,
                str_operation=lambda x: re.sub(r'\s*\n\s*', '', x).strip(),
            ),
            Quantity('charge', r'Charge =\s*([\-\+\d]+)', dtype=int),
            Quantity('spin_target', r'Multiplicity =\s*([\-\+\d]+)', dtype=int),
            Quantity(
                'lattice_vector',
                r'(?:TV|Tv)\s*0?\s*([\d\. ]+)',
                dtype=float,
                repeats=True,
            ),
            Quantity(
                'x_gaussian_atomic_masses',
                r'IAtWgt=([ \d\.]+)',
                dtype=float,
                repeats=True,
            ),
            Quantity(
                'system',
                (
                    r'((?:Standard|Z-Matrix|Input) orientation:[\s\S]+?)'
                    r'(?:Predicted change in Energy|PREDICTED CHANGE IN ENERGY|'
                    r'Z-Matrix orientation|Normal termination)'
                ),
                repeats=True,
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'calculation',
                            (
                                r'(<S\*\*2> of initial guess=[\s\S]+?'
                                r'(?:Initial guess read from the read\-write|\Z))'
                            ),
                            sub_parser=TextParser(quantities=calculation_quantities),
                            repeats=True,
                        )
                    ]
                    + calculation_quantities
                    + orientation_quantities
                ),
            ),
            Quantity(
                'iteration',
                r'ation Nr\.([\s\S]+?)(?:Iter|\n\n)',
                repeats=True,
                sub_parser=TextParser(quantities=calculation_quantities),
            ),
            Quantity('program_cpu_time', r'Job cpu time:(.*)\.\n', flatten=False),
            Quantity(
                'program_termination_date',
                r'Normal termination of Gaussian\s*\d+\s*at\s*(.*)\.\n',
                flatten=False,
            ),
        ]

        self._quantities = [
            Quantity(
                'program',
                r'\s*Gaussian\s*([0-9]+):\s*(\S+)\s*(\S+)\s*(\d+\-\w+\-\d+)',
                convert=False,
            ),
            Quantity('x_gaussian_chk_file', r'%[Cc]hk=([A-Za-z0-9.]*)', dtype=str),
            Quantity('x_gaussian_memory', r'%[Mm]em=([A-Za-z0-9.]*)', dtype=str),
            Quantity(
                'x_gaussian_number_of_processors',
                r'%[Nn][Pp]roc=([A-Za-z0-9.]*)',
                dtype=str,
            ),
            Quantity(
                'calc_type',
                r'\s-+\n\sGaussian ([\w\s]+)\n',
                str_operation=normalize_calc_type,
                flatten=False,
            ),
            Quantity(
                'run',
                # Capture each Gaussian job block ending in a "Normal termination" line.
                # Accepts files with or without a trailing newline.
                r'(-{10}\s*#[\s\S]+?Normal termination.*(?:\n|$))',
                repeats=True,
                sub_parser=TextParser(quantities=run_quantities),
            ),
        ]
