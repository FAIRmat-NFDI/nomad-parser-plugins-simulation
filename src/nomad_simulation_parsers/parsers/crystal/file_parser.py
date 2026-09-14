from typing import Any

import numpy as np
from nomad.units import ureg
from nomad_file_parser import Quantity, TextParser


def capture(regex: str) -> str:
    return r'(' + regex + r')'


FLT = r'-?(?:\d+\.?\d*|\d*\.?\d+)(?:E[\+-]?\d+)?'  # Floating point number
FLT_C = capture(FLT)  # Captures a floating point number
FLT_CRYSTAL_C = r'(-?\d+(?:.\d+)?\*\*-?.*\d+)'  # Crystal specific floating point syntax
WS = r'\s+'  # Series of white-space characters
INTEGER = r'-?\d+'  # Integer number
INTEGER_C = capture(INTEGER)  # Captures integer number
WORD = r'[a-zA-Z]+'  # A single alphanumeric word
WORD_C = capture(WORD)  # Captures a single alphanumeric word
BR = r'\r?\n'  # Newline that works for both Windows and Unix.
# Crystal can be run on a Windows machine as well.


def to_float(value: str) -> Any:
    """Transforms the Crystal-specific float notation into a floating point
    number.
    """
    base, exponent = value.split('**')
    base = int(base)
    exponent = int(''.join(exponent.split()))
    return pow(base, exponent)


class OutputParser(TextParser):
    def init_quantities(self):
        self._quantities = [
            # Header
            Quantity(
                'datetime',
                rf'(?:Date\:|date)\s+(.*?){BR}',
                str_operation=lambda x: x,
                repeats=False,
            ),
            Quantity(
                'hostname',
                rf'(?:Running on\:|hostname)\s+(.*?){BR}',
                str_operation=lambda x: x,
                repeats=False,
            ),
            Quantity(
                'os',
                rf'(?:system)\s+(.*?){BR}',
                str_operation=lambda x: x,
                repeats=False,
            ),
            Quantity(
                'user',
                rf'user\s+(.*?){BR}',
                str_operation=lambda x: x,
                repeats=False,
            ),
            Quantity(
                'input_path',
                rf'(?:Input data|input data in)\s+(.*?){BR}',
                str_operation=lambda x: x,
                repeats=False,
            ),
            Quantity(
                'output_path',
                rf'(?:Output\:|output data in)\s+(.*?){BR}',
                str_operation=lambda x: x,
                repeats=False,
            ),
            Quantity(
                'executable_path',
                rf'(?:Executable\:|crystal executable in)\s+(.*?){BR}',
                str_operation=lambda x: x,
                repeats=False,
            ),
            Quantity(
                'tmpdir',
                rf'(?:Temporary directory\:|temporary directory)\s+(.*?){BR}',
                str_operation=lambda x: x,
                repeats=False,
            ),
            Quantity(
                'system_type',
                r'(CRYSTAL|SLAB|POLYMER|HELIX|MOLECULE|EXTERNAL|DLVINPUT)',
                repeats=False,
            ),
            Quantity('calculation_type', r'(OPTGEOM|FREQCALC|ANHARM)', repeats=False),
            # Input
            Quantity(
                'dftd3',
                rf'(DFTD3{BR}[\s\S]*?END{BR})',
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'version',
                            r'(VERSION \d)',
                            str_operation=lambda x: x,
                            repeats=False,
                        ),
                    ]
                ),
                repeats=False,
            ),
            Quantity(
                'grimme',
                rf'(GRIMME{BR}[\s\S]*?END{BR})',
                repeats=False,
            ),
            Quantity(
                'dft',
                rf'(DFT{BR}[\w\s]*?END{BR})',
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'exchange',
                            rf'EXCHANGE{BR}(LDA|VBH|BECKE|PBE|PBESOL|mPW91|PWGGA|SOGGA|WCGGA)',
                            repeats=False,
                        ),
                        Quantity(
                            'correlation',
                            rf'CORRELAT{BR}(PZ|VBH|VWN|LYP|P86|PBE|PBESOL|PWGGA|PWLSD|WL)',
                            repeats=False,
                        ),
                        Quantity(
                            'exchange_correlation',
                            r'(SVWN|BLYP|PBEXC|PBESOLXC|SOGGAXC|B3PW|B3LYP|PBE0|PBESOL0|B1WC|WCILYP|B97H|PBE0-13|HYBRID|NONLOCAL|HSE06|HSESOL|HISS|RSHXLDA|wB97|wB97X|LC-WPBE|LC-WPBESOL|LC-WBLYP|M05-2X|M05|M062X|M06HF|M06L|M06|B2PLYP|B2GPPLYP|mPW2PLYP|DHYBRID)',
                            repeats=False,
                        ),
                    ]
                ),
                repeats=False,
            ),
            Quantity(
                'program_version',
                rf'{BR} \*\s+CRYSTAL([\d]+)\s+\*',
                repeats=False,
                dtype=str,
            ),
            Quantity(
                'distribution',
                rf'{BR} \*\s*({WORD} : \d+[\.\d+]*)',
                str_operation=lambda x: x,
                repeats=False,
            ),
            Quantity(
                'start_timestamp',
                rf' EEEEEEEEEE STARTING  DATE\s+(.*? TIME .*?){BR}',
                str_operation=lambda x: x,
                repeats=False,
            ),
            Quantity(
                'title',
                rf' EEEEEEEEEE STARTING  DATE.*?{BR}\s*(.*?){BR}{BR}',
                str_operation=lambda x: x,
                repeats=False,
            ),
            Quantity(
                'hamiltonian_type',
                r' (KOHN-SHAM HAMILTONIAN|HARTREE-FOCK HAMILTONIAN)',
                str_operation=lambda x: x,
                repeats=False,
            ),
            Quantity(
                'xc_out',
                r' \(EXCHANGE\)\[CORRELATION\] FUNCTIONAL:(\([\s\S]+?\)\[[\s\S]+?\])',
                str_operation=lambda x: x,
                repeats=False,
            ),
            Quantity(
                'hybrid_out',
                rf' HYBRID EXCHANGE - PERCENTAGE OF FOCK EXCHANGE\s+{FLT_C}',
                repeats=False,
            ),
            # Geometry optimization settings
            Quantity(
                'initial_trust_radius',
                rf' INITIAL TRUST RADIUS\s+{FLT_C}',
                repeats=False,
            ),
            Quantity(
                'maximum_trust_radius',
                rf' MAXIMUM TRUST RADIUS\s+{FLT_C}',
                repeats=False,
            ),
            Quantity(
                'maximum_gradient_component',
                rf' MAXIMUM GRADIENT COMPONENT\s+{FLT_C}',
                repeats=False,
            ),
            Quantity(
                'rms_gradient_component',
                rf' R\.M\.S\. OF GRADIENT COMPONENT\s+{FLT_C}',
                repeats=False,
            ),
            Quantity(
                'rms_displacement_component',
                rf' R\.M\.S\. OF DISPLACEMENT COMPONENTS\s+{FLT_C}',
                repeats=False,
            ),
            Quantity(
                'geometry_change',
                rf' MAXIMUM DISPLACEMENT COMPONENT\s+{FLT_C}',
                unit=ureg.bohr,
                repeats=False,
            ),
            Quantity(
                'energy_change',
                rf' THRESHOLD ON ENERGY CHANGE\s+{FLT_C}',
                unit=ureg.hartree,
                repeats=False,
            ),
            Quantity(
                'extrapolating_polynomial_order',
                rf' EXTRAPOLATING POLYNOMIAL ORDER{WS}{INTEGER_C}',
                repeats=False,
            ),
            Quantity(
                'max_steps',
                rf' MAXIMUM ALLOWED NUMBER OF STEPS\s+{INTEGER_C}',
                repeats=False,
            ),
            Quantity(
                'sorting_of_energy_points',
                rf'SORTING OF ENERGY POINTS\:\s+{WORD_C}',
                repeats=False,
            ),
            # System
            Quantity(
                'material_type',
                rf' ((?:MOLECULAR|SLAB) CALCULATION){BR}',
                str_operation=lambda x: x,
                repeats=False,
            ),
            Quantity(
                'crystal_family',
                rf' CRYSTAL FAMILY\s*:\s*([\s\S]+?)\s*{BR}',
                str_operation=lambda x: x,
                repeats=False,
            ),
            Quantity(
                'crystal_class',
                rf' CRYSTAL CLASS  \(GROTH - 1921\)\s*:\s*([\s\S]+?)\s*{BR}',
                str_operation=lambda x: x,
                repeats=False,
            ),
            Quantity(
                'space_group',
                rf' SPACE GROUP \(CENTROSYMMETRIC\)\s*:\s*([\s\S]+?)\s*{BR}',
                str_operation=lambda x: x,
                repeats=False,
            ),
            Quantity(
                'dimensionality',
                r' GEOMETRY FOR WAVE FUNCTION - DIMENSIONALITY OF THE SYSTEM\s+(\d)',
                repeats=False,
            ),
            Quantity(
                'lattice_parameters',
                rf' (?:PRIMITIVE CELL - CENTRING CODE\s*[\s\S]*?\s*VOLUME=\s*{FLT} - '
                rf'DENSITY\s*{FLT} g/cm\^3{BR}|PRIMITIVE CELL{BR})'
                rf'\s+A\s+B\s+C\s+ALPHA\s+BETA\s+GAMMA.*\s+'
                rf'{FLT_C}\s+{FLT_C}\s+{FLT_C}\s+{FLT_C}\s+{FLT_C}\s+{FLT_C}',
                shape=(6),
                dtype=np.float64,
                repeats=False,
            ),
            Quantity(
                'labels_positions',
                rf' ATOMS IN THE ASYMMETRIC UNIT\s+{INTEGER} - ATOMS IN THE UNIT CELL:'
                rf'\s+{INTEGER}{BR}'
                rf'\s+ATOM\s+X(?:/A|\(ANGSTROM\))\s+Y(?:/B|\(ANGSTROM\))\s+'
                rf'Z(?:/C|\(ANGSTROM\))\s*{BR}'
                rf' \*+?'
                rf'((?:\s+{INTEGER}\s+(?:T|F)\s+{INTEGER}\s+[\s\S]*?\s+{FLT}\s+{FLT}\s+'
                rf'{FLT}{BR})+)',
                shape=(-1, 7),
                dtype=str,
                repeats=False,
            ),
            Quantity(
                'labels_positions_raw',
                rf'AT\.IRR\.\s+AT\s+AT\.N\.\s+X\s+Y\s+Z\s*{BR}'
                rf'((?:\s+{INTEGER}\s+{INTEGER}\s+{INTEGER}\s+{FLT}\s+{FLT}\s+'
                rf'{FLT}{BR})+)',
                shape=(-1, 6),
                dtype=str,
            ),
            # Used to capture an edited geometry. Can contain
            # substitutions, supercells, deformations etc. in any order.
            Quantity(
                'system_edited',
                r' \*\s+GEOMETRY EDITING([\s\S]+?)T = ATOM BELONGING TO THE ASYMMETRIC '
                r'UNIT',
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'lattice_parameters',
                            rf'A\s+B\s+C\s+ALPHA\s+BETA\s+GAMMA.+'
                            rf'\s+{FLT_C}\s+{FLT_C}\s+{FLT_C}\s+{FLT_C}\s+{FLT_C}\s+'
                            rf'{FLT_C}',
                            shape=(6),
                            dtype=np.float64,
                            repeats=False,
                        ),
                        Quantity(
                            'labels_positions',
                            rf'\s+ATOM\s+X(?:/A|\(ANGSTROM\))\s+Y(?:/B|\(ANGSTROM\))\s+'
                            rf'Z(?:/C|\(ANGSTROM\))\s*{BR}'
                            rf' \*+?'
                            rf'((?:\s+{INTEGER}\s+(?:T|F)\s+{INTEGER}\s+[\s\S]*?\s+'
                            rf'{FLT}\s+{FLT}\s+{FLT}{BR})+)',
                            shape=(-1, 7),
                            dtype=str,
                            repeats=False,
                        ),
                        Quantity(
                            'labels_positions_nanotube',
                            rf'\s+ATOM\s+X/A\s+Y\(ANGSTROM\)\s+Z\(ANGSTROM\)\s+'
                            rf'R\(ANGS\)\s*{BR}'
                            rf' \*+?'
                            rf'((?:\s+{INTEGER}\s+(?:T|F)\s+{INTEGER}\s+[\s\S]*?\s+'
                            rf'{FLT}\s+{FLT}\s+{FLT}\s+{FLT}{BR})+)',
                            shape=(-1, 8),
                            dtype=str,
                            repeats=False,
                        ),
                    ]
                ),
                repeats=False,
            ),
            Quantity(
                'lattice_vectors_restart',
                rf' DIRECT LATTICE VECTOR COMPONENTS \(ANGSTROM\){BR}'
                + rf'\s+{FLT_C}\s+{FLT_C}\s+{FLT_C}{BR}'
                + rf'\s+{FLT_C}\s+{FLT_C}\s+{FLT_C}{BR}'
                + rf'\s+{FLT_C}\s+{FLT_C}\s+{FLT_C}{BR}',
                shape=(3, 3),
                dtype=np.float64,
                repeats=False,
            ),
            Quantity(
                'labels_positions_restart',
                rf'   ATOM N\.AT\.  SHELL    X\(A\)      Y\(A\)      Z\(A\)      EXAD  '
                rf'     N\.ELECT\.{BR}'
                rf' \*+?'
                rf'((?:\s+{INTEGER}\s+{INTEGER}\s+{WORD}\s+{INTEGER}\s+{FLT}\s+{FLT}'
                rf'\s+{FLT}\s+{FLT}\s+{FLT}{BR})+)',
                shape=(-1, 9),
                dtype=str,
                repeats=False,
            ),
            Quantity(
                'symmops',
                rf' NUMBER OF SYMMETRY OPERATORS\s*:\s*(\d){BR}',
                repeats=False,
            ),
            # Method
            Quantity(
                'basis_set',
                rf' \*+?'
                rf'{BR} LOCAL ATOMIC FUNCTIONS BASIS SET{BR}'
                rf' \*+?'
                rf'{BR}   ATOM   X\(AU\)   Y\(AU\)   Z\(AU\)  N. TYPE  EXPONENT'
                rf'  S COEF   P COEF   D/F/G COEF{BR}'
                rf'([\s\S]*?){BR} INFORMATION',
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'basis_sets',
                            rf'({BR}{WS}{INTEGER}{WS}{WORD}{WS}{FLT}{WS}{FLT}{WS}{FLT}'
                            rf'{BR}(?:(?:\s+(?:\d+-\s+)?\d+\s+(?:S|P|SP|D|F|G)\s*{BR}'
                            rf'[\s\S]*?(?:{WS}{FLT}(?:{WS})?{FLT}(?:{WS})?{FLT}'
                            rf'(?:{WS})?{FLT}{BR})+)+)?)',
                            sub_parser=TextParser(
                                quantities=[
                                    Quantity(
                                        'species',
                                        rf'{BR}({WS}{INTEGER}{WS}{WORD}{WS}{FLT}{WS}'
                                        rf'{FLT}{WS}{FLT}{BR})',
                                        repeats=False,
                                    ),
                                    Quantity(
                                        'shells',
                                        rf'(\s+(?:\d+-\s+)?\d+\s+(?:S|P|SP|D|F|G)\s*'
                                        rf'{BR}[\s\S]*?(?:{WS}{FLT}(?:{WS})?{FLT}'
                                        rf'(?:{WS})?{FLT}(?:{WS})?{FLT}{BR})+)',
                                        sub_parser=TextParser(
                                            quantities=[
                                                Quantity(
                                                    'shell_range',
                                                    r'(\s+(?:\d+-\s+)?\d+)',
                                                    str_operation=lambda x: ''.join(
                                                        x.split()
                                                    ),
                                                    repeats=False,
                                                ),
                                                Quantity(
                                                    'shell_type',
                                                    rf'((?:S|P|SP|D|F|G))\s*{BR}',
                                                    str_operation=lambda x: x.strip(),
                                                    repeats=False,
                                                ),
                                                Quantity(
                                                    'shell_coefficients',
                                                    rf'{WS}({FLT})(?:{WS})?({FLT})'
                                                    rf'(?:{WS})?({FLT})(?:{WS})?({FLT})'
                                                    rf'{BR}',
                                                    repeats=True,
                                                    dtype=np.float64,
                                                    shape=(4),
                                                ),
                                            ]
                                        ),
                                        repeats=True,
                                    ),
                                ]
                            ),
                            repeats=True,
                        ),
                    ]
                ),
                repeats=False,
            ),
            Quantity(
                'fock_ks_matrix_mixing',
                rf' INFORMATION \*+.*?\*+.*?\:\s+FOCK/KS MATRIX MIXING SET TO\s+'
                rf'{INTEGER_C}\s+\%{BR}',
                repeats=False,
            ),
            Quantity(
                'coulomb_bipolar_buffer',
                rf' INFORMATION \*+.*?\*+.*?\:\s+COULOMB BIPOLAR BUFFER SET TO\s+'
                rf'{FLT_C} Mb{BR}',
                repeats=False,
            ),
            Quantity(
                'exchange_bipolar_buffer',
                rf' INFORMATION \*+.*?\*+.*?\:\s+EXCHANGE BIPOLAR BUFFER SET TO\s+'
                rf'{FLT_C} Mb{BR}',
                repeats=False,
            ),
            Quantity(
                'toldee',
                rf' INFORMATION \*+ TOLDEE \*+\s*\*+ SCF TOL ON TOTAL ENERGY SET TO\s+'
                rf'{FLT_C}{BR}',
                repeats=False,
            ),
            Quantity(
                'n_atoms_per_cell',
                r' N\. OF ATOMS PER CELL\s+' + INTEGER_C,
                repeats=False,
            ),
            Quantity('n_shells', r' NUMBER OF SHELLS\s+' + INTEGER_C, repeats=False),
            Quantity('n_ao', r' NUMBER OF AO\s+' + INTEGER_C, repeats=False),
            Quantity(
                'n_electrons',
                r' N\. OF ELECTRONS PER CELL\s+' + INTEGER_C,
                repeats=False,
            ),
            Quantity(
                'n_core_electrons',
                r' CORE ELECTRONS PER CELL\s+' + INTEGER_C,
                repeats=False,
            ),
            Quantity(
                'n_symmops',
                r' N\. OF SYMMETRY OPERATORS\s+' + INTEGER_C,
                repeats=False,
            ),
            Quantity(
                'tol_coulomb_overlap',
                r' COULOMB OVERLAP TOL\s+\(T1\) ' + FLT_CRYSTAL_C,
                str_operation=to_float,
                repeats=False,
            ),
            Quantity(
                'tol_coulomb_penetration',
                r' COULOMB PENETRATION TOL\s+\(T2\) ' + FLT_CRYSTAL_C,
                str_operation=to_float,
                repeats=False,
            ),
            Quantity(
                'tol_exchange_overlap',
                r' EXCHANGE OVERLAP TOL\s+\(T3\) ' + FLT_CRYSTAL_C,
                str_operation=to_float,
                repeats=False,
            ),
            Quantity(
                'tol_pseudo_overlap_f',
                r' EXCHANGE PSEUDO OVP \(F\(G\)\)\s+\(T4\) ' + FLT_CRYSTAL_C,
                str_operation=to_float,
                repeats=False,
            ),
            Quantity(
                'tol_pseudo_overlap_p',
                r' EXCHANGE PSEUDO OVP \(P\(G\)\)\s+\(T5\) ' + FLT_CRYSTAL_C,
                str_operation=to_float,
                repeats=False,
            ),
            Quantity(
                'pole_order',
                r' POLE ORDER IN MONO ZONE\s+' + INTEGER_C,
                repeats=False,
            ),
            Quantity(
                'calculation_type',
                rf' TYPE OF CALCULATION \:\s+(.*?{BR}\s+.*?){BR}',
                str_operation=lambda x: ' '.join(x.split()),
                repeats=False,
            ),
            Quantity(
                'xc_functional',
                rf' \(EXCHANGE\)\[CORRELATION\] FUNCTIONAL:(\(.+\)\[.+\]){BR}',
                str_operation=lambda x: x,
                repeats=False,
            ),
            Quantity(
                'cappa',
                rf'CAPPA:IS1\s+{INTEGER_C};IS2\s+{INTEGER_C};IS3\s+{INTEGER_C}; K PTS '
                rf'MONK NET\s+{INTEGER_C}; SYMMOPS:\s*K SPACE\s+{INTEGER_C};G SPACE\s+'
                rf'{INTEGER_C}',
                repeats=False,
            ),
            Quantity(
                'scf_max_iteration',
                r' MAX NUMBER OF SCF CYCLES\s+' + INTEGER_C,
                repeats=False,
            ),
            Quantity(
                'convergenge_deltap',
                r'CONVERGENCE ON DELTAP\s+' + FLT_CRYSTAL_C,
                str_operation=to_float,
                repeats=False,
            ),
            Quantity(
                'weight_f',
                r'WEIGHT OF F\(I\) IN F\(I\+1\)\s+' + INTEGER_C,
                repeats=False,
            ),
            Quantity(
                'scf_threshold_energy_change',
                r'CONVERGENCE ON ENERGY\s+' + FLT_CRYSTAL_C,
                str_operation=to_float,
                repeats=False,
                unit=ureg.hartree,
            ),
            Quantity(
                'shrink',
                r'SHRINK\. FACT\.\(MONKH\.\)\s+('
                + INTEGER
                + WS
                + INTEGER
                + WS
                + INTEGER
                + r')',
                repeats=False,
            ),
            Quantity(
                'n_k_points_ibz',
                r'NUMBER OF K POINTS IN THE IBZ\s+' + INTEGER_C,
                repeats=False,
            ),
            Quantity(
                'shrink_gilat',
                r'SHRINKING FACTOR\(GILAT NET\)\s+' + INTEGER_C,
                repeats=False,
            ),
            Quantity(
                'n_k_points_gilat',
                r'NUMBER OF K POINTS\(GILAT NET\)\s+' + INTEGER_C,
                repeats=False,
            ),
            # SCF
            Quantity(
                'scf_block',
                r' CHARGE NORMALIZATION FACTOR([\s\S]*?) == SCF ENDED',
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'scf_iterations',
                            r'( CHARGE NORMALIZATION FACTOR[\s\S]*? '
                            r'(?:TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT PDIG|TTTTTTTTTTTTTTTTT'
                            r'TTTTTTTTTTTTT MPP_KSPA|== SCF ENDED))',
                            sub_parser=TextParser(
                                quantities=[
                                    Quantity(
                                        'charge_normalization_factor',
                                        rf' CHARGE NORMALIZATION FACTOR{WS}{FLT}{BR}',
                                        repeats=False,
                                    ),
                                    Quantity(
                                        'total_atomic_charges',
                                        rf' TOTAL ATOMIC CHARGES:{BR}(?:{WS}{FLT})+'
                                        rf'{BR}',
                                        repeats=False,
                                    ),
                                    Quantity(
                                        'QGAM',
                                        rf' TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT QGAM        '
                                        rf'TELAPSE{WS}{FLT}{WS}TCPU{WS}{FLT}{BR}',
                                        repeats=False,
                                    ),
                                    Quantity(
                                        'BIEL2',
                                        rf' TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT BIEL2       '
                                        rf' TELAPSE{WS}{FLT}{WS}TCPU{WS}{FLT}{BR}',
                                        repeats=False,
                                    ),
                                    Quantity(
                                        'energy_kinetic',
                                        rf' ::: KINETIC ENERGY\s+{FLT_C}{BR}',
                                        unit=ureg.hartree,
                                        repeats=False,
                                    ),
                                    Quantity(
                                        'energy_ee',
                                        rf' ::: TOTAL E-E\s+{FLT_C}{BR}',
                                        unit=ureg.hartree,
                                        repeats=False,
                                    ),
                                    Quantity(
                                        'energy_en_ne',
                                        rf' ::: TOTAL E-N \+ N-E\s+{FLT_C}{BR}',
                                        unit=ureg.hartree,
                                        repeats=False,
                                    ),
                                    Quantity(
                                        'energy_nn',
                                        rf' ::: TOTAL N-N\s+{FLT_C}{BR}',
                                        unit=ureg.hartree,
                                        repeats=False,
                                    ),
                                    Quantity(
                                        'virial_coefficient',
                                        rf' ::: VIRIAL COEFFICIENT\s+{FLT_C}{BR}',
                                        repeats=False,
                                    ),
                                    Quantity(
                                        'TOTENY',
                                        rf' TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT TOTENY      '
                                        rf'  TELAPSE{WS}{FLT}{WS}TCPU{WS}{FLT}{BR}',
                                        repeats=False,
                                    ),
                                    Quantity(
                                        'integrated_density',
                                        rf' NUMERICALLY INTEGRATED DENSITY{WS}{FLT}'
                                        rf'{BR}',
                                        repeats=False,
                                    ),
                                    Quantity(
                                        'NUMDFT',
                                        rf' TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT NUMDFT     '
                                        rf'   TELAPSE{WS}{FLT}{WS}TCPU{WS}{FLT}{BR}',
                                        repeats=False,
                                    ),
                                    Quantity(
                                        'energies',
                                        rf' CYC{WS}{INTEGER}{WS}ETOT\(AU\){WS}{FLT_C}'
                                        rf'{WS}DETOT{WS}{FLT_C}{WS}tst{WS}{FLT}{WS}PX'
                                        rf'{WS}{FLT}{BR}',
                                        repeats=False,
                                        dtype=np.float64,
                                        unit=ureg.hartree,
                                    ),
                                    Quantity(
                                        'FDIK',
                                        rf' TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT FDIK      '
                                        rf'  TELAPSE{WS}{FLT}{WS}TCPU{WS}{FLT}{BR}',
                                        repeats=False,
                                    ),
                                ]
                            ),
                            repeats=True,
                        ),
                    ]
                ),
                repeats=False,
            ),
            Quantity(
                'number_of_scf_iterations',
                rf' == SCF ENDED - CONVERGENCE ON (?:ENERGY|TESTER)\s+E\(AU\)\s*{FLT}'
                rf'\s*CYCLES\s+{INTEGER_C}',
                repeats=False,
            ),
            Quantity(
                'energy_total',
                rf' TOTAL ENERGY\((?:DFT|HF)\)\(AU\)\(\s*{INTEGER}\)\s*{FLT_C} DE\s*'
                rf'{FLT} (?:tester|tst)\s*{FLT}',
                unit=ureg.hartree,
                repeats=False,
            ),
            # Geometry optimization steps
            Quantity(
                'geo_opt',
                rf'( (?:COORDINATE AND CELL OPTIMIZATION|COORDINATE OPTIMIZATION) - '
                rf'POINT\s+1{BR}'
                rf'[\s\S]*?'
                rf' \*+?'
                rf'{BR}'
                rf'\s*\* OPT END - CONVERGED \* E\(AU\)\:\s+{FLT}\s+POINTS\s+{INTEGER})'
                rf'\s+\*{BR}',
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'geo_opt_step',
                            rf' (?:COORDINATE AND CELL OPTIMIZATION|COORDINATE '
                            rf'OPTIMIZATION) - POINT\s+{INTEGER}{BR}'
                            rf'([\s\S]*?)'
                            rf' ((?:TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT OPTI|\* OPT END).+)',
                            sub_parser=TextParser(
                                quantities=[
                                    Quantity(
                                        'lattice_parameters',
                                        rf' (?:PRIMITIVE CELL - CENTRING CODE [\s\S]*?'
                                        rf'VOLUME=\s*{FLT} - DENSITY\s*{FLT} g/cm\^3'
                                        rf'{BR}|PRIMITIVE CELL{BR})'
                                        rf'         A              B              C    '
                                        rf'       ALPHA      BETA       GAMMA\s*'
                                        rf'{FLT_C}\s+{FLT_C}\s+{FLT_C}\s+{FLT_C}\s+'
                                        rf'{FLT_C}\s+{FLT_C}{BR}',
                                        shape=(6),
                                        dtype=np.float64,
                                        repeats=False,
                                    ),
                                    Quantity(
                                        'labels_positions',
                                        rf'\s+ATOM\s+X(?:/A|\(ANGSTROM\))\s+'
                                        rf'Y(?:/B|\(ANGSTROM\))\s+Z(?:/C|\(ANGSTROM\))'
                                        rf'\s*{BR}'
                                        rf' \*+?'
                                        rf'((?:\s+{INTEGER}\s+(?:T|F)\s+{INTEGER}\s+'
                                        rf'[\s\S]*?\s+{FLT}\s+{FLT}\s+{FLT}{BR})+)',
                                        shape=(-1, 7),
                                        dtype=str,
                                        repeats=False,
                                    ),
                                    Quantity(
                                        'labels_positions_nanotube',
                                        rf'\s+ATOM\s+X/A\s+Y\(ANGSTROM\)\s+Z'
                                        rf'\(ANGSTROM\)\s+R\(ANGS\)\s*{BR}'
                                        rf' \*+?'
                                        rf'((?:\s+{INTEGER}\s+(?:T|F)\s+{INTEGER}\s+'
                                        rf'[\s\S]*?\s+{FLT}\s+{FLT}\s+{FLT}\s+{FLT}'
                                        rf'{BR})+)',
                                        shape=(-1, 8),
                                        dtype=str,
                                        repeats=False,
                                    ),
                                    Quantity(
                                        'energy',
                                        rf' TOTAL ENERGY\({WORD}\)\(AU\)\(\s*'
                                        rf'{INTEGER}\)\s*{FLT_C}',
                                        unit=ureg.hartree,
                                        repeats=False,
                                    ),
                                    Quantity(
                                        'time_physical',
                                        rf'OPT.+? TELAPSE\s+({FLT})',
                                    ),
                                ]
                            ),
                            repeats=True,
                        ),
                        Quantity(
                            'converged',
                            rf' \* OPT END - ([\s\S]*?) \* E\(AU\)\:\s+{FLT}\s+POINTS'
                            rf'\s+{INTEGER}',
                            repeats=False,
                        ),
                    ]
                ),
                repeats=False,
            ),
            # Band structure
            Quantity(
                'band_structure',
                rf' \*+?'
                rf'{BR}'
                rf' \*                                                                 '
                rf'            \*{BR}'
                rf' \*  BAND STRUCTURE                                                 '
                rf'            \*{BR}'
                rf'[\s\S]*?'
                rf' \*  FROM BAND\s+{INTEGER} TO BAND\s+{INTEGER}\s+\*{BR}'
                rf' \*  TOTAL OF\s+{INTEGER} K-POINTS ALONG THE PATH\s+\*{BR}'
                rf' \*                                                                 '
                rf'            \*{BR}'
                rf' \*+?'
                rf'{BR}'
                rf'([\s\S]*?'
                rf' ENERGY RANGE \(A\.U\.\)\s*{FLT} - \s*{FLT} EFERMI\s*{FLT_C}{BR})',
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'segments',
                            rf' (LINE\s+{INTEGER} \( {FLT} {FLT} {FLT}: {FLT} {FLT} '
                            rf'{FLT}\) IN TERMS OF PRIMITIVE LATTICE VECTORS{BR}'
                            rf'\s+{INTEGER} POINTS - SHRINKING_FACTOR\s*{INTEGER}{BR}'
                            rf' CARTESIAN COORD\.\s+\( {FLT} {FLT} {FLT}\):\( {FLT} '
                            rf'{FLT} {FLT}\) STEP\s+{FLT}{BR}{BR}{BR})',
                            sub_parser=TextParser(
                                quantities=[
                                    Quantity(
                                        'start_end',
                                        rf'LINE\s+{INTEGER} \( {FLT_C} {FLT_C} {FLT_C}:'
                                        rf' {FLT_C} {FLT_C} {FLT_C}\) IN TERMS OF '
                                        rf'PRIMITIVE LATTICE VECTORS{BR}',
                                        type=np.float64,
                                        shape=(2, 3),
                                        repeats=False,
                                    ),
                                    Quantity(
                                        'n_steps',
                                        rf'\s+{INTEGER_C} POINTS - ',
                                        repeats=False,
                                    ),
                                    Quantity(
                                        'shrinking_factor',
                                        rf'SHRINKING_FACTOR\s*{INTEGER_C}{BR}',
                                        repeats=False,
                                    ),
                                ]
                            ),
                            repeats=True,
                        ),
                        Quantity(
                            'fermi_energy',
                            rf' ENERGY RANGE \(A\.U\.\)\s*{FLT} - \s*{FLT} EFERMI\s*'
                            rf'{FLT_C}',
                            repeats=False,
                        ),
                    ]
                ),
                repeats=False,
            ),
            # DOS
            Quantity(
                'dos',
                rf' RESTART WITH NEW K POINTS NET{BR}'
                rf'([\s\S]+?'
                rf' TOTAL AND PROJECTED DENSITY OF STATES - FOURIER LEGENDRE METHOD{BR}'
                rf'[\s\S]+?)'
                rf' TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT DOSS        TELAPSE',
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'k_points',
                            rf' \*\*\* K POINTS COORDINATES (OBLIQUE COORDINATES IN '
                            rf'UNITS OF IS = {int}){BR}',
                            repeats=False,
                        ),
                        Quantity(
                            'highest_occupied',
                            rf' TOP OF VALENCE BANDS -    BAND\s*{INTEGER}; K\s*'
                            rf'{INTEGER}; EIG {FLT_C}\s*AU',
                            unit=ureg.hartree,
                            repeats=False,
                        ),
                        Quantity(
                            'lowest_unoccupied',
                            rf' BOTTOM OF VIRTUAL BANDS - BAND\s*{INTEGER}; K\s*'
                            rf'{INTEGER}; EIG\s*{FLT_C}\s*AU',
                            unit=ureg.hartree,
                            repeats=False,
                        ),
                    ]
                ),
                repeats=False,
            ),
            Quantity(
                'end_timestamp',
                rf' EEEEEEEEEE TERMINATION  DATE\s+(.*? TIME .*?){BR}',
                str_operation=lambda x: x,
                repeats=False,
            ),
            # Forces
            Quantity(
                'forces',
                rf' CARTESIAN FORCES IN HARTREE/BOHR \(ANALYTICAL\){BR}'
                rf'   ATOM                     X                   Y                   '
                rf'Z{BR}'
                rf'((?:'
                rf'{WS}'
                rf'{INTEGER}'
                rf'{WS}'
                rf'{INTEGER}'
                rf'{WS}'
                rf'{FLT}'
                rf'{WS}'
                rf'{FLT}'
                rf'{WS}'
                rf'{FLT}'
                rf'{BR})*)',
                shape=(-1, 5),
                dtype=str,
                repeats=False,
            ),
            Quantity(
                'end_timestamp',
                rf' EEEEEEEEEE TERMINATION  DATE\s+(.*? TIME .*?){BR}',
                str_operation=lambda x: x,
                repeats=False,
            ),
            Quantity('time_end', rf'END +TELAPSE +({FLT_C})', dtype=np.float64),
            # Filepaths
            Quantity(
                'f25_filepath1',
                rf'file fort\.25 saved as ([\s\S]+?){BR}',
                str_operation=lambda x: x,
                repeats=False,
            ),
            Quantity(
                'f25_filepath2',
                rf'BAND/MAPS/DOSS data for plotting fort.25 saved as ([\s\S]+?){BR}',
                str_operation=lambda x: x,
                repeats=False,
            ),
        ]


class F25Parser(TextParser):
    def init_quantities(self):
        self._quantities = [
            # Band structure energies
            Quantity(
                'segments',
                rf'(-\%-0BAND\s*{INTEGER}\s*{INTEGER}\s?{FLT}\s?{FLT}\s?{FLT}{BR}'
                rf'\s*{FLT}\s*{FLT}{BR}'
                rf'\s*{INTEGER}\s*{INTEGER}\s*{INTEGER}\s*{INTEGER}\s*{INTEGER}\s*'
                rf'{INTEGER}{BR}'
                rf'(?:\s*{FLT})+)',
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'first_row',
                            rf'-\%-0BAND\s*{INTEGER_C}\s*{INTEGER_C}\s?{FLT_C}\s?'
                            rf'{FLT_C}\s?{FLT_C}{BR}',
                            repeats=False,
                        ),
                        Quantity(
                            'second_row',
                            rf'\s?{FLT_C}\s?{FLT_C}{BR}',
                            repeats=False,
                        ),
                        Quantity(
                            'energies',
                            rf'\s*{INTEGER}\s*{INTEGER}\s*{INTEGER}\s*{INTEGER}\s*'
                            rf'{INTEGER}\s*{INTEGER}{BR}'
                            rf'((?:{FLT}\s?)+)',
                            str_operation=lambda x: x,
                            repeats=False,
                        ),
                    ]
                ),
                repeats=True,
            ),
            # DOS values
            Quantity(
                'dos',
                rf'(-\%-0DOSS\s*{INTEGER}\s*{INTEGER}\s?{FLT}\s?{FLT}\s?{FLT}{BR}'
                rf'\s*{FLT}\s?{FLT}{BR}'
                rf'\s*{INTEGER}\s*{INTEGER}\s*{INTEGER}\s*{INTEGER}\s*{INTEGER}'
                rf'\s*{INTEGER}{BR}'
                rf'(?:\s*{FLT})+)',
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'first_row',
                            rf'-\%-0DOSS\s*{INTEGER_C}\s*{INTEGER_C}\s?{FLT_C}\s?'
                            rf'{FLT_C}\s?{FLT_C}{BR}',
                            repeats=False,
                        ),
                        Quantity(
                            'second_row',
                            rf'\s?{FLT_C}\s?{FLT_C}{BR}',
                            repeats=False,
                        ),
                        Quantity(
                            'values',
                            rf'\s*{INTEGER}\s*{INTEGER}\s*{INTEGER}\s*{INTEGER}\s*'
                            rf'{INTEGER}\s*{INTEGER}{BR}'
                            rf'((?:\s*{FLT})+)',
                            str_operation=lambda x: x,
                            repeats=False,
                        ),
                    ]
                ),
                repeats=False,
            ),
        ]
