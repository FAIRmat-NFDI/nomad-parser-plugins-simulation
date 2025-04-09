import numpy as np
from nomad.parsing.file_parser.text_parser import Quantity, TextParser

from ..common import RE_FLOAT, general_quantities, header_quantities, tail_quantities


class EPWFileParser(TextParser):
    def init_quantities(self) -> None:
        self._quantities = [
            Quantity(
                'header',
                r'(Program EPW[\s\S]+?WALL\s+\-{50})',
                sub_parser=TextParser(
                    quantities=header_quantities
                    + general_quantities
                    + [
                        Quantity('restart', r'RESTART \- (RESTART)', dtype=str),
                    ]
                ),
            ),
            Quantity(
                'irreducible_q_point',
                r'(irreducible q point # +\d+\s+\=+[\s\S]+?\={50})',
                repeats=True,
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'n_symmetries',
                            r'Symmetries of small group of q\: (\d+)',
                            dtype=int,
                        ),
                        Quantity(
                            'n_q_star',
                            r'Number of q in the star *= *(\d+)',
                            dtype=int,
                        ),
                        Quantity(
                            'q_star',
                            rf'List of q in the star\:\s+((?:\d+ +'
                            rf'{RE_FLOAT} +{RE_FLOAT} +{RE_FLOAT}\s+)+)',
                            dtype=np.float64,
                        ),
                    ]
                ),
            ),
            Quantity(
                'n_ws_vectors_electrons',
                r'Number of WS vectors for electrons *(\d+)',
                dtype=int,
            ),
            Quantity(
                'n_ws_vectors_phonons',
                r'Number of WS vectors for phonons *(\d+)',
                dtype=int,
            ),
            Quantity(
                'n_ws_vectors_electron_phonon',
                r'Number of WS vectors for electron\-phonon *(\d+)',
                dtype=int,
            ),
            Quantity(
                'n_max_cores',
                r'Maximum number of cores for efficient parallelization *(\d+)',
                dtype=int,
            ),
            Quantity(
                'use_ws',
                r'Results may improve by using use_ws == \.(TRUE)\.',
                str_operation=lambda x: True,
                dtype=bool,
            ),
            Quantity(
                'q_mesh',
                r'Using uniform q\-mesh: *(\d+) +(\d+) +(\d+)',
                dtype=np.int32,
            ),
            Quantity(
                'n_q_mesh',
                r'Size of q point mesh for interpolation: *(\d+)',
                dtype=int,
            ),
            Quantity(
                'k_mesh',
                r'Using uniform k\-mesh: *(\d+) +(\d+) +(\d+)',
                dtype=np.int32,
            ),
            Quantity(
                'n_k_mesh',
                r'Size of k point mesh for interpolation: *(\d+)',
                dtype=int,
            ),
            Quantity(
                'n_max_kpoints_per_pool',
                r'Max number of k points per pool: *(\d+)',
                dtype=int,
            ),
            Quantity(
                'e_fermi_coarse_grid',
                rf'Fermi energy coarse grid = *({RE_FLOAT}) eV',
                dtype=float,
                unit='eV',
            ),
            Quantity(
                'n_electrons',
                rf'The Fermi level will be determined with *({RE_FLOAT}) electrons',
                dtype=float,
            ),
            Quantity(
                'e_fermi',
                rf'Fermi energy is calculated from the fine k\-mesh: '
                rf'Ef = *({RE_FLOAT}) eV',
                dtype=float,
                unit='eV',
            ),
            Quantity(
                'self_energy_migdal_approximation',
                r'(Phonon \(Imaginary\) Self\-Energy in the Migdal Approximation\s+'
                r'\=+[\s\S]+?\={50})',
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'fermi_surface_thickness',
                            rf'Fermi Surface thickness = *({RE_FLOAT}) eV',
                            dtype=float,
                            unit='eV',
                        ),
                        Quantity(
                            'golden_rule_t',
                            rf'Golden Rule strictly enforced with T = *({RE_FLOAT}) eV',
                            dtype=float,
                            unit='eV',
                        ),
                        Quantity(
                            'gaussian_broadening',
                            rf'Gaussian Broadening: *({RE_FLOAT}) eV',
                            dtype=float,
                            unit='eV',
                        ),
                        Quantity('n_gauss', r'gauss *= *(\d+)', dtype=int),
                        Quantity(
                            'dos_ef',
                            rf'DOS = *({RE_FLOAT}) states/spin/eV/Unit Cell',
                            dtype=float,
                        ),
                        Quantity(
                            'e_fermi',
                            rf'at Ef= *({RE_FLOAT}) eV',
                            dtype=float,
                            unit='eV',
                        ),
                        Quantity(
                            'self_energy',
                            r'(ismear = +\d+ iq =.+\s+\-+[\s\S]+?\-{50})',
                            repeats=True,
                            sub_parser=TextParser(
                                quantities=[
                                    Quantity('ismear', r'ismear = *(\d+)', dtype=int),
                                    Quantity('iq', r'iq = *(\d+)', dtype=int),
                                    Quantity(
                                        'coord',
                                        rf'coord\.: +'
                                        rf'{RE_FLOAT} +{RE_FLOAT} +{RE_FLOAT}',
                                        dtype=np.float64,
                                    ),
                                    Quantity('wt', rf'wt: +{RE_FLOAT}', dtype=float),
                                    Quantity(
                                        'temp',
                                        rf'Temp: *({RE_FLOAT}) *K',
                                        dtype=float,
                                        unit='K',
                                    ),
                                    Quantity(
                                        'lambda_gamma_omega',
                                        rf'lambda___\( *\d+ *\)= *({RE_FLOAT}) *'
                                        rf'gamma___= *({RE_FLOAT}) meV *omega= *'
                                        rf'({RE_FLOAT}) meV',
                                        dtype=np.float64,
                                        repeats=True,
                                    ),
                                    Quantity(
                                        'lambda_gamma_omega_tr',
                                        rf'lambda_tr\( *\d+ *\)= *({RE_FLOAT}) *'
                                        rf'gamma_tr= *({RE_FLOAT}) meV *omega= *'
                                        rf'({RE_FLOAT}) meV',
                                        dtype=np.float64,
                                        repeats=True,
                                    ),
                                    Quantity(
                                        'lambda_tot',
                                        rf'lambda___\( *tot *\)= *({RE_FLOAT})',
                                        dtype=float,
                                    ),
                                    Quantity(
                                        'lambda_tot_tr',
                                        rf'lambda_tr\( *tot *\)= *({RE_FLOAT})',
                                        dtype=float,
                                    ),
                                ]
                            ),
                        ),
                    ]
                ),
            ),
            Quantity(
                'eliashberg_spectral_function_migdal_approximation',
                r'(Eliashberg Spectral Function in the Migdal Approximation\s+'
                r'\=+[\s\S]+?\={50})',
                sub_parser=TextParser(
                    quantities=[
                        Quantity('lambda', rf'lambda : *({RE_FLOAT})', dtype=float),
                        Quantity(
                            'lambda_tr', rf'lambda_tr : *({RE_FLOAT})', dtype=float
                        ),
                        Quantity('logavg', rf'logavg = *({RE_FLOAT})', dtype=float),
                        Quantity('l_a2f', rf'l_a2f = *({RE_FLOAT})', dtype=float),
                        Quantity(
                            'mu_tc',
                            rf'mu = *({RE_FLOAT}) Tc = *({RE_FLOAT}) K',
                            dtype=np.float64,
                            repeats=True,
                        ),
                        Quantity(
                            'timing',
                            rf' +(.+?) +\: +({RE_FLOAT})s CPU +({RE_FLOAT})s WALL '
                            rf'\( +(\d+) calls',
                            dtype=np.float64,
                            repeats=True,
                            str_operation=lambda x: x.rsplit(' ', 3),
                        ),
                    ]
                ),
            ),
        ] + tail_quantities
