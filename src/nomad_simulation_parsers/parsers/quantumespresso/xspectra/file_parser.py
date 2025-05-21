import numpy as np
from nomad.parsing.file_parser.text_parser import Quantity, TextParser

from ..common import (
    RE_FLOAT,
    RE_N,
    general_quantities,
    header_quantities,
    tail_quantities,
)


def to_val(val_in: str) -> bool | str:
    val = [v.strip().lower() for v in val_in.strip().split(':')]
    val[1] = val[1] == 'true' if val[1] in ['true', 'false'] else val[1]
    return val


def to_error(val_in: str) -> tuple[int, float]:
    val = val_in.split()
    return int(val[0]), float(val[1])


class XSpectraFileParser(TextParser):
    def init_quantities(self) -> None:
        k_calculation_quantities = [
            Quantity(
                'k_point',
                rf'{RE_FLOAT},* +{RE_FLOAT},* +{RE_FLOAT}',
                dtype=np.float64,
            ),
            Quantity(
                'norm_initial_vector',
                rf'[Nn]orm.+?vector[=:] *{RE_FLOAT}',
                dtype=float,
            ),
            Quantity(
                'converged',
                r'(=\> CONVERGED)',
                str_operation=lambda x: True,
            ),
            Quantity(
                'converged',
                r'(not converged)',
                str_operation=lambda x: False,
            ),
            Quantity(
                'n_iter_error',
                rf'iter +(\d+) with error= *({RE_FLOAT})',
                str_operation=to_error,
            ),
            Quantity(
                'n_iter_error',
                rf'final error after *(\d+) *iterations\: ({RE_FLOAT})',
                str_operation=to_error,
            ),
        ]

        self._quantities = [
            Quantity(
                'header',
                r'(Program XSpectra[\s\S]+?)Starting',
                sub_parser=TextParser(
                    quantities=header_quantities
                    + general_quantities
                    + [
                        Quantity(
                            'xspectra_calculation',
                            r'calculation\: *(.+)',
                            flatten=False,
                            dtype=str,
                        ),
                        Quantity(
                            'xspectra_xepsilon',
                            rf'xepsilon +\[.+?\]\: +'
                            rf'({RE_FLOAT}) +({RE_FLOAT}) +({RE_FLOAT})',
                            dtype=np.float64,
                        ),
                        Quantity(
                            'xspectra_xonly_plot',
                            r'xonly_plot\: *(\S)',
                            str_operation=lambda x: x == 'T',
                        ),
                        Quantity(
                            'xspectra_filecore',
                            r'filecore \(core-wavefunction file\): *(\S+)',
                            dtype=str,
                        ),
                        Quantity(
                            'xspectra_main_plot_parameters',
                            rf'main plot parameters\:\s+([\s\S]+?){RE_N} *{RE_N}',
                            sub_parser=TextParser(
                                quantities=[
                                    Quantity(
                                        'key_val',
                                        r'(\w+) *\[*.*\]*(\:) *(\S+)',
                                        repeats=True,
                                        str_operation=to_val,
                                    )
                                ]
                            ),
                        ),
                        Quantity('homo', rf'ehomo *\[eV\]: *{RE_FLOAT}', dtype=float),
                        Quantity('lumo', rf'elumo *\[eV\]: *{RE_FLOAT}', dtype=float),
                        Quantity(
                            'fermi_energy', rf'ef *\[eV\]: *{RE_FLOAT}', dtype=float
                        ),
                        Quantity(
                            'potential_file',
                            r'The potential is recalculated from file :\s+(\S+)',
                            dtype=str,
                        ),
                    ]
                ),
            ),
            Quantity(
                'xanes',
                r'Starting XANES calculation([\s\S]+?xanes +: +)',
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'algorithm',
                            r'\s*Method of calculation based on the\s*([a-zA-Z\s]*) '
                            r'algorithm',
                            repeats=False,
                        ),
                        Quantity(
                            'step_1',
                            r'(Begin STEP 1 [\s\S]+?End STEP 1)',
                            sub_parser=TextParser(
                                quantities=[
                                    Quantity(
                                        'k_calculation',
                                        rf'((?:k\-point *# *1:|k\-point *:*1\s)'
                                        rf'[\s\S]+?){RE_N} *{RE_N}',
                                        repeats=True,
                                        sub_parser=TextParser(
                                            quantities=k_calculation_quantities
                                        ),
                                    )
                                ]
                            ),
                        ),
                        Quantity(
                            'step_2',
                            r'(Begin STEP 2 [\s\S]+?End STEP 2)',
                            sub_parser=TextParser(
                                quantities=[
                                    Quantity(
                                        'energy_zero',
                                        rf'(?:energy\-zero of the spectrum|xe0) \[eV\]'
                                        rf': *({RE_FLOAT})',
                                        dtype=float,
                                        unit='eV',
                                    ),
                                    Quantity(
                                        'xemin',
                                        rf'xemin \[eV\]: *({RE_FLOAT})',
                                        dtype=float,
                                        unit='eV',
                                    ),
                                    Quantity(
                                        'xemax',
                                        rf'xemax \[eV\]: *({RE_FLOAT})',
                                        dtype=float,
                                        unit='eV',
                                    ),
                                    Quantity(
                                        'xnepoint', r'xnepoint: *(\d+)', dtype=int
                                    ),
                                    Quantity(
                                        'broadening_parameter',
                                        rf'constant broadening parameter \[eV\]:'
                                        rf' *{RE_FLOAT}',
                                        dtype=float,
                                    ),
                                    Quantity(
                                        'energy_core_level',
                                        rf'Core level energy \[eV\]: *{RE_FLOAT}',
                                        dtype=float,
                                        unit='eV',
                                    ),
                                    Quantity(
                                        'file',
                                        r'Cross-section successfully written in (\S+)',
                                    ),
                                ]
                            ),
                        ),
                    ]
                ),
            ),
        ] + tail_quantities
