import numpy as np
from nomad.parsing.file_parser.text_parser import Quantity, TextParser

from ..common import (
    RE_FLOAT,
    RE_N,
    calculation_quantities,
    general_quantities,
    header_quantities,
    scf_iteration_quantities,
)


class PhononFileParser(TextParser):
    def init_quantities(self) -> None:
        representation_quantities = [
            Quantity('number', r'tion # *(\d+)', dtype=int),
            Quantity(
                'modes',
                r'modes* #([\d ]+)',
                str_operation=lambda x: x.strip().split(),
                dtype=np.int32,
            ),
            Quantity(
                'scf',
                r'Self\-consistent Calculation([\s\S]+?)'
                r'End of self-consistent calculation',
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'iteration',
                            r'( # *\d+[\s\S]+?(?:\n *iter|End))',
                            repeats=True,
                            sub_parser=TextParser(quantities=scf_iteration_quantities),
                        )
                    ]
                ),
            ),
            Quantity(
                'converged',
                r'(Convergence has been achieved)',
                str_operation=lambda x: True,
            ),
        ]

        self._quantities = [
            Quantity(
                'header',
                r'(Program PHONON[\s\S]+?)Calculation',
                sub_parser=TextParser(
                    quantities=header_quantities + general_quantities
                ),
            ),
            Quantity(
                'calculation',
                rf'(q = +{RE_FLOAT} +{RE_FLOAT} +{RE_FLOAT} *{RE_N}[\s\S]+?)'
                rf'(?:Calculation of|\Z)',
                repeats=True,
                sub_parser=TextParser(
                    quantities=general_quantities
                    + calculation_quantities
                    + [
                        Quantity(
                            'q',
                            r'q = +({re_f} +{re_f} +{re_f})',
                            dtype=np.dtype(np.float64),
                        ),
                        Quantity(
                            'dynamical_matrix',
                            r'(Computing dynamical matrix[\s\S]+?)(?:Diagonalizing|\Z)',
                            sub_parser=TextParser(
                                quantities=general_quantities
                                + calculation_quantities
                                + [
                                    Quantity(
                                        'representation',
                                        r'(tion # *\d+ mode[\s\S]+?)(?:Represent|\Z)',
                                        repeats=True,
                                        sub_parser=TextParser(
                                            quantities=representation_quantities
                                        ),
                                    ),
                                ]
                            ),
                        ),
                        Quantity(
                            'frequencies',
                            rf'freq \( +\d+\) += +{RE_FLOAT} \[THz\] += +'
                            rf'({RE_FLOAT}) \[cm\-1\]',
                            repeats=True,
                            dtype=np.float64,
                        ),
                    ]
                ),
            ),
        ]
