import re
from typing import Any

import numpy as np
from nomad_file_parser import FileParser, Quantity, TextParser

RE_FLOAT = r'[-+]?\d+\.\d*(?:[Ee][-+]\d+)?'
RE_N = r'\n'


class OutParser(TextParser):
    def init_quantities(self):
        self._quantities = [
            Quantity('program_version', r'^LOBSTER *v([\d\.]+) *', repeats=False),
            Quantity(
                'datetime',
                r'starting on host \S* on '
                r'(\d{4}-\d\d-\d\d\sat\s\d\d:\d\d:\d\d)\s[A-Z]{3,4}',
                repeats=False,
                flatten=False,
            ),
            Quantity(
                'x_lobster_code',
                r'detecting used PAW program... (.*)',
                repeats=False,
                flatten=False,
            ),
            Quantity(
                'x_lobster_basis',
                r'setting up local basis functions\.\.\.\s*(?:WARNING.*\s*)*\s*'
                r'((?:[a-zA-Z]{1,2}\s+\(.+\)(?:\s+\d\S+)+\s+)+)',
                repeats=False,
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'x_lobster_basis_species',
                            r'([a-zA-Z]+){1,2}\s+\(([^)]+)\)((?:\s+\d\S+)+)\s+',
                            repeats=True,
                        )
                    ]
                ),
            ),
            Quantity(
                'spilling',
                r'((?:spillings|abs. )[\s\S]*?charge\s*spilling:\s*\d+\.\d+%)',
                repeats=True,
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'abs_total_spilling',
                            r'abs.\s*total\s*spilling:\s*(\d+\.\d+)%',
                            repeats=False,
                        ),
                        Quantity(
                            'abs_charge_spilling',
                            r'abs.\s*charge\s*spilling:\s*(\d+\.\d+)%',
                            repeats=False,
                        ),
                    ]
                ),
            ),
            Quantity('finished', r'finished in (\d)', repeats=False),
        ]


class _COXPCARParser(TextParser):
    def init_quantities(self):

        def get_data(block: str) -> np.ndarray:
            return np.array(
                [line.strip().split() for line in block.splitlines() if line.strip()],
                dtype=np.float64,
            )

        self._quantities = [
            Quantity(
                'coxp_pairs',
                r'(No\.\d+[^\n]*)',
                repeats=True,
            ),
            Quantity(
                'coxp_lines',
                rf'((?:{RE_N}\s*{RE_FLOAT} +{RE_FLOAT}.+)+)',
                dtype=np.float64,
                repeats=False,
                str_operation=get_data,
            ),
        ]


# Alternative file parser for performance
class COXPCARParser(FileParser):
    pair_re = re.compile(rf'No\.(\d+)\:(.+?)\-\>(.+?)\(({RE_FLOAT})\)')

    def parse(self, key=None):
        if self._results is None:
            self._results = {}

        coxp_pairs = []
        n_header = -1
        with self.open_mainfile_obj() as f:
            for line in f:
                n_header += 1
                match = self.pair_re.match(
                    line.decode('utf-8') if isinstance(line, bytes) else line
                )
                if not match and coxp_pairs:
                    break
                if match:
                    coxp_pairs.append(match.groups())

        bond_pairs = []
        unique = []
        for pair in coxp_pairs:
            if pair[0] not in unique:
                unique.append(pair[0])
                bond_pairs.append(pair)

        self._results['bond_pairs'] = bond_pairs
        self._results['coxp_pairs'] = coxp_pairs

        # Use numpy to read data
        data = np.loadtxt(self.mainfile, skiprows=n_header)
        n_spin = 2 if data.shape[1] == (1 + 4 + (len(coxp_pairs) * 4)) else 1

        self._results['energy'] = data[:, 0]
        data = np.reshape(
            data.T[1:, :], (n_spin, 2 * (len(coxp_pairs) + 1), data.shape[0])
        )
        self._results['total_coxp'] = data[:, 0]
        self._results['total_icoxp'] = data[:, 1]
        self._results['pair_coxp'] = np.transpose(data[:, 2::2], (1, 0, 2))
        self._results['pair_icoxp'] = np.transpose(data[:, 3::2], (1, 0, 2))


class CHARGEParser(TextParser):
    kinds = ['mulliken', 'loewdin']

    def init_quantities(self):
        def to_data(block: str) -> dict[str, Any]:
            data = np.array(
                [line.strip().split() for line in block.splitlines() if line.strip()]
            ).T
            dct = dict(
                indices=data[0].astype(int),
                symbols=data[1],
            )
            for n, kind in enumerate(self.kinds):
                dct[kind] = data[2 + n].astype(float)
            return dct

        self._quantities = [
            Quantity(
                'charges',
                rf'((?:\d+ +[A-Z]\w+ +{RE_FLOAT} +{RE_FLOAT}\s+)+)',
                str_operation=to_data,
            ),
            Quantity('total', rf'total +({RE_FLOAT} +{RE_FLOAT})', dtype=np.float64),
        ]


class ICOXPLISTParser(TextParser):
    def init_quantities(self):

        def get_data(block: str) -> np.ndarray:
            return np.transpose(
                [line.strip().split() for line in block.splitlines() if line.strip()]
            )

        self._quantities = [
            Quantity(
                'data',
                rf'((?:\d+ +\w+ +\w+ +{RE_FLOAT} +.+\s*)+)',
                repeats=True,
                str_operation=get_data,
                convert=False,
            )
        ]

    def reset(self):
        super().reset()
        self.init_quantities()
