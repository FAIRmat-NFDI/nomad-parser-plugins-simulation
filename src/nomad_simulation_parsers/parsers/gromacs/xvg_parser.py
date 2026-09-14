import numpy as np
from nomad_file_parser.text_parser import Quantity, TextParser


class GromacsXvgParser(TextParser):
    def init_quantities(self):
        self._quantities = [
            Quantity('title', r'@\s+title\s+\"(.+?)\"', flatten=False),
            Quantity('xaxis', r'xaxis\s+label\s+\"(.+?)\"', flatten=False),
            Quantity('yaxis', r'yaxis\s+label\s+\"(.+?)\"', flatten=False),
            Quantity(
                'column_headers',
                r'@\s+s\d{1,2}\s+legend\s+\"(.+?)\"',
                repeats=True,
                flatten=False,
            ),
        ]

    def parse(self, key=None):
        super().parse(key)
        # TODO extend DataTextParser so it takes in kwarg comments
        self._results['column_vals'] = np.loadtxt(self.mainfile, comments=['@', '#'])
