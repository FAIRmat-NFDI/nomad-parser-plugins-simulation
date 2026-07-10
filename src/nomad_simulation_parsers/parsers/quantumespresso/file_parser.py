from nomad_file_parser.text_parser import Quantity, TextParser


class QuantumEspressoFileParser(TextParser):
    def init_quantities(self) -> None:
        self._quantities = [
            Quantity(
                'program',
                r'(Program\s*\w+\s*v[\S\s]+?(?:JOB DONE|\Z))',
                repeats=True,
                flatten=False,
            ),
        ]
