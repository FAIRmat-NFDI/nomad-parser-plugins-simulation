from nomad.parsing.file_parser.text_parser import TextParser


class TurbomoleOutReader(TextParser):
    """
    Minimal stub for the Turbomole text reader.
    TODO:
        Add ParsedQuantity definitions for geometry, energies, SCF info and so on.
    """

    def __init__(self) -> None:
        super().__init__()

    def parse(self, filepath: str) -> None:
        super().parse(filepath)
