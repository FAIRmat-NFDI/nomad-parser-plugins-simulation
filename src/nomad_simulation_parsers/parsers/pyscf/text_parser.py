from nomad.parsing.file_parser.text_parser import TextParser


class PySCFOutReader(TextParser):
    """
    Minimal stub for the PySCF text reader.

    TODO:
        - Add ParsedQuantity definitions for energies, SCF, geometry, etc.
    """

    def __init__(self) -> None:
        super().__init__()

    def parse(self, filepath: str) -> None:
        super().parse(filepath)
