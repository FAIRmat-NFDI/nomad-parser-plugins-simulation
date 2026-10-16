from nomad.datamodel import EntryArchive

from nomad_simulation_parsers.parsers.quantumespresso.parser import (
    QuantumEspressoArchiveWriter,
)
from nomad_simulation_parsers.schema_packages.quantumespresso import xspectra

from ..parser import MainfileTextParser
from .file_parser import XSpectraFileParser


class XspectraMainfileParser(MainfileTextParser):
    pass


class XSpectraArchiveWriter(QuantumEspressoArchiveWriter):
    schema = xspectra
    mainfile_parser = XspectraMainfileParser(text_parser=XSpectraFileParser())

    def parse_program(self, archive: EntryArchive, index: int) -> None:
        super().parse_program(archive, index)
