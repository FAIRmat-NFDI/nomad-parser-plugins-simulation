from nomad.datamodel import EntryArchive

from nomad_simulation_parsers.parsers.quantumespresso.parser import (
    QuantumEspressoArchiveWriter,
)
from nomad_simulation_parsers.schema_packages.quantumespresso import common, xspectra

from ..parser import MainfileTextParser
from .file_parser import XSpectraFileParser


class XspectraMainfileParser(MainfileTextParser):
    pass


class XSpectraArchiveWriter(QuantumEspressoArchiveWriter):
    schema = xspectra
    mainfile_parser = XspectraMainfileParser(text_parser=XSpectraFileParser())

    def parse_program(self, archive: EntryArchive, index: int) -> None:
        self.simulation_parser.annotation_key = common.OUT_KEY
        super().parse_program(archive, index)
