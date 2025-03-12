from nomad.datamodel import EntryArchive

from nomad_simulation_parsers.parsers.quantumespresso.parser import (
    QuantumEspressoArchiveWriter,
)
from nomad_simulation_parsers.schema_packages.quantumespresso import phonon

from ..parser import MainfileParser
from .file_parser import PhononFileParser


class PhononMainfileParser(MainfileParser):
    pass


class PhononArchiveWriter(QuantumEspressoArchiveWriter):
    schema = phonon
    mainfile_parser = PhononMainfileParser(text_parser=PhononFileParser())

    def parse_program(self, archive: EntryArchive, index: int) -> None:
        self.simulation_parser.annotation_key = 'out'
        super().parse_program(archive, index)
