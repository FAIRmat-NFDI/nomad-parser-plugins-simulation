from nomad.datamodel import EntryArchive

from nomad_simulation_parsers.parsers.quantumespresso.parser import (
    QuantumEspressoArchiveWriter,
)
from nomad_simulation_parsers.schema_packages.quantumespresso import phonon

from ..parser import MainfileTextParser
from .file_parser import PhononFileParser


class PhononMainfileParser(MainfileTextParser):
    pass


class PhononArchiveWriter(QuantumEspressoArchiveWriter):
    schema = phonon
    mainfile_parser = PhononMainfileParser(text_parser=PhononFileParser())

    def parse_program(self, archive: EntryArchive, index: int) -> None:
        super().parse_program(archive, index)
