from nomad.datamodel import EntryArchive

from nomad_simulation_parsers.parsers.quantumespresso.parser import (
    QuantumEspressoArchiveWriter,
)
from nomad_simulation_parsers.schema_packages.quantumespresso import pwscf

from ..parser import MainfileParser
from .file_parser import PWSCFFileParser


class PWSCFMainfileParser(MainfileParser):
    def get_version(self, name_version: list[str]):
        return ' '.join(name_version[1:]).lstrip('v.')


class PWSCFArchiveWriter(QuantumEspressoArchiveWriter):
    schema = pwscf
    mainfile_parser = PWSCFMainfileParser(text_parser=PWSCFFileParser())

    def parse_program(self, archive: EntryArchive, index: int) -> None:
        self.simulation_parser.annotation_key = 'out'
        super().parse_program(archive, index)
