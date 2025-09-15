from nomad.datamodel import EntryArchive
from nomad.utils import get_logger

from nomad_simulation_parsers.parsers.quantumespresso.parser import (
    QuantumEspressoArchiveWriter,
)
from nomad_simulation_parsers.schema_packages.quantumespresso import common, phonon

from ..parser import MainfileTextParser
from .file_parser import PhononFileParser

LOGGER = get_logger(__name__)


class PhononMainfileParser(MainfileTextParser):
    @property
    def logger(self):
        return LOGGER


class PhononArchiveWriter(QuantumEspressoArchiveWriter):
    schema = phonon
    mainfile_parser = PhononMainfileParser(text_parser=PhononFileParser())

    def parse_program(self, archive: EntryArchive, index: int) -> None:
        self.simulation_parser.annotation_key = common.OUT_KEY
        super().parse_program(archive, index)
