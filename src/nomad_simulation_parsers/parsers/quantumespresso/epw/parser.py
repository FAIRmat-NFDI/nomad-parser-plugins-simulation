from nomad.datamodel import EntryArchive
from nomad.utils import get_logger
from nomad_simulation_parsers.parsers.quantumespresso.parser import (
    QuantumEspressoArchiveWriter,
)
from nomad_simulation_parsers.schema_packages.quantumespresso import epw

from ..parser import MainfileParser
from .file_parser import EPWFileParser

LOGGER = get_logger(__name__)


class EPWMainfileParser(MainfileParser):
    @property
    def logger(self):
        return LOGGER


class EPWArchiveWriter(QuantumEspressoArchiveWriter):
    schema = epw
    mainfile_parser = EPWMainfileParser(text_parser=EPWFileParser())

    def parse_program(self, archive: EntryArchive, index: int) -> None:
        self.simulation_parser.annotation_key = 'out'
        super().parse_program(archive, index)
