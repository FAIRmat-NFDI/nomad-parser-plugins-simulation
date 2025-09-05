from nomad.datamodel import EntryArchive
from nomad.utils import get_logger

from nomad_simulation_parsers.schema_packages.quantumespresso import gipaw

from ..parser import MainfileParser, QuantumEspressoArchiveWriter
from .file_parser import GIPAWFileParser

LOGGER = get_logger(__name__)


class GIPAWMainfileParser(MainfileParser):
    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER


class GIPAWArchiveWriter(QuantumEspressoArchiveWriter):
    schema = gipaw
    mainfile_parser = GIPAWMainfileParser(text_parser=GIPAWFileParser())

    def parse_program(self, archive: EntryArchive, index: int) -> None:
        self.simulation_parser.annotation_key = 'out'
        super().parse_program(archive, index)
