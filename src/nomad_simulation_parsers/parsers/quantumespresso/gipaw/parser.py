from nomad.datamodel import EntryArchive
from nomad.utils import get_logger

from nomad_simulation_parsers.schema_packages.quantumespresso import gipaw

from ..parser import MainfileTextParser, MainfileXMLParser, QuantumEspressoArchiveWriter
from .file_parser import GIPAWFileParser

LOGGER = get_logger(__name__)


class GIPAWMainfileTextParser(MainfileTextParser):
    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER


class GIPAWMainfileXMLParser(MainfileXMLParser):
    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER


class GIPAWArchiveWriter(QuantumEspressoArchiveWriter):
    schema = gipaw
    _text_parser = GIPAWMainfileTextParser(text_parser=GIPAWFileParser())
    _xml_parser = GIPAWMainfileXMLParser()

    def parse_program(self, archive: EntryArchive, index: int) -> None:
        super().parse_program(archive, index)
