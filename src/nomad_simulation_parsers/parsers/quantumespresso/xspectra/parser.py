from nomad.datamodel import EntryArchive
from nomad.utils import get_logger

from nomad_simulation_parsers.parsers.quantumespresso.parser import (
    QuantumEspressoArchiveWriter,
)
from nomad_simulation_parsers.schema_packages.quantumespresso import common, xspectra

from ..parser import MainfileParser
from .file_parser import XSpectraFileParser

LOGGER = get_logger(__name__)


class XspectraMainfileParser(MainfileParser):
    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER


class XSpectraArchiveWriter(QuantumEspressoArchiveWriter):
    schema = xspectra
    mainfile_parser = XspectraMainfileParser(text_parser=XSpectraFileParser())

    def parse_program(self, archive: EntryArchive, index: int) -> None:
        self.simulation_parser.annotation_key = common.OUT_KEY
        super().parse_program(archive, index)
