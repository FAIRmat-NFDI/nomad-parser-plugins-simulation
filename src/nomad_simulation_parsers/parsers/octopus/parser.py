from importlib import reload

from nomad.datamodel import EntryArchive
from nomad.parsing.file_parser import ArchiveWriter
from nomad.parsing.file_parser.mapping_parser import MetainfoParser, TextParser
from nomad.parsing.parser import MatchingParser
from nomad.utils import get_logger
from nomad_simulations.schema_packages.general import Program, Simulation
from structlog.stdlib import BoundLogger

from nomad_simulation_parsers.schema_packages import octopus

from .file_parsers import OutParser

LOGGER = get_logger(__name__)


class OctopusMainfileParser(TextParser):
    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER

    def get_header(self, header: list[list[str]]) -> dict[str, str]:
        return {key: val for key, val in header.get('options', [])}


class OctopusMetainfoParser(MetainfoParser):
    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER


class OctopusArchiveWriter(ArchiveWriter):
    mainfile_parser = OctopusMainfileParser(text_parser=OutParser())
    archive_parser = OctopusMetainfoParser()

    def write_to_archive(self) -> None:
        # Reload the octopus package to update the mapping annotations
        reload(octopus)

        self.mainfile_parser.filepath = self.mainfile
        self.archive.data = Simulation(program=Program(name='Octopus'))

        self.archive_parser.data_object = self.archive.data
        self.archive_parser.annotation_key = 'out'

        self.mainfile_parser.convert(self.archive_parser)


class OctopusParser(MatchingParser):
    archive_writer = OctopusArchiveWriter()

    """
    Main parser interface to NOMAD.
    """

    def parse(
        self,
        mainfile: str,
        archive: EntryArchive,
        logger: BoundLogger,
        child_archives: dict[str, EntryArchive] = {},
    ) -> None:
        self.archive_writer.write(mainfile, archive, logger, child_archives)
