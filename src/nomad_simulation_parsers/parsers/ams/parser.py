from importlib import reload

from nomad.datamodel import EntryArchive
from nomad.parsing.file_parser import ArchiveWriter
from nomad.parsing.file_parser.mapping_parser import MetainfoParser, TextParser
from nomad.parsing.parser import MatchingParser
from nomad.utils import get_logger
from nomad_simulations.schema_packages.general import Program, Simulation
from structlog.stdlib import BoundLogger

from nomad_simulation_parsers.schema_packages import ams

from .file_parser import OutParser

LOGGER = get_logger(__name__)


class MainfileParser(TextParser):
    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER


# TODO temporary fix for structlog unable to propagate logger
class AMSMetainfoParser(MetainfoParser):
    @property
    def logger(self):
        return LOGGER


class AMSArchiveWriter(ArchiveWriter):
    mainfile_parser = MainfileParser(text_parser=OutParser())
    metainfo_parser = AMSMetainfoParser()

    def write_to_archive(self):
        # reload schema package to use correct annotations
        reload(ams)

        self.metainfo_parser.annotation_key = 'out'
        self.archive.data = Simulation(program=Program(name='AMS'))
        self.metainfo_parser.data_object = self.archive.data

        self.mainfile_parser.filepath = self.mainfile
        self.mainfile_parser.convert(self.metainfo_parser)


class AMSParser(MatchingParser):
    """
    Main parse interface to NOMAD.
    """

    archive_writer = AMSArchiveWriter()

    def parse(
        self,
        mainfile: str,
        archive: EntryArchive,
        logger: BoundLogger = None,
        child_archives: dict[str, EntryArchive] = {},
    ):
        self.archive_writer.write(mainfile, archive, logger, child_archives)
