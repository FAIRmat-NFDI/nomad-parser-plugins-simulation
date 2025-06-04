from nomad.parsing.parser import MatchingParser
from nomad.parsing.file_parser import ArchiveWriter
from nomad.datamodel import EntryArchive
from structlog.stdlib import BoundLogger




class GPAWArchiveWriter(ArchiveWriter):
    def write_to_archive(self):
        pass



class GPAWParser(MatchingParser):
    """
    Main parser interface to NOMAD.
    """
    archive_writer = GPAWArchiveWriter()

    def parse(self, mainfile: str, archive: EntryArchive, logger: BoundLogger, child_archives: dict[str, EntryArchive] = {}):
        self.archive_writer.write(mainfile, archive, logger, child_archives)
