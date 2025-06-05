from nomad.datamodel import EntryArchive
from nomad.parsing.file_parser import ArchiveWriter
from nomad.parsing.parser import MatchingParser
from nomad_simulations.schema_packages.general import Program, Simulation
from structlog.stdlib import BoundLogger


class LammpsArchiveWriter(ArchiveWriter):
    def write_to_archive(self):
        self.archive.data = Simulation(program=Program(name='LAMMPS'))
        # TODO extend


class LammpsParser(MatchingParser):
    """
    Main parser interface to NOMAD.
    """

    archive_writer = LammpsArchiveWriter()

    def parse(
        self,
        mainfile: str,
        archive: EntryArchive,
        logger: BoundLogger = None,
        child_archives: dict[str, EntryArchive] = None,
    ):
        self.archive_writer.write(mainfile, archive, logger, child_archives)