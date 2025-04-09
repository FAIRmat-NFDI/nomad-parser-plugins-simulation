import os
import re
from collections.abc import Iterable
from datetime import datetime
from importlib import reload
from types import ModuleType

from nomad.config import config
from nomad.datamodel import EntryArchive
from nomad.datamodel.metainfo.workflow import Link, TaskReference
from nomad.parsing import MatchingParser
from nomad.parsing.file_parser import ArchiveWriter
from nomad.parsing.file_parser.mapping_parser import MetainfoParser, TextParser
from nomad.utils import get_logger
from nomad_simulations.schema_packages.general import Program, Simulation
from nomad_simulations.schema_packages.workflow import (
    SerialWorkflow,
    SimulationWorkflow,
    SinglePoint,
)
from structlog.stdlib import BoundLogger

from nomad_simulation_parsers.parsers.utils.general import search_files
from nomad_simulation_parsers.schema_packages.quantumespresso import common

from .file_parser import QuantumEspressoFileParser

LOGGER = get_logger(__name__)


# TODO temporary fix for structlog unable to propagate logger
class QuantumEspressoMetainfoParser(MetainfoParser):
    @property
    def logger(self):
        return LOGGER


class MainfileParser(TextParser):
    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER

    def get_version(self, name_version: list[str]):
        return ' '.join(name_version[1:]).lstrip('v.')


class QuantumEspressoArchiveWriter(ArchiveWriter):
    """
    Wrapper for the program-specific archive writer.
    """

    schema: ModuleType = common
    simulation_parser = QuantumEspressoMetainfoParser()
    mainfile_parser = MainfileParser(text_parser=QuantumEspressoFileParser())

    def parse_program(self, archive: EntryArchive, index: int) -> None:
        # reload the schema annotations
        reload(self.schema)
        self.simulation_parser.data_object = Simulation(
            program=Program(name='Quantum Espresso')
        )
        # convert
        self.mainfile_parser.convert(self.simulation_parser)
        # set the parsed data to archive
        archive.data = self.simulation_parser.data_object

    def parse_workflow(self) -> None:
        # multi run file
        multirun_workflow_archive = self.child_archives.get('workflow_multirun')
        if multirun_workflow_archive is not None:
            multirun_workflow_archive.workflow2 = SerialWorkflow(
                tasks=[TaskReference(task=self.archive.workflow2)]
            )
            for key, child_archive in self.child_archives.items():
                if key.startswith('workflow'):
                    continue
                multirun_workflow_archive.workflow2.tasks.append(
                    TaskReference(task=child_archive.workflow2)
                )

        # mainfiles in the same upload
        generic_workflow_archive = self.child_archives.get('workflow_generic')
        if generic_workflow_archive is not None:
            from nomad.app.v1.models import MetadataRequired
            from nomad.search import search

            parent_archive = multirun_workflow_archive or self.archive
            # add current archive workflow to generic workflow tasks
            generic_workflow_archive.workflow2 = SimulationWorkflow(
                tasks=[TaskReference(task=parent_archive.workflow2)]
            )

            upload_id = self.archive.metadata.upload_id
            metadata = search(
                owner='visible',
                user_id=self.archive.metadata.main_author.user_id,
                query={'upload_id': upload_id},
                required=MetadataRequired(
                    include=['entry_id', 'mainfile', 'parser_name']
                ),
            ).data
            parent_file = self.mainfile.split('raw/')[-1]
            parent_dir = os.path.dirname(parent_file)
            for result in metadata:
                parser_name = result.get('parser_name')
                # include only qe calculations
                if 'quantumespresso' not in parser_name:
                    continue
                mainfile = result.get('mainfile')
                if not mainfile or mainfile == parent_file:
                    # skip the current mainfile
                    continue
                entry_id = result.get('entry_id')
                if not entry_id:
                    continue
                # link only entries in the same directory or sub-directories
                if mainfile.startswith(parent_dir):
                    entry_archive: EntryArchive = self.archive.m_context.load_archive(
                        entry_id, upload_id, None
                    )
                    # add workflow to generic workflow tasks
                    generic_workflow_archive.workflow2.tasks.append(
                        TaskReference(task=entry_archive.workflow2)
                    )
                    # add parent scf as input to task
                    if entry_archive.workflow2:
                        entry_archive.workflow2.inputs.append(
                            Link(section=parent_archive.workflow2)
                        )

    def write_to_archive(self) -> None:
        def load_writer(header: str) -> QuantumEspressoArchiveWriter:
            if 'pwscf' in header:
                from .pwscf.parser import PWSCFArchiveWriter

                return PWSCFArchiveWriter()
            if 'epw' in header:
                from .epw.parser import EPWArchiveWriter

                return EPWArchiveWriter()
            if 'phonon' in header:
                from .phonon.parser import PhononArchiveWriter

                return PhononArchiveWriter()
            if 'xspectra' in header:
                from .xspectra.parser import XSpectraArchiveWriter

                return XSpectraArchiveWriter()
            return None

        # set up mainfile parser
        self.mainfile_parser.filepath = self.mainfile

        if not self.mainfile_parser.data_object.get('program'):
            return

        for n, program in enumerate(
            self.mainfile_parser.data_object.get('program', [])
        ):
            writer = load_writer(program[:50].lower())
            if writer is None:
                self.logger.error('Parser not found for program.')
                continue
            writer.mainfile_parser.data_object.mainfile = self.mainfile
            # parse only the relevant program
            writer.mainfile_parser.data_object._file_handler = program.encode()

            # write the first program to the main archive, the rest to child archives
            program_name = writer.mainfile_parser.data_object.get('header', {}).get(
                'program_name_version', ['']
            )[0]
            archive = (
                self.archive
                if n == 0
                else self.child_archives.get(f'{n} {program_name}')
            )
            if archive is None:
                self.logger.error('Archive not found for program.')
                continue
            writer.parse_program(archive, n)
            archive.workflow2 = SinglePoint()

        self.parse_workflow()


def sort_qe_files(filenames: list[str]) -> list[tuple[str, datetime]]:
    """
    Sort QE mainfiles based on execution time.
    """
    sorted_files = []
    re_pattern = re.compile(r'starts on *(\w+) *at *([\d ]+\:[\d ]+\:[\d ]+)')
    for name in filenames:
        with open(name) as f:
            head = f.read(config.process.parser_matching_size)
            match = re_pattern.search(head)
            if not match:
                continue
            sorted_files.append(
                (
                    name,
                    datetime.strptime(
                        ''.join(match.groups()).replace(' ', ''), '%d%b%Y%H:%M:%S'
                    ),
                )
            )

    sorted_files.sort(key=lambda x: x[1])
    return sorted_files


class QuantumEspressoParser(MatchingParser):
    """
    Common parser for Quantum Espresso mainfiles including
    PWSCF, Phonon, EPW and XSpectra.
    """

    archive_writer = QuantumEspressoArchiveWriter()

    def is_mainfile(
        self,
        filename: str,
        mime: str,
        buffer: bytes,
        decoded_buffer: str,
        compression: str = None,
    ) -> bool | Iterable[str]:
        is_mainfile = super().is_mainfile(
            filename, mime, buffer, decoded_buffer, compression
        )
        if is_mainfile:
            children = []
            programs = []
            program_re = re.compile(r'Program +(\w+)')
            with open(filename) as f:
                for line in f:
                    match = program_re.search(line)
                    if not match:
                        continue
                    programs.append(f'{len(programs)} {match.group(1)}')
            if 'pwscf' in programs[0].lower():
                self.level = 2
                # search all qe mainfiles in the directory and sub directories
                qe_files = search_files(
                    '*.out', os.path.dirname(filename), include_all=True
                )
                if len(qe_files) > 1:
                    sorted_files = sort_qe_files(qe_files)
                    if sorted_files and sorted_files[0][0] == filename:
                        children.append('workflow_generic')

            if len(programs) > 1:
                # create separate entries for each program instance
                children.extend(['workflow_multirun', *programs[1:]])

            self.creates_children = len(children) > 0

            return children or True

        return is_mainfile

    def parse(
        self,
        mainfile: str,
        archive: EntryArchive,
        logger: BoundLogger,
        child_archives: dict[str, EntryArchive] = {},
    ) -> None:
        print('PPPP')
        self.archive_writer.write(mainfile, archive, logger, child_archives)
