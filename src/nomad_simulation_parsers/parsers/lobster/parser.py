import os
from datetime import datetime
from typing import Any

import ase
import ase.io
from nomad.datamodel import EntryArchive
from nomad.parsing import MatchingParser
from nomad.parsing.file_parser import ArchiveWriter
from nomad.parsing.file_parser.mapping_parser import (
    MappingParser,
    MetainfoParser,
    TextParser,
)
from nomad.utils import get_logger
from nomad_simulations.schema_packages.general import Program, Simulation
from structlog.stdlib import BoundLogger

from nomad_simulation_parsers.parsers.utils.general import search_files
from nomad_simulation_parsers.schema_packages import lobster

from .file_parser import ICOXPLISTParser, OutParser

LOGGER = get_logger(__name__)


class LobsterMainfileParser(TextParser):
    text_parser = OutParser()

    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER

    def to_unix_time(self, datetime_str: str):
        return datetime.strptime(datetime_str, '%Y-%m-%d at %H:%M:%S').timestamp()

    def to_basis_set(self, basis_used: str):
        # checks necessary as LOBSTER 5.1.1 writes basis names now in lower case
        basis_used_lower = basis_used.lower()
        if basis_used_lower == 'pbevaspfit2015':
            positions = [3, 7]
        elif basis_used_lower in ['bunge', 'koga']:
            positions = [0]
        char_list = list(basis_used)
        # Loop through the positions and capitalize them if within bounds
        for pos in positions:
            if 0 <= pos < len(char_list):
                char_list[pos] = char_list[pos].upper()
        # Join the list back into a string
        return ''.join(char_list)


class LobsterStructureParser(MappingParser):
    code_name: str = ''

    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER

    def load_file(self) -> ase.Atoms:
        code_name = self.code_name.lower()
        basedir = os.path.dirname(self.filepath)
        if code_name == 'vasp':
            structure_format = 'vasp'
            files = search_files('CONTCAR*', basedir)
        elif code_name == 'quantum espresso':
            structure_format = 'espresso-in'
            files = search_files('*scf.in', basedir)
        else:
            files = []

        if not files:
            raise self.logger.error('No Lobster structure file found.')
        return ase.io.read(files[0], format=structure_format)

    def to_dict(self):
        if self.data_object is not None:
            return self.data_object.todict()
        return {}

    def from_dict(self, data_dict: dict[str, Any]):
        pass

    def get_atoms(self):
        return [
            dict(number=number, symbol=self.data_object.symbols[n])
            for n, number in enumerate(self.data.get('numbers', []))
        ]


class LobsterICOXPLISTParser(TextParser):
    text_parser = ICOXPLISTParser()
    sources = ['ICOHPLIST', 'ICOOPLIST', 'ICOBILIST']
    version: str = ''

    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER

    def load_file(self):
        self.text_parser.version = self.version
        return super().load_file()


class LobsterMetainfoParser(MetainfoParser):
    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER


class LobsterArchiveWriter(ArchiveWriter):
    code_name = 'LOBSTER'
    mainfile_parser = LobsterMainfileParser()
    metainfo_parser = LobsterMetainfoParser()
    structure_parser = LobsterStructureParser()
    icoxplist_parser = LobsterICOXPLISTParser()

    def write_to_archive(self):
        self.archive.data = Simulation(program=Program(name=self.code_name))
        self.metainfo_parser.data_object = self.archive.data

        # parser mainfile
        self.metainfo_parser.annotation_key = lobster.OUT_KEY
        self.mainfile_parser.filepath = self.mainfile
        self.mainfile_parser.convert(self.metainfo_parser)

        # parse structure
        self.metainfo_parser.annotation_key = lobster.STRUCTURE_KEY
        self.structure_parser.code_name = self.mainfile_parser.data.get(
            'x_lobster_code'
        )
        self.structure_parser.filepath = self.mainfile
        self.structure_parser.convert(self.metainfo_parser)

        # parse ICOXPLIST
        self.metainfo_parser.annotation_key = lobster.ICOXPLIST_KEY
        self.icoxplist_parser.filepath = self.mainfile
        self.icoxplist_parser.version = self.mainfile_parser.data.get('program_version')
        self.icoxplist_parser.convert(self.metainfo_parser)


class LobsterParser(MatchingParser):
    def parse(self, mainfile: str, archive: EntryArchive, logger: BoundLogger) -> None:
        archive_writer = LobsterArchiveWriter()
        archive_writer.write(mainfile, archive, logger)
