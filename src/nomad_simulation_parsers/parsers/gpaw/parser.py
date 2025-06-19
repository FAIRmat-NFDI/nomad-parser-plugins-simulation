from importlib import reload
from typing import Any

from nomad.datamodel import EntryArchive
from nomad.parsing.file_parser import ArchiveWriter
from nomad.parsing.file_parser.mapping_parser import MappingParser, MetainfoParser
from nomad.parsing.parser import MatchingParser
from nomad.utils import get_logger
from nomad_simulation_parsers.schema_packages import gpaw
from nomad_simulations.schema_packages.general import Simulation
from structlog.stdlib import BoundLogger

from .gpw_parser import GPWFileParser

LOGGER = get_logger(__name__)


class GPWParser(MappingParser):
    file_parser = GPWFileParser()

    @property
    def logger(self):
        return LOGGER

    def to_dict(self):
        if self.data_object is None:
            return {}
        self.data_object.parse()
        return self.data_object.results

    def from_dict(self):
        pass

    def load_file(self):
        self.file_parser.mainfile = self.filepath
        return self.file_parser

    def get_energies(self) -> dict[str, Any]:
        energies = {}
        for key in [
            'total',
            'free',
            'XC',
            'kinetic_electronic',
            'correction_entropy',
        ]:
            value = self.file_parser.apply_unit(
                self.file_parser.parser.get_parameter(f'energy_{key}'), 'energyunit'
            )
            if value is None:
                continue
            if key == 'total':
                energies[key] = value
            else:
                energies.setdefault('contributions', []).append(
                    dict(name=key, value=value)
                )
        return energies

    def get_forces(self) -> dict[str, Any]:
        forces = {}
        energyunit = self.file_parser.apply_unit(1, 'energyunit').units
        lengthunit = self.file_parser.apply_unit(1, 'lengthunit').units
        value = self.file_parser.parser.get_array('atom_forces_free')
        if value is not None:
            forces['value'] = value * energyunit / lengthunit
        return forces

    def get_eigenvalues(self) -> list[dict[str, Any]]:
        eigenvalues = self.file_parser.apply_unit(
            self.file_parser.parser.get_array('eigenvalues'), 'energyunit'
        )
        occupations = self.file_parser.parser.get_array('occupation')
        kpoints = self.file_parser.parser.get_array('kpoints')
        return [
            dict(eigenvalues=eigenvalue, occupations=occupations[n], kpoints=kpoints)
            for n, eigenvalue in enumerate(eigenvalues)
        ]


class GPAWMetainfoParser(MetainfoParser):
    @property
    def logger(self):
        return LOGGER


class GPAWArchiveWriter(ArchiveWriter):
    mainfile_parser = GPWParser()
    archive_parser = GPAWMetainfoParser()

    def write_to_archive(self):
        # reload schema annotations
        reload(gpaw)

        self.mainfile_parser.filepath = self.mainfile
        self.archive_parser.annotation_key = 'gpw'
        self.archive_parser.data_object = Simulation()

        self.mainfile_parser.convert(self.archive_parser)
        self.archive.data = self.archive_parser.data_object


class GPAWParser(MatchingParser):
    """
    Main parser interface to NOMAD.
    """

    archive_writer = GPAWArchiveWriter()

    def parse(
        self,
        mainfile: str,
        archive: EntryArchive,
        logger: BoundLogger,
        child_archives: dict[str, EntryArchive] = {},
    ):
        self.archive_writer.write(mainfile, archive, logger, child_archives)
