from typing import Any

from nomad.datamodel import EntryArchive
from nomad.parsing.file_parser import ArchiveWriter
from nomad.parsing.file_parser.mapping_parser import MappingParser, MetainfoParser
from nomad.parsing.parser import MatchingParser
from nomad.utils import get_logger
from nomad_simulations.schema_packages.general import Simulation
from nomad_simulations.schema_packages.workflow.general import EnergyConvergenceTarget
from nomad_simulations.schema_packages.workflow.single_point import (
    SinglePoint,
    SinglePointMethod,
)
from structlog.stdlib import BoundLogger

from nomad_simulation_parsers.schema_packages import gpaw

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

    def get_scf_steps(self) -> dict[str, Any]:
        code_specific_quantities = {}
        converged = self.file_parser.parser.get_parameter('converged')
        if converged is not None:
            code_specific_quantities['converged'] = bool(converged)

        energy_error = self.file_parser.parser.get_parameter('energyerror')
        if energy_error is not None:
            code_specific_quantities['energyerror'] = float(energy_error)

        if code_specific_quantities:
            return {'code_specific_quantities': code_specific_quantities}
        return {}

    def build_workflow(self):
        workflow = SinglePoint()
        workflow.method = SinglePointMethod()
        energy_error = self.file_parser.parser.get_parameter('energyerror')
        if energy_error is not None:
            workflow.method.convergence_targets = [
                EnergyConvergenceTarget(
                    threshold=self.file_parser.apply_unit(energy_error, 'energyunit'),
                    threshold_type='absolute',
                )
            ]
        return workflow


class GPAWMetainfoParser(MetainfoParser):
    @property
    def logger(self):
        return LOGGER


class GPAWArchiveWriter(ArchiveWriter):
    mainfile_parser = GPWParser()
    archive_parser = GPAWMetainfoParser()

    def write_to_archive(self):
        self.mainfile_parser.filepath = self.mainfile
        self.archive_parser.annotation_key = gpaw.GPW_KEY
        self.archive_parser.data_object = Simulation()

        self.mainfile_parser.convert(self.archive_parser)
        self.archive.data = self.archive_parser.data_object
        self.archive.workflow2 = self.mainfile_parser.build_workflow()

        self.mainfile_parser.close()
        self.archive_parser.close()


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
