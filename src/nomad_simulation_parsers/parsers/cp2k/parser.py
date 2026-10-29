from typing import Any

import numpy as np
from nomad.datamodel import EntryArchive
from nomad.parsing import MatchingParser
from nomad.units import ureg
from nomad_file_parser import ArchiveWriter
from nomad_file_parser.mapping_parser import MetainfoParser, TextParser
from nomad_simulations.schema_packages.general import Program, Simulation
from nomad_simulations.schema_packages.workflow.single_point import SinglePoint
from structlog.stdlib import BoundLogger

from nomad_simulation_parsers.schema_packages import cp2k

from .file_parser import CP2KOutParser


class CP2KOutputParser(TextParser):
    def get_version(self, source: dict[str, Any]) -> str | None:
        return dict(source.get('cp2k') or []).get('version_string:')

    def get_run_type(self, source: dict[str, Any]) -> str | None:
        return dict(source.get('global') or []).get('run_type')

    def get_systems(self, source: dict[str, Any]) -> list[dict[str, Any]]:
        quickstep = source.get('quickstep') or source.get('qs_dftb')
        if quickstep is None:
            self.logger.warning('No CP2K Quickstep output found')
            return []

        coordinates = quickstep.get('atomic_coordinates')
        if coordinates is None:
            return []

        system = {
            'positions': coordinates.positions,
            'atoms': [
                {'label': label, 'number': coordinates.atomic_numbers[label]}
                for label in coordinates.labels
            ],
        }
        lattice_vectors = source.get('lattice_vectors')
        if lattice_vectors is not None and np.shape(lattice_vectors) == (3, 3):
            system['lattice_vectors'] = np.asarray(lattice_vectors) * ureg.angstrom
        return [system]

    def get_outputs(self, source: dict[str, Any]) -> list[dict[str, Any]]:
        quickstep = source.get('quickstep') or source.get('qs_dftb')
        if quickstep is None:
            return []

        run_type = self.get_run_type(source)
        if run_type not in ('ENERGY', 'ENERGY_FORCE'):
            if run_type:
                self.logger.warning(
                    'CP2K run type is not yet mapped', run_type=run_type
                )
            return []

        single_point = quickstep.get('single_point')
        if single_point is None:
            return []

        output = {}
        energy = single_point.get('energy_total')
        if energy is not None:
            output['energy'] = energy
        forces = single_point.get('atom_forces')
        if forces is not None:
            output['forces'] = forces * ureg.hartree / ureg.bohr
        iterations = single_point.get('iteration') or []
        if iterations:
            output['scf_steps'] = {
                'energies_total': [step['energy_total'] for step in iterations],
                'delta_energies_total': [
                    abs(step['energy_change']) for step in iterations
                ],
            }
        return [output] if output else []


class CP2KArchiveWriter(ArchiveWriter):
    output_parser = CP2KOutputParser(text_parser=CP2KOutParser())
    archive_parser = MetainfoParser()

    def write_to_archive(self):
        self.output_parser.logger = self.logger
        self.archive_parser.logger = self.logger
        self.output_parser.filepath = self.mainfile
        self.archive_parser.annotation_key = cp2k.OUT_KEY
        self.archive_parser.data_object = Simulation(program=Program(name='CP2K'))
        self.output_parser.convert(self.archive_parser)
        self.archive.data = self.archive_parser.data_object
        if self.output_parser.get_run_type(self.output_parser.data) in (
            'ENERGY',
            'ENERGY_FORCE',
        ):
            self.archive.workflow2 = SinglePoint()


class CP2KParser(MatchingParser):
    archive_writer = CP2KArchiveWriter()

    def parse(
        self,
        mainfile: str,
        archive: EntryArchive,
        logger: BoundLogger,
        child_archives: dict[str, EntryArchive] | None = None,
    ):
        self.archive_writer.write(mainfile, archive, logger, child_archives)
