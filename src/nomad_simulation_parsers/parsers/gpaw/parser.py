from typing import Any

import numpy as np
from nomad.datamodel import EntryArchive
from nomad.parsing.parser import MatchingParser
from nomad.utils import get_logger
from nomad_file_parser import ArchiveWriter
from nomad_file_parser.mapping_parser import MappingParser, MetainfoParser
from nomad_simulations.schema_packages.general import Simulation
from nomad_simulations.schema_packages.workflow.general import EnergyConvergenceTarget
from nomad_simulations.schema_packages.workflow.single_point import (
    SinglePoint,
    SinglePointMethod,
)
from structlog.stdlib import BoundLogger

from nomad_simulation_parsers.parsers.utils.general import (
    calculate_band_gap_from_occupations,
)
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
            value = self.data_object.apply_unit(
                self.data_object.parser.get_parameter(f'energy_{key}'), 'energyunit'
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
        energyunit = self.data_object.apply_unit(1, 'energyunit').units
        lengthunit = self.data_object.apply_unit(1, 'lengthunit').units
        value = self.data_object.parser.get_array('atom_forces_free')
        if value is not None:
            forces['value'] = value * energyunit / lengthunit
        return forces

    def get_eigenvalues(self) -> list[dict[str, Any]]:
        reference_energy = self.get_reference_energy()
        eigenvalues = self.data_object.apply_unit(
            self.data_object.parser.get_array('eigenvalues'), 'energyunit'
        )
        occupations = self.data_object.parser.get_array('occupation')

        if eigenvalues is None or occupations is None:
            return []

        ndim_non_spin_polarized = 2  # shape: (n_kpoints, n_bands)
        n_spin_channels = 2

        # eigenvalues and occupations should have shape [n_spin, n_kpoints, n_bands]
        # or [n_kpoints, n_bands] for non-spin-polarized
        if eigenvalues.ndim == ndim_non_spin_polarized:
            # Non-spin-polarized: reshape to add spin dimension
            eigenvalues = eigenvalues[np.newaxis, :, :]
            occupations = occupations[np.newaxis, :, :]

        n_spin = eigenvalues.shape[0]
        n_bands = (
            eigenvalues.shape[2]
            if eigenvalues.ndim > ndim_non_spin_polarized
            else eigenvalues.shape[1]
        )

        data = []
        for spin_idx in range(n_spin):
            entry = dict(
                value=eigenvalues[spin_idx],  # 2D: [n_kpoints, n_bands]
                occupation=occupations[spin_idx],  # 2D: [n_kpoints, n_bands]
                n_levels=n_bands,
                highest_occupied=reference_energy,
            )
            # Only add spin_channel if there are multiple spins
            if n_spin == n_spin_channels:
                entry['spin_channel'] = spin_idx
            data.append(entry)

        return data

    def get_reference_energy(self):
        fermi_level = self.data_object.parser.get_parameter('fermilevel')
        if fermi_level is None:
            return None

        fermi_values = np.asarray(fermi_level, dtype=float).reshape(-1)
        if fermi_values.size == 0:
            return None

        return self.data_object.apply_unit(float(fermi_values[0]), 'energyunit')

    def get_band_structures(self) -> list[dict[str, Any]]:
        reference_energy = self.get_reference_energy()
        band_paths = self.data_object.parser.get_array('band_paths')

        band_structures = []
        for band_path in band_paths or []:
            eigenvalues = band_path.get('eigenvalues')
            if eigenvalues is None:
                continue
            band_structures.append(
                dict(
                    value=self.data_object.apply_unit(eigenvalues, 'energyunit'),
                    highest_occupied=reference_energy,
                )
            )
        return band_structures

    def get_band_gaps(self) -> list[dict[str, Any]]:
        """Calculate band gaps from eigenvalues using common utility."""
        band_gaps = []
        for spin_channel, eigenvalue_data in enumerate(self.get_eigenvalues()):
            energies = eigenvalue_data.get('value')
            occupations = eigenvalue_data.get('occupation')

            # Use common utility for band gap calculation (handles units automatically)
            gap_result = calculate_band_gap_from_occupations(
                energies, occupations, spin_channel=spin_channel
            )
            if gap_result is not None:
                band_gaps.append(gap_result)

        return band_gaps

    def get_scf_steps(self) -> dict[str, Any]:
        code_specific_quantities = {}
        converged = self.data_object.parser.get_parameter('converged')
        if converged is not None:
            code_specific_quantities['converged'] = bool(converged)

        energy_error = self.data_object.parser.get_parameter('energyerror')
        if energy_error is not None:
            code_specific_quantities['energyerror'] = float(energy_error)

        if code_specific_quantities:
            return {'code_specific_quantities': code_specific_quantities}
        return {}

    def build_workflow(self):
        workflow = SinglePoint()
        workflow.method = SinglePointMethod()
        energy_error = self.data_object.parser.get_parameter('energyerror')
        if energy_error is not None:
            workflow.method.convergence_targets = [
                EnergyConvergenceTarget(
                    threshold=self.data_object.apply_unit(energy_error, 'energyunit'),
                    threshold_type='absolute',
                )
            ]
        return workflow

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
