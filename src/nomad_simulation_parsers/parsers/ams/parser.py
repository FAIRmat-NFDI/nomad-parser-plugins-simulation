import os
from typing import Any

import numpy as np
from nomad.datamodel import EntryArchive
from nomad.parsing.parser import MatchingParser
from nomad.units import ureg
from nomad.utils import get_logger
from nomad_file_parser import ArchiveWriter
from nomad_file_parser.mapping_parser import MetainfoParser, TextParser
from nomad_simulations.schema_packages.general import Program, Simulation
from nomad_simulations.schema_packages.workflow.general import (
    EnergyConvergenceTarget,
    ForceConvergenceTarget,
)
from nomad_simulations.schema_packages.workflow.geometry_optimization import (
    GeometryOptimization,
    GeometryOptimizationMethod,
)
from nomad_simulations.schema_packages.workflow.molecular_dynamics import (
    MolecularDynamics,
)
from nomad_simulations.schema_packages.workflow.single_point import (
    SinglePoint,
    SinglePointMethod,
)
from structlog.stdlib import BoundLogger

from nomad_simulation_parsers.parsers.utils.general import (
    calculate_band_gap_from_occupations,
    link_outputs_to_model_systems,
    search_files,
)
from nomad_simulation_parsers.schema_packages import ams

from .file_parser import OutParser
from .file_parser import RKFParser as RKFTextParser

LOGGER = get_logger(__name__)
MIN_TUPLE_FIELDS = 3


class MainfileParser(TextParser):
    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER

    def get_xc_functionals(self, source: dict[str, Any]) -> list[str]:
        xc_functionals = []
        for xc_type in ['LDA', 'GGA', 'MGGA']:
            functionals = source.get(xc_type, '').split()
            kind = ['XC'] if len(functionals) == 1 else ['X', 'C']
            for n, functional in enumerate(functionals):
                xc_functionals.append(
                    f'{xc_type}_{kind[n]}_{functional.rstrip("x").rstrip("c").upper()}'
                )
        return xc_functionals

    def get_periodic_boundary_conditions(
        self, source: dict[str, Any] | Any
    ) -> list[bool] | None:
        non_periodic = [False, False, False]
        lattice_vectors = (
            source.get('lattice_vectors') if isinstance(source, dict) else source
        )
        if lattice_vectors is None:
            return non_periodic
        if hasattr(lattice_vectors, 'magnitude'):
            lattice_vectors = lattice_vectors.magnitude

        vectors = np.asarray(lattice_vectors)
        if vectors.size == 0:
            return non_periodic

        if vectors.ndim == 1:
            # Accept flattened lattice payloads.
            n_vectors = min(vectors.shape[0] // 3, 3)
        else:
            n_vectors = min(vectors.shape[0], 3)

        return [idx < n_vectors for idx in range(3)]

    def get_contributions(self, source: dict[str, Any]) -> list[dict[str, Any]]:
        return [
            dict(name=key, value=val) for key, val in source.items() if val is not None
        ]

    def get_eigenvalues(
        self, source: dict[str, Any] | list[np.ndarray]
    ) -> list[dict[str, np.ndarray]]:
        if source is None:
            return []

        if isinstance(source, dict):
            energies = source.get('energies', [])
            occupations = source.get('occupations', [])
        else:
            if not isinstance(source, tuple | list) or len(source) < MIN_TUPLE_FIELDS:
                return []
            energies = []
            occupations = source[2]

        nspin = max(len(energies), len(occupations))
        eigenvalues = [dict() for _ in range(nspin)]
        for n, energy in enumerate(energies):
            eigenvalues[n]['eigenvalues'] = energy
        for n, occupations in enumerate(occupations):
            eigenvalues[n]['occupations'] = occupations
        return [eig for eig in eigenvalues if eig]

<<<<<<< HEAD
    def get_band_gaps(self, source: Any) -> list[dict[str, Any]]:  # noqa: PLR0912
        if source is None:
            return []

        if hasattr(source, 'get'):
            value = source.get('value')
            if value is None:
                homo = source.get('energy_highest_occupied')
                lumo = source.get('energy_lowest_unoccupied')
                if homo is not None and lumo is not None:
                    value = lumo - homo
            if value is None:
                return []

            if hasattr(value, 'magnitude'):
                if value.magnitude < 0:
                    value = 0 * value.units
            else:
                value = max(0.0, value)

            band_gap = {'value': value}
            spin_channel = source.get('spin_channel')
            if spin_channel is not None:
                band_gap['spin_channel'] = spin_channel
            return [band_gap]

        # Fallback for parsed tuple/list payload used by AMS band-energy ranges.
        # Expected form: (energies, ..., occupations) with spin channels separated.
        if not isinstance(source, tuple | list) or len(source) < MIN_TUPLE_FIELDS:
            return []

        energies_channels = source[0] if isinstance(source[0], list) else []
        occupations_channels = source[2] if isinstance(source[2], list) else []
        if not energies_channels or not occupations_channels:
            return []

        band_gaps = []
        for spin_channel, energies_entry in enumerate(energies_channels):
            if spin_channel >= len(occupations_channels):
                continue
            occupation_data = occupations_channels[spin_channel]
            if not isinstance(occupation_data, list) or not occupation_data:
                continue
            occupations = np.asarray(occupation_data[0], dtype=float)
            energies = np.asarray(energies_entry, dtype=float)
            if (
                energies.size == 0
                or occupations.size == 0
                or energies.shape != occupations.shape
            ):
                continue

            # Use common utility for band gap calculation
            gap = calculate_band_gap_from_occupations(energies, occupations)
            if gap is not None:
                # Override spin channel from utility with AMS channel index
                gap['spin_channel'] = spin_channel
                band_gaps.append(gap)
        return band_gaps

    def get_dos(self, source: dict[str, Any]) -> list[dict[str, Any]]:
        if source is None:
            return []

        if not isinstance(source, dict):
            return []

        dos = source.get('dos')
        if dos is None:
            return []

        dos_dimensions = 2
        minimum_dos_columns = 2
        dos = np.asarray(dos)
        if dos.ndim != dos_dimensions or dos.shape[1] < minimum_dos_columns:
            return []

        energies = dos[:, 0]
        return [dict(energies=energies, value=values) for values in dos[:, 1:].T]

=======
>>>>>>> 2bebd7c (Schema update convergence targets (#150))
    def get_scf_steps(self, source: dict[str, Any]) -> dict[str, Any]:
        self_consistency = source.get('self_consistency', {})
        energy_change = self_consistency.get('energy_change')
        if energy_change is None:
            return {}

        delta_energies_total = [abs(value) for value in energy_change]
        scf_steps = {'delta_energies_total': delta_energies_total}

        scf_options = source.get('scf_options')
        code_specific_quantities = {}
<<<<<<< HEAD
        if scf_options is not None:
            if not isinstance(scf_options, dict):
                return scf_steps

=======
        if hasattr(scf_options, 'get'):
>>>>>>> 2bebd7c (Schema update convergence targets (#150))
            n_scf_steps_max = scf_options.get('x_ams_ncyclx')
            convrg = scf_options.get('x_ams_convrg')
            if n_scf_steps_max is not None:
                code_specific_quantities['n_scf_steps_max'] = int(n_scf_steps_max)
            if convrg is not None:
                code_specific_quantities['convrg'] = float(convrg)

        if code_specific_quantities:
            scf_steps['code_specific_quantities'] = code_specific_quantities

        return scf_steps

    def _get_scf_energy_threshold(self, source: dict[str, Any]):
        scf_options = source.get('scf_options')
<<<<<<< HEAD
        if scf_options is None:
            return None
        convrg = (
            scf_options.get('x_ams_convrg') if hasattr(scf_options, 'get') else None
        )
=======
        if not hasattr(scf_options, 'get'):
            return None
        convrg = scf_options.get('x_ams_convrg')
>>>>>>> 2bebd7c (Schema update convergence targets (#150))
        if convrg is None:
            return None
        return float(convrg) * ureg.hartree

    def build_workflow(self, source: dict[str, Any]):
        if (geometry := source.get('geometry_optimization')) is not None:
            workflow = GeometryOptimization()
            workflow.method = GeometryOptimizationMethod()
            targets = []
            force_thr = geometry.get('convergence_tolerance_force_maximum')
            if force_thr is not None:
                targets.append(
                    ForceConvergenceTarget(
                        threshold=force_thr,
                        threshold_type='maximum',
                    )
                )
            energy_thr = geometry.get('convergence_tolerance_energy_difference')
            if energy_thr is not None:
                targets.append(
                    EnergyConvergenceTarget(
                        threshold=energy_thr,
                        threshold_type='absolute',
                    )
                )
            if targets:
                workflow.method.convergence_targets = targets

            scf_threshold = self._get_scf_energy_threshold(geometry)
            if scf_threshold is not None:
                workflow.method.single_point_convergence_targets = [
                    EnergyConvergenceTarget(
                        threshold=scf_threshold,
                        threshold_type='absolute',
                    )
                ]
            return workflow

        if source.get('molecular_dynamics') is not None:
            return MolecularDynamics()

        workflow = SinglePoint()
        workflow.method = SinglePointMethod()
        single_point = source.get('single_point', source)
        scf_threshold = self._get_scf_energy_threshold(single_point)
        if scf_threshold is not None:
            workflow.method.convergence_targets = [
                EnergyConvergenceTarget(
                    threshold=scf_threshold,
                    threshold_type='absolute',
                )
            ]
        return workflow


class RKFParser(MainfileParser):
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
    rkf_parser = RKFParser(text_parser=RKFTextParser())

    def write_to_archive(self):
        self.metainfo_parser.annotation_key = ams.OUT_KEY
        self.archive.data = Simulation(program=Program(name='AMS'))
        self.metainfo_parser.data_object = self.archive.data

        rkf_files = search_files('ams.rkf', os.path.dirname(self.mainfile))
        self.parser = self.mainfile_parser
        self.parser.filepath = self.mainfile
        if rkf_files:
            if len(rkf_files) > 1:
                self.logger.warning('Multiple ams.rkf files found.')
            self.parser = self.rkf_parser
            self.parser.filepath = rkf_files[0]
            self.parser.data_object.parse()

        self.parser.convert(self.metainfo_parser)
<<<<<<< HEAD
        link_outputs_to_model_systems(self.archive.data)

=======
>>>>>>> 2bebd7c (Schema update convergence targets (#150))
        self.archive.workflow2 = self.parser.build_workflow(self.parser.data)


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
