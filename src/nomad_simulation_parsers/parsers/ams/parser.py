import os
from typing import Any

import numpy as np
from nomad.datamodel import EntryArchive
from nomad.parsing.file_parser import ArchiveWriter
from nomad.parsing.file_parser.mapping_parser import MetainfoParser, TextParser
from nomad.parsing.parser import MatchingParser
from nomad.units import ureg
from nomad.utils import get_logger
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

from nomad_simulation_parsers.parsers.utils.general import search_files
from nomad_simulation_parsers.schema_packages import ams

from .file_parser import OutParser
from .file_parser import RKFParser as RKFTextParser

LOGGER = get_logger(__name__)


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

    def get_contributions(self, source: dict[str, Any]) -> list[dict[str, Any]]:
        return [
            dict(name=key, value=val) for key, val in source.items() if val is not None
        ]

    def get_eigenvalues(
        self, source: dict[str, Any] | list[np.ndarray]
    ) -> list[dict[str, np.ndarray]]:
        energies = source.get('energies', []) if isinstance(source, dict) else []
        occupations = (
            source.get('occupations', []) if isinstance(source, dict) else source[2]
        )
        nspin = max(len(energies), len(occupations))
        eigenvalues = [dict() for _ in range(nspin)]
        for n, energy in enumerate(energies):
            eigenvalues[n]['eigenvalues'] = energy
        for n, occupations in enumerate(occupations):
            eigenvalues[n]['occupations'] = occupations
        return [eig for eig in eigenvalues if eig]

    def get_band_gaps(self, source: dict[str, Any]) -> list[dict[str, Any]]:
        if not hasattr(source, 'get'):
            return []

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

    def get_dos(self, source: dict[str, Any]) -> list[dict[str, Any]]:
        if not hasattr(source, 'get'):
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

    def get_scf_steps(self, source: dict[str, Any]) -> dict[str, Any]:
        self_consistency = source.get('self_consistency', {})
        energy_change = self_consistency.get('energy_change')
        if energy_change is None:
            return {}

        delta_energies_total = [abs(value) for value in energy_change]
        scf_steps = {'delta_energies_total': delta_energies_total}

        scf_options = source.get('scf_options')
        code_specific_quantities = {}
        if hasattr(scf_options, 'get'):
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
        if not hasattr(scf_options, 'get'):
            return None
        convrg = scf_options.get('x_ams_convrg')
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
