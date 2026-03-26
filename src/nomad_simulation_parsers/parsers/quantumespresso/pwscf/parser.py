from typing import Any

import numpy as np
from nomad.datamodel import EntryArchive
from nomad.units import ureg
from nomad.utils import get_logger
from nomad_simulations.schema_packages.workflow import (
    GeometryOptimization,
    MolecularDynamics,
    SinglePoint,
)
from nomad_simulations.schema_packages.workflow.general import (
    EnergyConvergenceTarget,
    ForceConvergenceTarget,
)
from nomad_simulations.schema_packages.workflow.geometry_optimization import (
    GeometryOptimizationMethod,
)
from nomad_simulations.schema_packages.workflow.single_point import SinglePointMethod

from nomad_simulation_parsers.parsers.quantumespresso.parser import (
    QuantumEspressoArchiveWriter,
)
from nomad_simulation_parsers.schema_packages.quantumespresso import pwscf

from ..parser import MainfileTextParser, MainfileXMLParser
from .file_parser import PWSCFFileParser

LOGGER = get_logger(__name__)


class PWSCFMainfileTextParser(MainfileTextParser):
    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER

    def get_force_contributions(self, source: dict[str, Any]) -> list[dict[str, Any]]:
        keys = ['dispersion']
        return [
            dict(name=key, value=source[f'forces_{key}'])
            for key in keys
            if source.get(f'forces_{key}') is not None
        ]

    def get_eigenvalues(self, source: dict[str, Any]) -> list[dict[str, Any]]:
        eigenvalues = source.get('band_energies')
        if eigenvalues is None:
            return []
        n_spin = self.get_n_spin_channels()
        n_eigs = len(eigenvalues[0])
        n_bands = np.size(eigenvalues) // int(n_spin * n_eigs)
        eigenvalues = np.reshape(eigenvalues, (n_spin, n_bands, n_eigs)) * ureg.eV
        results = [dict(eigenvalues=eig) for n, eig in enumerate(eigenvalues)]
        occupations = source.get('occupation_numbers')
        if occupations is not None:
            occupations = np.reshape(occupations, (n_spin, n_bands, n_eigs))
            for n, occ in enumerate(occupations):
                results[n]['occupations'] = occ
        return results

    def get_configurations(self, source: dict[str, Any]) -> list[dict[str, Any]]:
        methods = {
            'self_consistent': 'single_point',
            'bandstructure': 'single_point',
            'bfgs_geometry_optimization': 'geometry_optimization',
            'molecular_dynamics': 'molecular_dynamics',
            'langevin_dynamics': 'langevin_dynamics',
            'damped_dynamics': 'geometry_optimization',
            'vcs_wentzcovitch_damped_minimization': 'geometry_optimization',
        }

        configurations = []
        for key in methods:
            config = source.get(key)
            if config is None:
                continue
            configurations.append(config.get('self_consistent', config))
        return configurations

    def _resolve_scf_source(self, source: Any) -> Any:
        if isinstance(source, list):
            return source[-1] if source else None
        return source

    def get_scf_steps(self, source: Any) -> dict[str, Any]:  # noqa: PLR0912
        section = self._resolve_scf_source(source)
        if section is None or not hasattr(section, 'get'):
            return {}
        iterations = section.get('iteration', [])
        if not iterations:
            return {}

        energies_total = []
        delta_energies_total = []
        durations = []
        thresholds = []
        ddv_scf = []

        for iteration in iterations:
            energies = iteration.get('energies', {})
            energy_total = energies.get('energy_total')
            if energy_total is not None:
                energies_total.append(energy_total)
            energy_accuracy = energies.get('energy_total_accuracy_estimate')
            if energy_accuracy is not None:
                delta_energies_total.append(abs(energy_accuracy))
            duration = iteration.get('time')
            if duration is not None:
                durations.append(float(duration))
            threshold = iteration.get('threshold')
            if threshold is not None:
                thresholds.append(float(threshold))
            dscf = iteration.get('ddv_scf')
            if dscf is not None:
                ddv_scf.append(float(dscf))

        if not energies_total:
            return {}

        if not delta_energies_total and len(energies_total) > 1:
            for idx in range(1, len(energies_total)):
                delta_energies_total.append(
                    abs(energies_total[idx] - energies_total[idx - 1])
                )

        scf_steps = {'energies_total': energies_total}
        if delta_energies_total:
            scf_steps['delta_energies_total'] = delta_energies_total
        if len(durations) == len(energies_total):
            scf_steps['durations'] = durations
        if thresholds or ddv_scf:
            scf_steps['code_specific_quantities'] = {
                'threshold': thresholds,
                'ddv_scf': ddv_scf,
            }
        return scf_steps

    def build_workflow(self):
        if self.data.get('bfgs_geometry_optimization') is not None:
            workflow = GeometryOptimization()
            workflow.method = GeometryOptimizationMethod()
            header = self.data.get('header', {})
            force_thr = header.get('forc_conv_thr')
            energy_thr = header.get('etot_conv_thr')
            targets = []
            if force_thr is not None:
                targets.append(
                    ForceConvergenceTarget(
                        threshold=float(force_thr) * ureg.rydberg / ureg.bohr,
                        threshold_type='maximum',
                    )
                )
            if energy_thr is not None:
                targets.append(
                    EnergyConvergenceTarget(
                        threshold=float(energy_thr) * ureg.rydberg,
                        threshold_type='absolute',
                    )
                )
            if targets:
                workflow.method.convergence_targets = targets
            scf_thr = header.get('scf_threshold_energy_change')
            if scf_thr is not None:
                workflow.method.single_point_convergence_targets = [
                    EnergyConvergenceTarget(
                        threshold=scf_thr,
                        threshold_type='absolute',
                    )
                ]
            return workflow

        if self.data.get('molecular_dynamics') is not None:
            return MolecularDynamics()

        workflow = SinglePoint()
        workflow.method = SinglePointMethod()
        scf_thr = self.data.get('header', {}).get('scf_threshold_energy_change')
        if scf_thr is not None:
            workflow.method.convergence_targets = [
                EnergyConvergenceTarget(
                    threshold=scf_thr,
                    threshold_type='absolute',
                )
            ]
        return workflow


class PWSCFMainfileXMLParser(MainfileXMLParser):
    def get_configurations(self, source: dict[str, Any]):
        keys = ['input', 'output']
        return [source[key] for key in keys if source.get(key) is not None]

    def get_scf_steps(self, source: dict[str, Any]) -> dict[str, Any]:
        steps = self.data.get('step', [])
        if not isinstance(steps, list) or not steps:
            return {}

        energies_total = []
        delta_energies_total = []
        scf_steps_count = []
        for step in steps:
            total_energy = step.get('total_energy', {}).get('etot')
            if total_energy is not None:
                energies_total.append(self.apply_unit(total_energy, name='energy'))
            scf_conv = step.get('scf_conv', {})
            scf_error = scf_conv.get('scf_error')
            if scf_error is not None:
                delta_energies_total.append(self.apply_unit(scf_error, name='energy'))
            n_scf_steps = scf_conv.get('n_scf_steps')
            if n_scf_steps is not None:
                scf_steps_count.append(int(n_scf_steps))

        if not energies_total:
            return {}

        if not delta_energies_total and len(energies_total) > 1:
            for idx in range(1, len(energies_total)):
                delta_energies_total.append(
                    abs(energies_total[idx] - energies_total[idx - 1])
                )

        scf_steps = {'energies_total': energies_total}
        if delta_energies_total:
            scf_steps['delta_energies_total'] = delta_energies_total
        if scf_steps_count:
            scf_steps['code_specific_quantities'] = {'n_scf_steps': scf_steps_count}
        return scf_steps

    def build_workflow(self):
        control = self.data.get('input', {}).get('control_variables', {})
        calc = control.get('calculation')
        relax_labels = {'relax', 'vc-relax'}
        md_labels = {'md', 'vc-md'}

        workflow = None
        if calc in relax_labels:
            workflow = GeometryOptimization()
            workflow.method = GeometryOptimizationMethod()
            targets = []
            force_thr = control.get('forc_conv_thr')
            energy_thr = control.get('etot_conv_thr')
            if force_thr is not None:
                targets.append(
                    ForceConvergenceTarget(
                        threshold=float(force_thr) * ureg.rydberg / ureg.bohr,
                        threshold_type='maximum',
                    )
                )
            if energy_thr is not None:
                targets.append(
                    EnergyConvergenceTarget(
                        threshold=float(energy_thr) * ureg.rydberg,
                        threshold_type='absolute',
                    )
                )
            if targets:
                workflow.method.convergence_targets = targets
        elif calc in md_labels:
            workflow = MolecularDynamics()
        else:
            workflow = SinglePoint()
            workflow.method = SinglePointMethod()

        scf_thr = self.data.get('input', {}).get('electron_control', {}).get('conv_thr')
        if scf_thr is not None:
            sp_target = EnergyConvergenceTarget(
                threshold=float(scf_thr) * ureg.rydberg,
                threshold_type='absolute',
            )
            if isinstance(workflow, GeometryOptimization):
                workflow.method.single_point_convergence_targets = [sp_target]
            elif workflow.method is not None:
                workflow.method.convergence_targets = [sp_target]
        return workflow


class PWSCFArchiveWriter(QuantumEspressoArchiveWriter):
    schema = pwscf
    _text_parser = PWSCFMainfileTextParser(text_parser=PWSCFFileParser())
    _xml_parser = PWSCFMainfileXMLParser()

    def parse_program(self, archive: EntryArchive, index: int) -> None:
        super().parse_program(archive, index)
        archive.workflow2 = self.mainfile_parser.build_workflow()
