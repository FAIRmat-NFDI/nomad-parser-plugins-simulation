from typing import Any
import os

import numpy as np
from nomad.datamodel import EntryArchive
from nomad.units import ureg
from nomad.utils import get_logger
from nomad_simulations.schema_packages import outputs as simulation_outputs
from nomad_simulations.schema_packages import variables as simulation_variables
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
from nomad_simulation_parsers.parsers.utils.general import search_files
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
        section = self._resolve_scf_source(source)
        if section is None or not hasattr(section, 'get'):
            return []

        eigenvalues = section.get('band_energies')
        if eigenvalues is None:
            return []
        n_spin = self.get_n_spin_channels()
        n_eigs = len(eigenvalues[0])
        n_bands = np.size(eigenvalues) // int(n_spin * n_eigs)
        eigenvalues = np.reshape(eigenvalues, (n_spin, n_bands, n_eigs)) * ureg.eV
        results = [
            dict(eigenvalues=eig, n_levels=eig.shape[-1]) for n, eig in enumerate(eigenvalues)
        ]
        occupations = section.get('occupation_numbers')
        if occupations is not None:
            occupations = np.reshape(occupations, (n_spin, n_bands, n_eigs))
            for n, occ in enumerate(occupations):
                results[n]['occupations'] = occ
        else:
            fermi_energy = section.get('fermi_energy')
            if fermi_energy is not None:
                fermi_array = np.asarray(fermi_energy, dtype=float).reshape(-1)
                if fermi_array.size == 1:
                    fermi_values = [fermi_array[0]] * n_spin
                else:
                    fermi_values = [
                        fermi_array[min(i, fermi_array.size - 1)] for i in range(n_spin)
                    ]
                for n, eig in enumerate(eigenvalues):
                    results[n]['occupations'] = (
                        eig.magnitude <= fermi_values[n]
                    ).astype(float)
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
        header = source.get('header', {}) if isinstance(source, dict) else {}
        for key in methods:
            config = source.get(key)
            if config is None:
                continue
            sc_config = config.get('self_consistent', config)
            if isinstance(sc_config, list):
                if not sc_config:
                    continue
                sec = sc_config[-1]
            else:
                sec = sc_config

            if isinstance(sec, dict):
                if sec.get('simulation_cell') is None and header.get('simulation_cell') is not None:
                    sec = sec.copy()
                    sec['simulation_cell'] = header.get('simulation_cell')
                if sec.get('labels_positions') is None and header.get('labels_positions') is not None:
                    if sec is sc_config:
                        sec = sec.copy()
                    sec['labels_positions'] = header.get('labels_positions')

            payload_target = sec if isinstance(sec, dict) else getattr(sec, 'data', None)
            if isinstance(payload_target, dict):
                payload_target['electronic_eigenvalues'] = self.get_eigenvalues(sec)
                payload_target['electronic_band_structures'] = self.get_band_structures(
                    sec
                )
                payload_target['electronic_dos'] = self.get_dos(sec)

            configurations.append(sec)
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

    def get_reference_energy(self, source: dict[str, Any]):
        section = self._resolve_scf_source(source)
        if section is None or not hasattr(section, 'get'):
            return None

        homo_lumo = section.get('homo_lumo')
        if homo_lumo is not None:
            homo_vals = np.asarray(homo_lumo, dtype=float).reshape(-1)
            if homo_vals.size > 0:
                return homo_vals[0] * ureg.eV

        fermi = section.get('fermi_energy')
        if fermi is not None:
            fermi_vals = np.asarray(fermi, dtype=float).reshape(-1)
            if fermi_vals.size > 0:
                return fermi_vals[-1] * ureg.eV

        return None

    def get_band_structures(self, source: dict[str, Any]) -> list[dict[str, Any]]:
        reference_energy = self.get_reference_energy(source)
        band_structures = []
        for eigenvalues in self.get_eigenvalues(source):
            band_structures.append(
                dict(
                    value=eigenvalues.get('eigenvalues'),
                    occupation=eigenvalues.get('occupations'),
                    n_levels=eigenvalues.get('n_levels'),
                    spin_channel=eigenvalues.get('spin_channel'),
                    highest_occupied=reference_energy,
                )
            )
        return band_structures

    def get_dos(self, source: dict[str, Any]) -> list[dict[str, Any]]:
        if not hasattr(self, '_cached_dos_payload'):
            self._cached_dos_payload = None
            mainfile = getattr(self, 'filepath', None)
            if not isinstance(mainfile, str) or not mainfile:
                return []
            dos_files = search_files(
                pattern='*.dos', basedir=os.path.dirname(mainfile)
            )
            for dos_file in dos_files:
                try:
                    data = np.loadtxt(dos_file, comments='#')
                except Exception:
                    continue
                if data is None:
                    continue
                if data.ndim == 1 and data.size >= 2:
                    data = data.reshape(1, -1)
                if data.ndim == 2 and data.shape[1] >= 2:
                    energies = data[:, 0] * ureg.eV
                    values = np.abs(data[:, 1]) / ureg.eV
                    self._cached_dos_payload = dict(
                        value=values,
                        energies=dict(points=energies),
                    )
                    break

        if self._cached_dos_payload is None:
            return []

        result = dict(self._cached_dos_payload)
        reference_energy = self.get_reference_energy(source)
        if reference_energy is not None:
            result['energies_origin'] = reference_energy
        return [result]

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

        if not archive.data or not archive.data.outputs:
            return

        # TODO(mapping-migration): remove this parser-side electronic fallback once
        # PWSCF fixtures are fully covered by mapping-driven population only.
        # Mapping attempts tried in this iteration:
        # 1) Outputs mappings via ('get_eigenvalues'|'get_band_structures'|'get_dos', ['.@'])
        # 2) Precomputing payload keys in get_configurations + direct '.electronic_*' mappings
        # Both attempts left electronic sections empty on current fixtures
        # (tests/parsers/test_quantumespresso_parser.py), so we keep this as a
        # temporary compatibility fallback.
        if hasattr(self.mainfile_parser, 'get_eigenvalues'):
            configurations = self.mainfile_parser.get_configurations(
                self.mainfile_parser.data
            )
            for i, output in enumerate(archive.data.outputs):
                if i >= len(configurations):
                    break
                if output.electronic_eigenvalues:
                    continue

                try:
                    eigenvalues = self.mainfile_parser.get_eigenvalues(configurations[i])
                except Exception:
                    continue
                if not eigenvalues:
                    continue

                output.electronic_eigenvalues = [
                    simulation_outputs.ElectronicEigenvalues(
                        value=entry.get('eigenvalues'),
                        occupation=entry.get('occupations'),
                        n_levels=entry.get('n_levels'),
                        spin_channel=entry.get('spin_channel'),
                    )
                    for entry in eigenvalues
                ]

            for i, output in enumerate(archive.data.outputs):
                if output.electronic_band_structures:
                    continue

                if not output.electronic_eigenvalues:
                    continue

                config = configurations[i] if i < len(configurations) else None
                reference_energy = (
                    self.mainfile_parser.get_reference_energy(config)
                    if config is not None
                    else None
                )

                band_structures = []
                for eigenvalues in output.electronic_eigenvalues:
                    band_structure = simulation_outputs.ElectronicBandStructure(
                        value=eigenvalues.value,
                        occupation=eigenvalues.occupation,
                        n_levels=eigenvalues.n_levels,
                        spin_channel=eigenvalues.spin_channel,
                    )

                    if reference_energy is not None:
                        band_structure.highest_occupied = reference_energy

                    band_structures.append(band_structure)

                if band_structures:
                    output.electronic_band_structures = band_structures

                if output.electronic_dos and reference_energy is not None:
                    for dos in output.electronic_dos:
                        if dos.energies_origin is None:
                            dos.energies_origin = reference_energy

        dos_files = search_files(pattern='*.dos', basedir=os.path.dirname(self.mainfile))
        if not dos_files:
            return

        dos_data = None
        for dos_file in dos_files:
            try:
                data = np.loadtxt(dos_file, comments='#')
            except Exception:
                continue
            if data is None:
                continue
            if data.ndim == 1 and data.size >= 2:
                data = data.reshape(1, -1)
            if data.ndim == 2 and data.shape[1] >= 2:
                dos_data = data
                break

        if dos_data is None:
            return

        energies = dos_data[:, 0] * ureg.eV
        values = np.abs(dos_data[:, 1]) / ureg.eV
        for output in archive.data.outputs:
            if output.electronic_dos:
                continue
            output.electronic_dos = [
                simulation_outputs.ElectronicDensityOfStates(
                    value=values,
                    energies=simulation_variables.Energy2(points=energies),
                )
            ]

        if hasattr(self.mainfile_parser, 'get_configurations'):
            configurations = self.mainfile_parser.get_configurations(
                self.mainfile_parser.data
            )
            for i, output in enumerate(archive.data.outputs):
                if i >= len(configurations) or not output.electronic_dos:
                    continue
                reference_energy = self.mainfile_parser.get_reference_energy(
                    configurations[i]
                )
                if reference_energy is None:
                    continue
                for dos in output.electronic_dos:
                    if dos.energies_origin is None:
                        dos.energies_origin = reference_energy
