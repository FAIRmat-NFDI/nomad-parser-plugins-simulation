import os
from typing import Any

import numpy as np
from nomad.datamodel import EntryArchive
from nomad.units import ureg
from nomad.utils import get_logger
from nomad_file_parser.mapping_parser import TextParser
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
from nomad_simulation_parsers.schema_packages.quantumespresso import common, pwscf

from ..parser import MainfileTextParser, MainfileXMLParser
from .file_parser import PWSCFDOSTextParser, PWSCFFileParser

LOGGER = get_logger(__name__)
MIN_DOS_COLUMNS = 2
DOS_ARRAY_NDIM = 2


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
            dict(
                eigenvalues=eig,
                n_levels=eig.shape[-1],
                spin_channel=n if n_spin > 1 else None,
            )
            for n, eig in enumerate(eigenvalues)
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

    _configuration_methods = {
        'self_consistent': 'single_point',
        'bandstructure': 'single_point',
        'bfgs_geometry_optimization': 'geometry_optimization',
        'molecular_dynamics': 'molecular_dynamics',
        'langevin_dynamics': 'langevin_dynamics',
        'damped_dynamics': 'geometry_optimization',
        'vcs_wentzcovitch_damped_minimization': 'geometry_optimization',
    }

    def get_configurations(self, source: dict[str, Any]) -> list[dict[str, Any]]:
        methods = self._configuration_methods

        configurations = []
        for key in methods:
            config = source.get(key)
            if config is None:
                continue
            sc_config = config.get('self_consistent', config)
            sec_config = sc_config if isinstance(sc_config, list) else [sc_config]
            configurations.extend(sec_config)
        return configurations

    def get_configuration_forces(self, source: dict[str, Any]) -> list[list[Any]]:
        """Per-configuration force series, index-aligned with get_configurations.

        Multi-step runs (geometry optimization, molecular dynamics) carry one
        force array per SCF step; single-point configurations yield an empty
        series since their forces are covered by the step-resolved outputs.
        """
        forces = []
        for key in self._configuration_methods:
            config = source.get(key)
            if config is None:
                continue
            sc_config = config.get('self_consistent', config)
            if isinstance(sc_config, list):
                if not sc_config:
                    continue
                forces.append(
                    [
                        step.get('forces')
                        for step in sc_config
                        if hasattr(step, 'get') and step.get('forces') is not None
                    ]
                )
            else:
                forces.append([])
        return forces

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

    def get_reference_energy(self, source: dict[str, Any]) -> Any | None:
        section = self._resolve_scf_source(source)
        if section is None or not hasattr(section, 'get'):
            return None

        section = section.get('self_consistent', section)
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

    # DOS payload cache: `get_dos` runs once per configuration, but the .dos
    # sidecar files only need to be read once per mainfile.
    _dos_cache_key: str | None = None
    _cached_dos_payload: dict[str, Any] | None = None

    def get_dos(self, source: dict[str, Any]) -> list[dict[str, Any]]:
        mainfile = getattr(self, 'filepath', None)
        if not isinstance(mainfile, str) or not mainfile:
            return []
        if self._dos_cache_key != mainfile:
            self._dos_cache_key = mainfile
            self._cached_dos_payload = None
            dos_files = search_files(pattern='*.dos', basedir=os.path.dirname(mainfile))
            for dos_file in dos_files:
                try:
                    data = np.loadtxt(dos_file, comments='#')
                except Exception:
                    continue
                if data is None:
                    continue
                if data.ndim == 1 and data.size >= MIN_DOS_COLUMNS:
                    data = data.reshape(1, -1)
                if data.ndim == DOS_ARRAY_NDIM and data.shape[1] >= MIN_DOS_COLUMNS:
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


class DOSParser(TextParser):
    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER

    # TODO: fix to prevent creation of contributions sections.
    # clarify the purpose of ElectronicDensityOfStates.contributions
    # and projected_dos
    def get_dos_contributions(self) -> list[dict[str, Any]]:
        return []

    def load_file(self):
        self.text_parser.mainfile = None
        for dos_file in search_files(
            pattern='*.dos', basedir=os.path.dirname(self.filepath)
        ):
            self.text_parser.mainfile = dos_file
            if self.text_parser.energies is not None:
                break
        return self.text_parser


class PWSCFArchiveWriter(QuantumEspressoArchiveWriter):
    schema = pwscf
    _text_parser = PWSCFMainfileTextParser(text_parser=PWSCFFileParser())
    _xml_parser = PWSCFMainfileXMLParser()
    dos_parser = DOSParser(text_parser=PWSCFDOSTextParser())

    def parse_program(self, archive: EntryArchive, index: int) -> None:  # noqa: PLR0912, PLR0915
        super().parse_program(archive, index)
        archive.workflow2 = self.mainfile_parser.build_workflow()

        if not archive.data:
            return

        # parse dos
        self.dos_parser.filepath = self._mainfile_parser.filepath
        self.simulation_parser.annotation_key = common.DOS_KEY
        self.dos_parser.convert(self.simulation_parser)
        if archive.data.outputs and archive.data.outputs[-1].electronic_dos:
            # parse reference energy from out_file
            self.simulation_parser.data_object = archive.data.outputs[
                -1
            ].electronic_dos[-1]
            self.simulation_parser.annotation_key = common.DOS_OUT_KEY
            self.mainfile_parser.convert(self.simulation_parser)
            # reset to archive
            self.simulation_parser.data_object = archive
