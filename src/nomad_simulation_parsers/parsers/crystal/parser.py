import datetime
import os
import re
from typing import Any

import numpy as np
import pint
from ase.data import chemical_symbols
from nomad import atomutils
from nomad.datamodel import EntryArchive
from nomad.parsing import MatchingParser
from nomad.parsing.file_parser import ArchiveWriter
from nomad.parsing.file_parser.mapping_parser import MetainfoParser, TextParser
from nomad.units import ureg
from nomad_simulations.schema_packages.general import Program, Simulation
from nomad_simulations.schema_packages.workflow.general import EnergyConvergenceTarget
from nomad_simulations.schema_packages.workflow.geometry_optimization import (
    GeometryOptimization,
    GeometryOptimizationMethod,
)
from nomad_simulations.schema_packages.workflow.single_point import (
    SinglePoint,
    SinglePointMethod,
)
from structlog.stdlib import BoundLogger

from nomad_simulation_parsers.schema_packages import crystal

from .file_parser import F25Parser, OutputParser


class CrystalOutputParser(TextParser):
    libxc_map = {
        'PBEXC': ['GGA_C_PBE', 'GGA_X_PBE'],
        'PBE0': ['HYB_GGA_XC_PBEH'],
        'B3LYP': ['HYB_GGA_XC_B3LYP'],
        'HSE06': ['HYB_GGA_XC_HSE06'],
        'M06': ['MGGA_C_M06', 'HYB_MGGA_X_M06'],
        'M05-2X': ['HYB_MGGA_XC_M05_2X'],
        'LC-WPBE': ['HYB_GGA_XC_LRC_WPBE'],
        'PBE': ['GGA_X_PBE', 'GGA_C_PBE'],
        'PBESOL': ['GGA_X_PBE_SOL', 'GGA_C_PBE_SOL'],
        'BECKE': ['GGA_X_B88'],
        'LDA': ['LDA_X'],
        'PWGGA': ['GGA_X_PW91', 'GGA_C_PW91'],
        'PZ': ['LDA_C_PZ'],
        'WFN': ['LDA_C_VWN'],
    }

    @property
    def logger(self):
        pass

    def to_unix_time(self, value: str) -> float | None:
        """Transforms the Crystal-specific float notation into a floating point
        number.
        """
        if value is None:
            return None

        value = value.strip()
        date_time_obj = datetime.datetime.strptime(value, '%d %m %Y TIME %H:%M:%S.%f')
        return date_time_obj.timestamp()

    def get_lattice_vectors(self, source: dict[str, Any]) -> pint.Quantity:
        lattice_vectors = source.get(
            'lattice_parameters', source.get('lattice_vectors_restart')
        )
        if lattice_vectors.shape == (6,):
            lattice_vectors = atomutils.cellpar_to_cell(lattice_vectors, degrees=True)
        return lattice_vectors * ureg.angstrom

    def get_positions(self, source: dict[str, Any]) -> pint.Quantity:
        labels_positions = source.get(
            'labels_positions_nanotube', source.get('labels_positions')
        )

        if labels_positions is None:
            labels_positions = self.data.get('labels_positions_restart')

        if labels_positions is None:
            return labels_positions

        positions = labels_positions[:, 4:7].astype(np.float64)

        dimensionality = self.data.get('dimensionality')
        lattice_vectors = self.get_lattice_vectors(source)

        scaled_pos = np.zeros((len(positions), 3), dtype=np.float64)
        scaled_pos[:, :dimensionality] = positions[:, :dimensionality]
        if lattice_vectors is not None:
            cart_pos = atomutils.to_cartesian(scaled_pos, lattice_vectors.magnitude)
            cart_pos[:, dimensionality:] = positions[:, dimensionality:]
        else:
            cart_pos = scaled_pos

        return cart_pos * ureg.angstrom

    def get_atoms(self, source: dict[str, Any]) -> list[str]:
        labels_positions = source.get(
            'labels_positions_nanotube', source.get('labels_positions')
        )
        labels = []
        numbers = []
        if labels_positions is not None:
            labels = labels_positions[:, 3]
            numbers = labels_positions[:, 2]
        else:
            labels_positions = self.data.get('labels_positions_restart')
            if labels_positions is not None:
                labels = labels_positions[:, 2]
                numbers = labels_positions[:, 1]

        def normalize_label(label: str) -> str:
            norm = re.match(r'([a-z][a-z]?).*', label.lower())
            # unknown specie
            # TODO not possible to define ghost atom
            unknown = None
            if norm:
                label = norm.group(1).capitalize()
                return label if label in chemical_symbols[1:] else unknown
            return unknown

        def normalize_number(number: Any, normalized_label: str | None) -> int | None:
            try:
                raw = int(float(number))
            except Exception:
                return None

            # Legacy CRYSTAL parser semantics: NAT atomic numbers are mapped with modulo 100.
            # Example: 238 -> 38 (Sr), and ghost atoms remain 0.
            normalized = raw % 100
            if normalized == 0:
                return 0

            if 0 < normalized < len(chemical_symbols):
                return normalized

            return raw

        atoms = []
        for n, label in enumerate(labels):
            normalized_label = normalize_label(label)
            normalized_number = normalize_number(numbers[n], normalized_label)
            if normalized_label is None and normalized_number is not None:
                if 0 < normalized_number < len(chemical_symbols):
                    normalized_label = chemical_symbols[normalized_number]
            atoms.append(dict(label=normalized_label, number=normalized_number))

        return atoms or None

    def get_xc_functionals(self, source: dict[str, Any]) -> list[dict[str, Any]]:
        xc_functionals = set()
        for key in ['exchange', 'correlation', 'exchange_correlation']:
            xc_functionals.update(self.libxc_map.get(source.get(key), []))
        return [dict(name=xc) for xc in xc_functionals]

    def get_outputs(self, source: dict[str, Any]) -> list[dict[str, Any]]:
        outputs = []
        for n, step in enumerate(source.get('geo_opt', {}).get('geo_opt_step', [])):
            outputs.append(dict(energy=step.get('energy')))

        output0 = outputs[0] if outputs else dict()
        output0.setdefault('energy', source.get('energy_total'))
        forces = source.get('forces')
        if forces is not None:
            output0.setdefault(
                'forces', forces[:, 2:].astype(float) * ureg.hartree / ureg.bohr
            )
        scf_steps = self.get_scf_steps(source)
        if scf_steps:
            output0.setdefault('scf_steps', scf_steps)
        if not outputs and output0:
            outputs.append(output0)

        return outputs

    def get_scf_steps(self, source: dict[str, Any]) -> dict[str, Any]:
        scf_iterations = source.get('scf_block', {}).get('scf_iterations', [])
        if not scf_iterations:
            return {}

        energies_total = []
        delta_energies_total = []
        charge_normalization_factor = []
        for scf_step in scf_iterations:
            energies = scf_step.get('energies')
            if energies is not None and len(energies) > 0:
                energies_total.append(energies[0])
                if len(energies) > 1:
                    delta_energies_total.append(abs(energies[1]))

            charge_norm = scf_step.get('charge_normalization_factor')
            if charge_norm is not None:
                charge_normalization_factor.append(float(charge_norm))

        if not energies_total:
            return {}

        scf_steps = {'energies_total': energies_total}
        if delta_energies_total:
            scf_steps['delta_energies_total'] = delta_energies_total

        code_specific_quantities = {}
        n_scf_steps = source.get('number_of_scf_iterations')
        n_scf_steps_max = source.get('scf_max_iteration')
        if n_scf_steps is not None:
            code_specific_quantities['n_scf_steps'] = int(n_scf_steps)
        if n_scf_steps_max is not None:
            code_specific_quantities['n_scf_steps_max'] = int(n_scf_steps_max)
        if len(charge_normalization_factor) == len(energies_total):
            code_specific_quantities['charge_normalization_factor'] = (
                charge_normalization_factor
            )
        if code_specific_quantities:
            scf_steps['code_specific_quantities'] = code_specific_quantities
        return scf_steps

    def build_workflow(self, source: dict[str, Any]):
        scf_threshold = source.get('scf_threshold_energy_change')
        if source.get('geo_opt') is not None:
            workflow = GeometryOptimization()
            workflow.method = GeometryOptimizationMethod()

            energy_change = source.get('energy_change')
            if energy_change is not None:
                workflow.method.convergence_targets = [
                    EnergyConvergenceTarget(
                        threshold=energy_change,
                        threshold_type='absolute',
                    )
                ]

            if scf_threshold is not None:
                workflow.method.single_point_convergence_targets = [
                    EnergyConvergenceTarget(
                        threshold=scf_threshold,
                        threshold_type='absolute',
                    )
                ]
            return workflow

        workflow = SinglePoint()
        workflow.method = SinglePointMethod()
        if scf_threshold is not None:
            workflow.method.convergence_targets = [
                EnergyConvergenceTarget(
                    threshold=scf_threshold,
                    threshold_type='absolute',
                )
            ]
        return workflow

    def get_systems(self, source: dict[str, Any]) -> list[dict[str, Any]]:
        def get_pbc(system_source: dict[str, Any]) -> list[bool] | None:
            lattice_vectors = self.get_lattice_vectors(system_source)
            if lattice_vectors is None:
                return None
            dimensionality = int(self.data.get('dimensionality', 3) or 3)
            dimensionality = max(0, min(3, dimensionality))
            return [axis < dimensionality for axis in range(3)]

        initial = source.get('system_edited', source)
        systems = [
            # initial system
            dict(
                positions=self.get_positions(initial),
                atoms=self.get_atoms(initial),
                lattice_vectors=self.get_lattice_vectors(initial),
                periodic_boundary_conditions=get_pbc(initial),
            )
        ]
        # skip first step same as initial
        for step in source.get('geo_opt', {}).get('geo_opt_step', [])[1:]:
            systems.append(
                dict(
                    positions=self.get_positions(step),
                    atoms=self.get_atoms(step),
                    lattice_vectors=self.get_lattice_vectors(step),
                    periodic_boundary_conditions=get_pbc(step),
                )
            )
        return systems

class CrystalF25Parser(TextParser):
    @property
    def logger(self):
        pass

    @staticmethod
    def to_array(cols: int, rows: int, values: str) -> np.ndarray:
        """Transforms the Crystal-specific f25 array syntax into a numpy array."""
        values = values.replace('\n', '').replace('\r', '')
        values = [values[n : n + 12] for n in range(0, len(values), 12)]
        # does not seem to work
        # values = textwrap.wrap(values, 12)
        return np.array(values, dtype=np.float64).reshape((rows, cols))

    def get_dos(self, source: dict[str, Any]) -> list[dict[str, Any]]:
        first_row = source['first_row']
        cols, rows = (int(first_row[n]) for n in range(2))
        de = first_row[3]
        second_row = source['second_row']
        start_energy = second_row[1]
        dos_values = self.to_array(cols, rows, source['values']).T
        return [
            dict(
                energies=(start_energy + np.arange(rows) * de) * ureg.hartree,
                values=dos_values[n],
            )
            for n in range(len(dos_values))
        ]

    def get_band_structures(self, source: list[dict[str, Any]]) -> list[dict[str, Any]]:
        if not source:
            return []

        band_structures = []
        for segment in source:
            first_row = segment.get('first_row')
            energies = segment.get('energies')
            if first_row is None or energies is None:
                continue

            cols, rows = (int(first_row[n]) for n in range(2))
            values = self.to_array(cols, rows, energies)
            band_structures.append(dict(value=values[None, :]))

        return band_structures


class CrystalMetainfoParser(MetainfoParser):
    @property
    def logger(self):
        pass


class CrystalArchiveWriter(ArchiveWriter):
    output_parser = CrystalOutputParser(text_parser=OutputParser())
    f25_parser = CrystalF25Parser(text_parser=F25Parser())
    archive_parser = CrystalMetainfoParser()

    def write_to_archive(self):
        # main output file
        self.archive_parser.annotation_key = crystal.OUT_KEY
        self.archive_parser.data_object = Simulation(program=Program(name='Crystal'))

        self.output_parser.filepath = self.mainfile
        self.output_parser.convert(self.archive_parser)

        self.archive.data = self.archive_parser.data_object
        outputs_before_f25 = list(self.archive.data.outputs or [])
        self.archive.workflow2 = self.output_parser.build_workflow(
            self.output_parser.data
        )

        f25_filepath = self.output_parser.data.get(
            'f25_filepath1', self.output_parser.data.get('f25_filepath2')
        )
        if f25_filepath:
            # parser f25 file
            self.archive_parser.annotation_key = crystal.F25_KEY

            self.f25_parser.filepath = os.path.join(
                os.path.dirname(self.mainfile), os.path.basename(f25_filepath)
            )

            self.f25_parser.convert(self.archive_parser)
            outputs_after_f25 = self.archive.data.outputs or []
            if outputs_before_f25:
                if not outputs_after_f25:
                    self.archive.data.outputs = outputs_before_f25
                else:
                    for idx, output in enumerate(outputs_before_f25):
                        if idx >= len(outputs_after_f25):
                            outputs_after_f25.append(output)
                            continue
                        if (
                            outputs_after_f25[idx].scf_steps is None
                            and output.scf_steps is not None
                        ):
                            outputs_after_f25[idx].scf_steps = output.scf_steps


class CrystalParser(MatchingParser):
    archive_writer = CrystalArchiveWriter()

    def parse(
        self,
        mainfile: str,
        archive: EntryArchive,
        logger: BoundLogger,
        child_arhives: dict[str, EntryArchive] = {},
    ):
        self.archive_writer.write(mainfile, archive, logger, child_arhives)
