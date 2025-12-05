#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD.
# See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

from collections.abc import Iterable
from typing import Any

import numpy as np
from nomad.parsing.file_parser import ArchiveWriter
from nomad.utils import get_logger

# from runschema.method import Interaction, Model
from nomad_simulations.schema_packages import workflow
from nomad_simulations.schema_packages.atoms_state import ParticleState
from nomad_simulations.schema_packages.general import Simulation
from nomad_simulations.schema_packages.model_system import (
    AlternativeRepresentation,
    ModelSystem,
)

# nomad-simulations
from nomad_simulations.schema_packages.outputs import (
    TotalEnergy,
    TotalForce,
    TrajectoryOutputs,
)
from nomad_simulations.schema_packages.properties.energies import BaseEnergy
from nomad_simulations.schema_packages.properties.forces import BaseForce

# from runschema.method import Interaction, Model


class MDParser(ArchiveWriter):
    def __init__(self, **kwargs) -> None:
        self.info: dict[str, Any] = {}
        self.cum_max_atoms: int = 2500000
        self.logger = get_logger(__name__)
        self._trajectory_steps: list[int] = []
        self._thermodynamics_steps: list[int] = []
        self._trajectory_steps_sampled: list[int] = []
        self._steps: list[int] = []
        super().__init__(**kwargs)
        self._logger = get_logger(__name__)

    @property
    def logger(self):
        return self._logger

    @logger.setter
    def logger(self, value):
        self._logger = value

    @property
    def steps(self) -> list[int]:
        """
        Returns the set of trajectory and thermodynamics steps.
        """
        if not self._steps:
            self._steps = list(set(self.trajectory_steps + self.thermodynamics_steps))
            self._steps.sort()
        return self._steps

    @property
    def trajectory_steps(self) -> list[int]:
        """
        Returns the sampled trajectory steps.
        """
        if not self._trajectory_steps_sampled:
            self._trajectory_steps_sampled = [
                step
                for n, step in enumerate(self._trajectory_steps)
                if n % self.archive_sampling_rate == 0
            ]
        return self._trajectory_steps_sampled

    @trajectory_steps.setter
    def trajectory_steps(self, value: list[int]):
        self._trajectory_steps = list(set(value))
        self._trajectory_steps.sort()
        self.info['n_frames'] = len(self._trajectory_steps)
        self._trajectory_steps_sampled = []

    @property
    def thermodynamics_steps(self) -> list[int]:
        """
        Returns the thermodynamics steps.
        """
        # TODO is it necessary to sample thermodynamics steps
        return self._thermodynamics_steps

    @thermodynamics_steps.setter
    def thermodynamics_steps(self, value: list[int]):
        self._thermodynamics_steps = list(set(value))
        self._thermodynamics_steps.sort()

    @property
    def n_atoms(self) -> int:
        return np.amax(self.info.get('n_atoms', [0]))

    @n_atoms.setter
    def n_atoms(self, value: Iterable | int):
        self.info['n_atoms'] = [value] if not isinstance(value, Iterable) else value

    @property
    def archive_sampling_rate(self) -> int:
        """
        Returns the sampling rate of saved thermodynamics data and trajectory.
        """
        if self.info.get('archive_sampling_rate') is None:
            n_frames = self.info.get('n_frames', len(self._trajectory_steps))
            n_atoms = np.amax(self.n_atoms)
            if not n_atoms or not n_frames:
                self.info['archive_sampling_rate'] = 1
            else:
                cum_atoms = n_atoms * n_frames
                self.info['archive_sampling_rate'] = (
                    1
                    if cum_atoms <= self.cum_max_atoms
                    else -(-cum_atoms // self.cum_max_atoms)
                )
        return self.info.get('archive_sampling_rate')

    def parse(self, *args, **kwargs):
        self.info = {}
        self.trajectory_steps = []
        self.thermodynamics_steps = []
        self._steps = []
        self._trajectory_steps_sampled = []
        super().parse(*args, **kwargs)

    def parse_trajectory_step(
        self,
        data: dict[str, Any],
        simulation: Simulation,
        model_system: ModelSystem = None,
        representation: AlternativeRepresentation = None,
    ) -> None:
        """
        Create a system section and write the provided data.
        """
        if simulation is None:
            return

        if (step := data.get('step')) is not None and step not in self.trajectory_steps:
            return
        if model_system is None:
            model_system = ModelSystem()
        if representation is None:
            representation = AlternativeRepresentation()

        lattice_vectors = data.pop('lattice_vectors')
        periodic_boundary_conditions = data.pop('periodic_boundary_conditions')
        particle_labels = data.pop('labels')
        dimensions = data.pop('dimensions')

        for label in particle_labels:
            particle_state = ParticleState(label=label)
            model_system.particle_states.append(particle_state)

        q_lattice = model_system.m_get_quantity_definition('lattice_vectors')
        model_system.m_set(q_lattice, lattice_vectors)

        q_pbc = model_system.m_get_quantity_definition('periodic_boundary_conditions')
        model_system.m_set(q_pbc, periodic_boundary_conditions)

        model_system.dimensionality = dimensions

        self.parse_section(data, model_system)
        simulation.model_system.append(model_system)

    def parse_output_step(
        self,
        data: dict[str, Any],
        simulation: Simulation,
        output: TrajectoryOutputs = None,
    ) -> bool:
        if self.archive is None:
            return False

        if (
            step := data.get('step')
        ) is not None and step not in self.thermodynamics_steps:
            return False

        if output is None:
            output = TrajectoryOutputs()

        energy_contributions = data.get('total_energies', {}).pop('contributions', {})
        force_contributions = data.get('total_forces', {}).pop('contributions', {})
        self.parse_section(data, output)
        try:
            system_ref_index = self.trajectory_steps.index(output.step)
            output.model_system_ref = simulation.model_system[system_ref_index]
        except Exception:
            self.logger.warning('Could not set system reference in parsing of outputs.')

        if energy_contributions:
            if len(output.total_energies) == 0:
                output.total_energies.append(TotalEnergy())

        for energy_dict in energy_contributions:
            energy = BaseEnergy()
            output.total_energies[-1].contributions.append(energy)
            self.parse_section(energy_dict, energy)

        if force_contributions:
            if len(output.total_forces) == 0:
                output.total_forces.append(TotalForce())

        for force_dict in force_contributions:
            force = BaseForce
            output.total_forces[-1].contributions.append(force)
            self.parse_section(force_dict, force)

        simulation.outputs.append(output)

        return True

    def parse_md_workflow(self, data: dict[str, Any]) -> None:
        """
        Create an md workflow section and write the provided data.
        """
        if self.archive is None:
            return

        sec_workflow = workflow.MolecularDynamics()
        self.parse_section(data, sec_workflow)
        self.archive.workflow2 = sec_workflow

    # TODO Adapt these interaction functions for the new schema
    # def parse_interactions(self, interactions: List[Dict], sec_model: MSection)
    # -> None:
    #     if not interactions:
    #         return

    #     def write_interaction_values(values):
    #         sec_interaction = Interaction()
    #         sec_model.contributions.append(sec_interaction)
    #         sec_interaction.type = current_type
    #         sec_interaction.n_atoms = max(
    #             [len(v) for v in values.get('atom_indices', [[0]])]
    #         )
    #         for key, val in values.items():
    #             quantity_def = sec_interaction.m_def.all_quantities.get(key)
    #             if quantity_def:
    #                 try:
    #                     sec_interaction.m_set(quantity_def, val)
    #                 except Exception:
    #                     self.logger.error('Error setting metadata.',
    #  data={'key': key})

    #     interactions.sort(key=lambda x: x.get('type'))
    #     current_type = interactions[0].get('type')
    #     interaction_values: Dict[str, Any] = {}
    #     for interaction in interactions:
    #         interaction_type = interaction.get('type')
    #         if current_type and current_type != interaction_type:
    #             write_interaction_values(interaction_values)
    #             current_type = interaction_type
    #             interaction_values = {}
    #         interaction_values.setdefault('n_interactions', 0)
    #         interaction_values['n_interactions'] += 1
    #         for key, val in interaction.items():
    #             if key == 'type':
    #                 continue
    #             interaction_values.setdefault(key, [])
    #             interaction_values[key].append(val)
    #     if interaction_values:
    #         write_interaction_values(interaction_values)

    # def parse_interactions_by_type(
    #     self, interactions_by_type: List[Dict], sec_model: Model
    # ) -> None:
    #     for interaction_type_dict in interactions_by_type:
    #         sec_interaction = Interaction()
    #         sec_model.contributions.append(sec_interaction)
    #         self.parse_section(interaction_type_dict, sec_interaction)
    #  # TODO Shift Gromacs and Lammps parsers to use this function as well if possible
