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
import numpy as np
from nomad.datamodel.data import ArchiveSection
from nomad.datamodel.metainfo.annotations import Mapper as MapperAnnotation
from nomad.metainfo import Quantity, SchemaPackage, Section, SubSection
from nomad_simulations.schema_packages import (
    atoms_state,
    general,
    model_system,
    outputs,
    physical_property,
    properties,
)

# from simulationworkflowschema import molecular_dynamics
from nomad_simulations.schema_packages.workflow import molecular_dynamics, trajectory

from nomad_simulation_parsers.schema_packages.utils import add_mapping_annotation

m_package = SchemaPackage()

HDF5_KEY = 'hdf5'

# Global list of standard H5MD observables to exclude from custom outputs
STANDARD_H5MD_OBSERVABLES = [
    'energies',
    'temperatures',
    'forces',
    'custom_forces',
    'radial_distribution_functions',
    'mean_squared_displacements',
    'diffusion_constants',
]


# SIMULATION --> archive.data


class ParamEntry(ArchiveSection):
    """
    Generic section defining a parameter name and value
    """

    name = Quantity(
        type=str,
        shape=[],
        description="""
        Name of the parameter.
        """,
    )

    value = Quantity(
        type=str,
        shape=[],
        description="""
        Value of the parameter as a string.
        """,
    )

    unit = Quantity(
        type=str,
        shape=[],
        description="""
        Unit of the parameter as a string.
        """,
    )

    description = Quantity(
        type=str,
        shape=[],
        description="""
        Further description of the attribute.
        """,
    )


class Author(ArchiveSection):
    """
    Contains the specifications of the program.
    """

    name = Quantity(
        type=str,
        shape=[],
        description="""
        Specifies the name of the author who generated the h5md file.
        """,
    )

    add_mapping_annotation(name, HDF5_KEY, r'."@name"')

    email = Quantity(
        type=str,
        shape=[],
        description="""
        Author's email.
        """,
    )

    add_mapping_annotation(email, HDF5_KEY, r'."@email"')


## class Program(general.Program):

add_mapping_annotation(general.Program.name, HDF5_KEY, r'."@name"')
add_mapping_annotation(general.Program.version, HDF5_KEY, r'."@version"')

# SIMULATION.MODEL_SYSTEM --> archive.data.model_system

# ! Only use generic ParticleState
# ! ModelSystem normalizer takes care of assigning to AtomState or CGBeadState
### class ParticleState(atoms_state.ParticleState):
add_mapping_annotation(
    atoms_state.ParticleState.m_def,
    HDF5_KEY,
    ('to_species_labels', ['.@'], dict(path='particles.all.species_label')),
)
add_mapping_annotation(atoms_state.ParticleState.label, HDF5_KEY, '.label')

### Geometric properties now on ModelSystem directly (handled in get_traj_data)

add_mapping_annotation(
    model_system.ModelSystem.lattice_vectors, HDF5_KEY, '.lattice_vectors'
)

add_mapping_annotation(
    model_system.ModelSystem.periodic_boundary_conditions, HDF5_KEY, '.boundary'
)


class ModelSystem(model_system.ModelSystem):
    """
    Model system used as an input for simulating the material.
    """

    m_def = Section(
        validate=False,
    )

    custom_system_attributes = (
        SubSection(  # TODO should this be called parameters or attributes or what?
            sub_section=ParamEntry.m_def,
            description="""
        Contains additional information about the (sub)system .
        """,
            repeats=True,
        )
    )


# TODO add annotaiton for custom_system_attributes
# add_mapping_annotations(ModelSystem.custom_system_attributes, HDF5_KEY, (
#         'get_custom_system_attributes',
#         ['.@'],
#         dict(path=''),
#     )
# )

add_mapping_annotation(ModelSystem.n_particles, HDF5_KEY, '.n_particles')

add_mapping_annotation(
    ModelSystem.bond_list,
    HDF5_KEY,
    ('get_top_system_quantity', ['.@'], dict(path='connectivity.bonds')),
)

add_mapping_annotation(
    ModelSystem.dimensionality,
    HDF5_KEY,
    ('get_top_system_quantity', ['.@'], dict(path='particles.all.box.@dimension')),
)

add_mapping_annotation(ModelSystem.positions, HDF5_KEY, '.positions')

add_mapping_annotation(ModelSystem.velocities, HDF5_KEY, '.velocities')

add_mapping_annotation(
    ModelSystem.sub_systems,
    HDF5_KEY,
    ('get_sub_systems', ['.@'], dict(path='connectivity')),
)

add_mapping_annotation(ModelSystem.name, HDF5_KEY, '.label')

add_mapping_annotation(ModelSystem.composition_formula, HDF5_KEY, '.formula')

add_mapping_annotation(ModelSystem.particle_indices, HDF5_KEY, '.indices')

add_mapping_annotation(ModelSystem.branch_label, HDF5_KEY, '.type')

# TODO add function to filter out valid Enums
# add_mapping_annotations(ModelSystem.type, HDF5_KEY, '.type')


### SUBSECTIONS


# TODO need to add ParticleCell and distinguish in the parser
#### model_system.cell --> AtomicCell


# SIMULATION.METHOD --> archive.data.method

# TODO Add method, including full FF example
# class ForceCalculations(runschema.method.ForceCalculations):
#     m_def = Section(
#         validate=False,
#         extends_base_section=True,
#     )

#     x_h5md_parameters = SubSection(
#         sub_section=ParamEntry.m_def,
#         description="""
#         Contains non-normalized force calculation parameters.
#         """,
#         repeats=True,
#     )


# class NeighborSearching(runschema.method.NeighborSearching):
#     m_def = Section(
#         validate=False,
#         extends_base_section=True,
#     )

#     x_h5md_parameters = SubSection(
#         sub_section=ParamEntry.m_def,
#         description="""
#         Contains non-normalized neighbor searching parameters.
#         """,
#         repeats=True,
#     )


# SIMULATION.OUTPUTS --> archive.data.outputs


class CustomProperty(physical_property.PhysicalProperty):
    """
    Section describing a general type of ouput property.
    """

    name = Quantity(
        type=str,
        shape=[],
        description="""
        Name of the parameter.
        """,
    )
    add_mapping_annotation(name, HDF5_KEY, '.name')

    value = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description="""
        Value **magnitude** of the property.
        The unit is defined in the `unit` attribute.
        """,
    )
    add_mapping_annotation(value, HDF5_KEY, '.value')

    unit = Quantity(
        type=str,
        shape=[],
        description="""
        Unit of the parameter as a string consistent with the UnitRegistry.pint module.
        """,
    )
    add_mapping_annotation(unit, HDF5_KEY, '.unit')

    description = Quantity(
        type=str,
        shape=[],
        description="""
        Further description of the property.
        """,
    )
    add_mapping_annotation(description, HDF5_KEY, '.description')

    def normalize(self, *args, **kwargs) -> None:
        """
        Let parent do all the work, but prevent it from overwriting parsed names.
        This is more surgical than overriding the whole method.
        """
        # Temporarily hide the class name to prevent overwrite
        original_name = self.m_def.name
        self.m_def.name = None

        # Call parent normalization
        super().normalize(*args, **kwargs)

        # Restore the class name for future use
        self.m_def.name = original_name


# Add mapping annotation for the inherited name field
add_mapping_annotation(CustomProperty.name, HDF5_KEY, '.name')


## class BaseEnergy(properties.energies.BaseEnergy):

# ? Does this note still apply?
# value annotation defined in TotalEnergy.value since they refer to the same quantity
# in this case, we make sure to return the corresponding value from
# the get_contributions function in the TotalEnergy.contributions annotation
add_mapping_annotation(
    properties.energies.BaseEnergy.contribution_type, HDF5_KEY, '.name'
)


## class TotalEnergy(properties.TotalEnergy):


add_mapping_annotation(
    properties.TotalEnergy.value,
    HDF5_KEY,
    (
        'get_output_data',
        ['.@'],
        dict(path='observables.energies.total', observable_type='configurational'),
    ),
)

### SUBSECTIONS

add_mapping_annotation(
    properties.TotalEnergy.contributions,
    HDF5_KEY,
    (
        'get_contributions',
        ['.@'],
        dict(path='observables.energies', exclude=['total']),
    ),
)


## class BaseForce(properties.forces.BaseForce):


add_mapping_annotation(properties.forces.BaseForce.contribution_type, HDF5_KEY, '.name')


## class TotalForce(properties.TotalForce):


add_mapping_annotation(
    properties.TotalForce.value,
    HDF5_KEY,
    (
        'get_output_data',
        ['.@'],
        dict(path='observables.forces.total', observable_type='configurational'),
    ),
)

### SUBSECTIONS

add_mapping_annotation(
    properties.TotalForce.contributions,
    HDF5_KEY,
    (
        'get_contributions',
        ['.@'],
        dict(path='observables.forces', exclude=['total']),
    ),
)


## class Temperature(properties.Temperature):


add_mapping_annotation(
    properties.Temperature.value,
    HDF5_KEY,
    (
        'get_output_data',
        ['.@'],
        dict(path='observables.temperatures', observable_type='configurational'),
    ),
)


class TrajectoryOutputs(outputs.TrajectoryOutputs):
    m_def = Section(
        validate=False,
    )

    custom_outputs = SubSection(
        sub_section=CustomProperty.m_def,
        description="""
        Contains other generic custom outputs that are not already defined.
        """,
        repeats=True,
    )

    total_forces = SubSection(sub_section=properties.TotalForce.m_def, repeats=True)


add_mapping_annotation(
    TrajectoryOutputs.m_def, HDF5_KEY, ('get_output_steps', ['observables'])
)

add_mapping_annotation(TrajectoryOutputs.step, HDF5_KEY, '.step')

add_mapping_annotation(TrajectoryOutputs.time, HDF5_KEY, '.time')

### SUBSECTIONS

add_mapping_annotation(TrajectoryOutputs.total_energies, HDF5_KEY, '.@')

add_mapping_annotation(TrajectoryOutputs.total_forces, HDF5_KEY, '.@')

add_mapping_annotation(TrajectoryOutputs.temperatures, HDF5_KEY, '.@')

add_mapping_annotation(
    TrajectoryOutputs.custom_outputs,
    HDF5_KEY,
    (
        'get_custom_outputs',
        ['.@'],
        dict(
            path='observables',
            exclude=[
                'energies',
                'temperatures',
                'custom_forces',
            ],  # TODO get the exclusion list automatically
            observable_type='configurational',
        ),
    ),
)


class Simulation(general.Simulation):
    m_def = Section(
        validate=False,
    )

    # TODO Not sure how we are dealing with versioning with H5MD-NOMAD
    x_h5md_version = Quantity(
        type=np.dtype(np.int32),
        shape=[2],
        description="""
        Specifies the version of the h5md schema being followed.
        """,
    )
    add_mapping_annotation(x_h5md_version, HDF5_KEY, r'h5md."@version"')

    x_h5md_author = SubSection(sub_section=Author.m_def)

    add_mapping_annotation(x_h5md_author, HDF5_KEY, 'h5md.author')

    x_h5md_creator = SubSection(sub_section=general.Program.m_def)

    add_mapping_annotation(x_h5md_creator, HDF5_KEY, 'h5md.creator')

    model_system = SubSection(sub_section=ModelSystem.m_def, repeats=True)


add_mapping_annotation(Simulation.m_def, HDF5_KEY, '@')

add_mapping_annotation(Simulation.program, HDF5_KEY, 'h5md.program')

add_mapping_annotation(
    Simulation.model_system, HDF5_KEY, ('get_traj_data', ['particles.all'])
)


# WORKFLOW --> archive.workflow2
h5md_path_md = 'parameters.workflow.molecular_dynamics'

# WORKFLOW.METHOD --> archive.workflow2.method


## class ThermostatParameters(molecular_dynamics.ThermostatParameters):


h5md_path_thermostat = f'{h5md_path_md}.thermostat_parameters'

add_mapping_annotation(
    molecular_dynamics.ThermostatParameters.thermostat_type,
    HDF5_KEY,
    (
        'map_value',
        [h5md_path_thermostat],
        dict(key='thermostat_type', enum_spec='lower'),
    ),
)

add_mapping_annotation(
    molecular_dynamics.ThermostatParameters.reference_temperature,
    HDF5_KEY,
    (
        'map_value',
        [h5md_path_thermostat],
        dict(key='reference_temperature'),
    ),
)

add_mapping_annotation(
    molecular_dynamics.ThermostatParameters.coupling_constant,
    HDF5_KEY,
    (
        'map_value',
        [h5md_path_thermostat],
        dict(key='coupling_constant'),
    ),
)

add_mapping_annotation(
    molecular_dynamics.ThermostatParameters.effective_mass,
    HDF5_KEY,
    (
        'map_value',
        [h5md_path_thermostat],
        dict(key='effective_mass'),
    ),
)

add_mapping_annotation(
    molecular_dynamics.ThermostatParameters.temperature_profile,
    HDF5_KEY,
    (
        'map_value',
        [h5md_path_thermostat],
        dict(key='temperature_profile', enum_spec='lower'),
    ),
)

add_mapping_annotation(
    molecular_dynamics.ThermostatParameters.reference_temperature_start,
    HDF5_KEY,
    (
        'map_value',
        [h5md_path_thermostat],
        dict(key='reference_temperature_start'),
    ),
)

add_mapping_annotation(
    molecular_dynamics.ThermostatParameters.reference_temperature_end,
    HDF5_KEY,
    (
        'map_value',
        [h5md_path_thermostat],
        dict(key='reference_temperature_end'),
    ),
)

add_mapping_annotation(
    molecular_dynamics.ThermostatParameters.temperature_update_frequency,
    HDF5_KEY,
    (
        'map_value',
        [h5md_path_thermostat],
        dict(key='temperature_update_frequency'),
    ),
)

add_mapping_annotation(
    molecular_dynamics.ThermostatParameters.temperature_update_delta,
    HDF5_KEY,
    (
        'map_value',
        [h5md_path_thermostat],
        dict(key='temperature_update_delta'),
    ),
)

add_mapping_annotation(
    molecular_dynamics.ThermostatParameters.temperature_update_factor,
    HDF5_KEY,
    (
        'map_value',
        [h5md_path_thermostat],
        dict(key='temperature_update_factor'),
    ),
)

add_mapping_annotation(
    molecular_dynamics.ThermostatParameters.step_start,
    HDF5_KEY,
    (
        'map_value',
        [h5md_path_thermostat],
        dict(key='step_start'),
    ),
)

add_mapping_annotation(
    molecular_dynamics.ThermostatParameters.step_end,
    HDF5_KEY,
    (
        'map_value',
        [h5md_path_thermostat],
        dict(key='step_end'),
    ),
)


## class BarostatParameters(molecular_dynamics.BarostatParameters):


h5md_path_barostat = f'{h5md_path_md}.barostat_parameters'

add_mapping_annotation(
    molecular_dynamics.BarostatParameters.barostat_type,
    HDF5_KEY,
    (
        'map_value',
        [h5md_path_barostat],
        dict(key='barostat_type', enum_spec='lower'),
    ),
)

add_mapping_annotation(
    molecular_dynamics.BarostatParameters.coupling_type,
    HDF5_KEY,
    (
        'map_value',
        [h5md_path_barostat],
        dict(key='coupling_type', enum_spec='lower'),
    ),
)

add_mapping_annotation(
    molecular_dynamics.BarostatParameters.reference_pressure,
    HDF5_KEY,
    (
        'map_value',
        [h5md_path_barostat],
        dict(key='reference_pressure'),
    ),
)

add_mapping_annotation(
    molecular_dynamics.BarostatParameters.coupling_constant,
    HDF5_KEY,
    (
        'map_value',
        [h5md_path_barostat],
        dict(key='coupling_constant'),
    ),
)

add_mapping_annotation(
    molecular_dynamics.BarostatParameters.compressibility,
    HDF5_KEY,
    (
        'map_value',
        [h5md_path_barostat],
        dict(key='compressibility'),
    ),
)

add_mapping_annotation(
    molecular_dynamics.BarostatParameters.pressure_profile,
    HDF5_KEY,
    (
        'map_value',
        [h5md_path_barostat],
        dict(key='pressure_profile', enum_spec='lower'),
    ),
)

add_mapping_annotation(
    molecular_dynamics.BarostatParameters.reference_pressure_start,
    HDF5_KEY,
    (
        'map_value',
        [h5md_path_barostat],
        dict(key='reference_pressure_start'),
    ),
)

add_mapping_annotation(
    molecular_dynamics.BarostatParameters.reference_pressure_end,
    HDF5_KEY,
    (
        'map_value',
        [h5md_path_barostat],
        dict(key='reference_pressure_end'),
    ),
)

add_mapping_annotation(
    molecular_dynamics.BarostatParameters.pressure_update_frequency,
    HDF5_KEY,
    (
        'map_value',
        [h5md_path_barostat],
        dict(key='pressure_update_frequency'),
    ),
)

add_mapping_annotation(
    molecular_dynamics.BarostatParameters.pressure_update_delta,
    HDF5_KEY,
    (
        'map_value',
        [h5md_path_barostat],
        dict(key='pressure_update_delta'),
    ),
)

add_mapping_annotation(
    molecular_dynamics.BarostatParameters.pressure_update_factor,
    HDF5_KEY,
    (
        'map_value',
        [h5md_path_barostat],
        dict(key='pressure_update_factor'),
    ),
)

add_mapping_annotation(
    molecular_dynamics.BarostatParameters.step_start,
    HDF5_KEY,
    (
        'map_value',
        [h5md_path_barostat],
        dict(key='step_start'),
    ),
)

add_mapping_annotation(
    molecular_dynamics.BarostatParameters.step_end,
    HDF5_KEY,
    (
        'map_value',
        [h5md_path_barostat],
        dict(key='step_end'),
    ),
)


## class ShearParameters(molecular_dynamics.ShearParameters):


h5md_path_shear = f'{h5md_path_md}.shear_parameters'

add_mapping_annotation(
    molecular_dynamics.ShearParameters.shear_type,
    HDF5_KEY,
    (
        'map_value',
        [h5md_path_shear],
        dict(key='shear_type', enum_spec='lower'),
    ),
)

add_mapping_annotation(
    molecular_dynamics.ShearParameters.shear_rate,
    HDF5_KEY,
    (
        'map_value',
        [h5md_path_shear],
        dict(key='shear_rate'),
    ),
)

add_mapping_annotation(
    molecular_dynamics.ShearParameters.step_start,
    HDF5_KEY,
    (
        'map_value',
        [h5md_path_shear],
        dict(key='step_start'),
    ),
)

add_mapping_annotation(
    molecular_dynamics.ShearParameters.step_end,
    HDF5_KEY,
    (
        'map_value',
        [h5md_path_shear],
        dict(key='step_end'),
    ),
)


## class FreeEnergyCalculationParameters(
##     molecular_dynamics.FreeEnergyCalculationParameters
## ):


h5md_path_FEC = f'{h5md_path_md}.free_energy_calculation_parameters'

# TODO Change this to fec_type in the schema
add_mapping_annotation(
    molecular_dynamics.FreeEnergyCalculationParameters.calc_type,
    HDF5_KEY,
    (
        'map_value',
        [h5md_path_FEC],
        dict(key='type', enum_spec='lower'),
    ),
)

add_mapping_annotation(
    molecular_dynamics.FreeEnergyCalculationParameters.current_lambda_index,
    HDF5_KEY,
    (
        'map_value',
        [h5md_path_FEC],
        dict(key='lambda_index'),
    ),
)

# TODO: These properties are not yet implemented in FreeEnergyCalculationParameters
# add_mapping_annotation(
#     molecular_dynamics.FreeEnergyCalculationParameters.atom_indices,
#     HDF5_KEY,
#     (
#         'map_value',
#         [h5md_path_FEC],
#         dict(key='atom_indices'),
#     ),
# )

# add_mapping_annotation(
#     molecular_dynamics.FreeEnergyCalculationParameters.initial_state_vdw,
#     HDF5_KEY,
#     (
#         'map_value',
#         [h5md_path_FEC],
#         dict(key='initial_state_vdw'),
#     ),
# )

# add_mapping_annotation(
#     molecular_dynamics.FreeEnergyCalculationParameters.final_state_vdw,
#     HDF5_KEY,
#     (
#         'map_value',
#         [h5md_path_FEC],
#         dict(key='final_state_vdw'),
#     ),
# )

# add_mapping_annotation(
#     molecular_dynamics.FreeEnergyCalculationParameters.initial_state_coloumb,
#     HDF5_KEY,
#     (
#         'map_value',
#         [h5md_path_FEC],
#         dict(key='initial_state_coloumb'),
#     ),
# )

# add_mapping_annotation(
#     molecular_dynamics.FreeEnergyCalculationParameters.final_state_coloumb,
#     HDF5_KEY,
#     (
#         'map_value',
#         [h5md_path_FEC],
#         dict(key='final_state_coloumb'),
#     ),
# )

# add_mapping_annotation(
#     molecular_dynamics.FreeEnergyCalculationParameters.initial_state_bonded,
#     HDF5_KEY,
#     (
#         'map_value',
#         [h5md_path_FEC],
#         dict(key='initial_state_bonded'),
#     ),
# )

# add_mapping_annotation(
#     molecular_dynamics.FreeEnergyCalculationParameters.final_state_bonded,
#     HDF5_KEY,
#     (
#         'map_value',
#         [h5md_path_FEC],
#         dict(key='final_state_bonded'),
#     ),
# )

### SUBSECTIONS

add_mapping_annotation(
    molecular_dynamics.FreeEnergyCalculationParameters.lambdas, HDF5_KEY, '@'
)


## class Lambdas(molecular_dynamics.Lambdas):
# TODO add lambdas to test data

# ? Not sure about where this info is going in h5md
h5md_path_lambdas = f'{h5md_path_FEC}.lambdas'

add_mapping_annotation(
    molecular_dynamics.Lambdas.interaction_type,
    HDF5_KEY,
    (
        'map_value',
        [h5md_path_lambdas],
        dict(key='type', enum_spec='lower'),
    ),
)

add_mapping_annotation(
    molecular_dynamics.Lambdas.values,
    HDF5_KEY,
    (
        'map_value',
        [h5md_path_lambdas],
        dict(key='value'),
    ),
)


## class MolecularDynamicsMethod(molecular_dynamics.MolecularDynamicsMethod):


add_mapping_annotation(
    molecular_dynamics.MolecularDynamicsMethod.thermodynamic_ensemble,
    HDF5_KEY,
    (
        'map_value',
        [h5md_path_md],
        dict(key='thermodynamic_ensemble', enum_spec='upper'),
    ),
)

add_mapping_annotation(
    molecular_dynamics.MolecularDynamicsMethod.integrator_type,
    HDF5_KEY,
    (
        'map_value',
        [h5md_path_md],
        dict(key='integrator_type', enum_spec='lower'),
    ),
)

add_mapping_annotation(
    molecular_dynamics.MolecularDynamicsMethod.integration_timestep,
    HDF5_KEY,
    (
        'map_value',
        [h5md_path_md],
        dict(key='integration_timestep'),
    ),
)

add_mapping_annotation(
    molecular_dynamics.MolecularDynamicsMethod.n_steps,
    HDF5_KEY,
    (
        'map_value',
        [h5md_path_md],
        dict(key='n_steps'),
    ),
)

add_mapping_annotation(
    molecular_dynamics.MolecularDynamicsMethod.coordinate_save_frequency,
    HDF5_KEY,
    (
        'map_value',
        [h5md_path_md],
        dict(key='coordinate_save_frequency'),
    ),
)

add_mapping_annotation(
    molecular_dynamics.MolecularDynamicsMethod.velocity_save_frequency,
    HDF5_KEY,
    (
        'map_value',
        [h5md_path_md],
        dict(key='velocity_save_frequency'),
    ),
)

add_mapping_annotation(
    molecular_dynamics.MolecularDynamicsMethod.force_save_frequency,
    HDF5_KEY,
    (
        'map_value',
        [h5md_path_md],
        dict(key='force_save_frequency'),
    ),
)

add_mapping_annotation(
    molecular_dynamics.MolecularDynamicsMethod.thermodynamics_save_frequency,
    HDF5_KEY,
    (
        'map_value',
        [h5md_path_md],
        dict(key='thermodynamics_save_frequency'),
    ),
)

### SUBSECTIONS

add_mapping_annotation(
    molecular_dynamics.MolecularDynamicsMethod.thermostat_parameters, HDF5_KEY, '@'
)

add_mapping_annotation(
    molecular_dynamics.MolecularDynamicsMethod.barostat_parameters, HDF5_KEY, '@'
)

add_mapping_annotation(
    molecular_dynamics.MolecularDynamicsMethod.shear_parameters, HDF5_KEY, '@'
)

add_mapping_annotation(
    molecular_dynamics.MolecularDynamicsMethod.free_energy_calculation_parameters,
    HDF5_KEY,
    '@',
)

# WORKFLOW.RESULTS --> archive.workflow2.results


class ConfigurationalProperty(trajectory.ConfigurationalProperty):
    value_magnitude = Quantity(
        type=np.float64,
        shape=['n_frames'],
        description="""
        Values of the property.
        """,
    )

    value_unit = Quantity(
        type=str,
        shape=[],
        description="""
        Unit of the property, using UnitRegistry() notation.
        """,
    )


class EnsembleProperty(molecular_dynamics.EnsembleProperty):
    bins_magnitude = Quantity(
        type=np.float64,
        shape=['n_bins'],
        description="""
        Values of the variable along which the property is calculated.
        """,
    )

    bins_unit = Quantity(
        type=str,
        shape=[],
        description="""
        Unit of the given bins, using UnitRegistry() notation.
        """,
    )

    value_magnitude = Quantity(
        type=np.float64,
        shape=['n_bins'],
        description="""
        Values of the property.
        """,
    )

    value_unit = Quantity(
        type=str,
        shape=[],
        description="""
        Unit of the property, using UnitRegistry() notation.
        """,
    )


class CorrelationFunction(molecular_dynamics.CorrelationFunction):
    value_magnitude = Quantity(
        type=np.float64,
        shape=['n_times'],
        description="""
        Values of the property.
        """,
    )

    value_unit = Quantity(
        type=str,
        shape=[],
        description="""
        Unit of the property, using UnitRegistry() notation.
        """,
    )


class MolecularDynamicsResults(molecular_dynamics.MolecularDynamicsResults):
    ensemble_properties = SubSection(sub_section=EnsembleProperty.m_def, repeats=True)

    correlation_functions = SubSection(
        sub_section=CorrelationFunction.m_def, repeats=True
    )


# ? These quantities from normalization?
#     finished_normally = Quantity(
#         type=bool,
#         shape=[],
#         description="""
#         Indicates if calculation terminated normally.
#         """,
#     )

#     n_steps = Quantity(
#         type=np.int32,
#         shape=[],
#         description="""
#         Number of trajectory steps""",
#     )

#     trajectory = Quantity(
#         type=Reference(System),
#         shape=['n_steps'],
#         description="""
#         Reference to the system of each step in the trajectory.
#         """,
#     )


class MolecularDynamics(molecular_dynamics.MolecularDynamics):
    results = SubSection(sub_section=MolecularDynamicsResults.m_def)


### SUBSECTIONS

# Custom Ensemble Properties mapping using get_custom_ensemble_outputs
MolecularDynamicsResults.ensemble_properties.m_annotations.setdefault('mapping', {})[
    'hdf5'
] = MapperAnnotation(
    mapper=(
        'get_custom_ensemble_outputs',
        ['observables'],
        dict(
            exclude=STANDARD_H5MD_OBSERVABLES,
            type_filter=['ensemble_average'],
        ),
    )
)

# Individual field mappings for custom EnsembleProperty
add_mapping_annotation(EnsembleProperty.label, HDF5_KEY, '.label')

add_mapping_annotation(EnsembleProperty.bins_magnitude, HDF5_KEY, '.bins')

add_mapping_annotation(EnsembleProperty.bins_unit, HDF5_KEY, '.bins_unit')

add_mapping_annotation(EnsembleProperty.value_magnitude, HDF5_KEY, '.value_magnitude')

add_mapping_annotation(EnsembleProperty.value_unit, HDF5_KEY, '.value_unit')

# Radial Distribution Functions mapping using get_output_data
MolecularDynamicsResults.radial_distribution_functions.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(
    mapper=(
        'get_output_data',
        ['observables.radial_distribution_functions'],
    )
)

# Individual field mappings for RadialDistributionFunction
molecular_dynamics.RadialDistributionFunction.label.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(mapper='.label')

molecular_dynamics.RadialDistributionFunction.bins.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(mapper='.bins')

molecular_dynamics.RadialDistributionFunction.value.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(mapper='.value')

# Custom Correlation Functions mapping using get_custom_ensemble_outputs
MolecularDynamicsResults.correlation_functions.m_annotations.setdefault('mapping', {})[
    'hdf5'
] = MapperAnnotation(
    mapper=(
        'get_custom_ensemble_outputs',
        ['observables'],
        dict(
            exclude=STANDARD_H5MD_OBSERVABLES,
            type_filter=['correlation_function'],
        ),
    )
)

# Individual field mappings for custom CorrelationFunction
add_mapping_annotation(CorrelationFunction.label, HDF5_KEY, '.label')

add_mapping_annotation(CorrelationFunction.times, HDF5_KEY, '.times')

add_mapping_annotation(
    CorrelationFunction.value_magnitude, HDF5_KEY, '.value_magnitude'
)

add_mapping_annotation(CorrelationFunction.value_unit, HDF5_KEY, '.value_unit')

add_mapping_annotation(CorrelationFunction.direction, HDF5_KEY, '.direction')

add_mapping_annotation(CorrelationFunction.n_times, HDF5_KEY, '.n_times')

# Mean Squared Displacements mapping using get_output_data
MolecularDynamicsResults.mean_squared_displacements.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(
    mapper=(
        'get_output_data',
        ['observables.mean_squared_displacements'],
    )
)

# Individual field mappings for MeanSquaredDisplacement
molecular_dynamics.MeanSquaredDisplacement.label.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(mapper='.label')

molecular_dynamics.MeanSquaredDisplacement.times.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(mapper='.times')

molecular_dynamics.MeanSquaredDisplacement.value.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(mapper='.value')

molecular_dynamics.MeanSquaredDisplacement.direction.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(mapper='.direction')

molecular_dynamics.MeanSquaredDisplacement.n_times.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(mapper='.n_times')

# TODO add error quantities
# molecular_dynamics.MeanSquaredDisplacement.errors.m_annotations.setdefault(
#     'mapping', {}
# )['hdf5'] = MapperAnnotation(mapper='.errors')

# Diffusion Constants mapping using get_output_data
MolecularDynamicsResults.diffusion_constants.m_annotations.setdefault('mapping', {})[
    'hdf5'
] = MapperAnnotation(
    mapper=(
        'get_output_data',
        ['observables.diffusion_constants'],
    )
)

# Individual field mappings for DiffusionConstant
molecular_dynamics.DiffusionConstant.label.m_annotations.setdefault('mapping', {})[
    'hdf5'
] = MapperAnnotation(mapper='.label')

molecular_dynamics.DiffusionConstant.value.m_annotations.setdefault('mapping', {})[
    'hdf5'
] = MapperAnnotation(mapper='.value')

MolecularDynamicsResults.radii_of_gyration.m_annotations.setdefault('mapping', {})[
    'hdf5'
] = MapperAnnotation(mapper='.@')

# ! multi-ensemble property!
MolecularDynamicsResults.free_energy_calculations.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(mapper='.@')


# Use our custom MolecularDynamics class with custom results
add_mapping_annotation(MolecularDynamics.m_def, HDF5_KEY, '@')

add_mapping_annotation(MolecularDynamics.method, HDF5_KEY, '@')

add_mapping_annotation(MolecularDynamics.results, HDF5_KEY, '@')

add_mapping_annotation(molecular_dynamics.MolecularDynamics.outputs, HDF5_KEY, '@')

try:
    m_package.__init_metainfo__()
except Exception:
    pass

# TODO Check parameters for enums and add enum_spec
