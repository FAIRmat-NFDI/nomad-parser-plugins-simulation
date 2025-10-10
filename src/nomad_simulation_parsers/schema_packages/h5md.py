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
from nomad.metainfo.data_type import m_float64
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

m_package = SchemaPackage()


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

    name.m_annotations.setdefault('mapping', {})['hdf5'] = MapperAnnotation(
        mapper='."@name"'
    )

    email = Quantity(
        type=str,
        shape=[],
        description="""
        Author's email.
        """,
    )

    email.m_annotations.setdefault('mapping', {})['hdf5'] = MapperAnnotation(
        mapper='."@email"'
    )


## class Program(general.Program):


general.Program.name.m_annotations.setdefault('mapping', {})['hdf5'] = MapperAnnotation(
    mapper='."@name"',
)

general.Program.version.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(
        mapper='."@version"',
    )
)

# SIMULATION.MODEL_SYSTEM --> archive.data.model_system

# ! Only use generic ParticleState
# ! ModelSystem normalizer takes care of assigning to AtomState or CGBeadState
### class ParticleState(atoms_state.ParticleState):

atoms_state.ParticleState.m_def.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(
        mapper=('to_species_labels', ['.@'], dict(path='particles.all.species_label'))
    )
)

atoms_state.ParticleState.label.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(mapper='.label')
)

### class AtomicCell(model_system.AtomicCell):

model_system.AtomicCell.m_def.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(mapper=('get_cell_data', ['.@']))
)

model_system.AtomicCell.lattice_vectors.m_annotations.setdefault('mapping', {})[
    'hdf5'
] = MapperAnnotation(mapper='.lattice_vectors')

model_system.AtomicCell.periodic_boundary_conditions.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(mapper='.boundary')


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
# ModelSystem.custom_system_attributes.m_annotations.setdefault('mapping', {})[
#     'hdf5'
# ] = MapperAnnotation(
#     mapper=(
#         'get_custom_system_attributes',
#         ['.@'],
#         dict(path=''),
#     )
# )

ModelSystem.n_particles.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(mapper='.n_particles')
)


ModelSystem.bond_list.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(
        mapper=('get_top_system_quantity', ['.@'], dict(path='connectivity.bonds'))
    )
)

ModelSystem.dimensionality.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(
        mapper=(
            'get_top_system_quantity',
            ['.@'],
            dict(path='particles.all.box.@dimension'),
        )
    )
)

ModelSystem.positions.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(mapper='.positions')
)

ModelSystem.velocities.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(mapper='.velocities')
)

ModelSystem.sub_systems.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(mapper=('get_sub_systems', ['.@'], dict(path='connectivity')))
)

ModelSystem.name.m_annotations.setdefault('mapping', {})['hdf5'] = MapperAnnotation(
    mapper='.label'
)

ModelSystem.composition_formula.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(mapper='.formula')
)

ModelSystem.particle_indices.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(mapper='.indices')
)

ModelSystem.branch_label.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(mapper='.type')
)

# ModelSystem.type.m_annotations.setdefault('mapping', {})['hdf5'] = MapperAnnotation(
#     mapper='.type'
# ) # TODO add function to filter out valid Enums


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

    value = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description="""
        Value **magnitude** of the property.
        The unit is defined in the `unit` attribute.
        """,
    )
    value.m_annotations.setdefault('mapping', {})['hdf5'] = MapperAnnotation(
        mapper='.value'
    )

    unit = Quantity(
        type=str,
        shape=[],
        description="""
        Unit of the parameter as a string consistent with the UnitRegistry.pint module.
        """,
    )
    unit.m_annotations.setdefault('mapping', {})['hdf5'] = MapperAnnotation(
        mapper='.unit'
    )

    description = Quantity(
        type=str,
        shape=[],
        description="""
        Further description of the property.
        """,
    )
    description.m_annotations.setdefault('mapping', {})['hdf5'] = MapperAnnotation(
        mapper='.description'
    )

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
CustomProperty.name.m_annotations.setdefault('mapping', {})['hdf5'] = MapperAnnotation(
    mapper='.name'
)


## class BaseEnergy(properties.energies.BaseEnergy):

# ? Does this note still apply?
# value annotation defined in TotalEnergy.value since they refer to the same quantity
# in this case, we make sure to return the corresponding value from
# the get_contributions function in the TotalEnergy.contributions annotation
properties.energies.BaseEnergy.contribution_type.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(mapper='.name')


## class TotalEnergy(properties.TotalEnergy):


properties.TotalEnergy.value.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(
        mapper=(
            'get_output_data',
            ['.@'],
            dict(path='observables.energies.total'),
        )
    )
)

### SUBSECTIONS

properties.TotalEnergy.contributions.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(
        mapper=(
            'get_contributions',
            ['.@'],
            dict(path='observables.energies', exclude=['total']),
        )
    )
)


## class BaseForce(properties.forces.BaseForce):


properties.forces.BaseForce.contribution_type.m_annotations.setdefault('mapping', {})[
    'hdf5'
] = MapperAnnotation(mapper='.name')


## class TotalForce(properties.TotalForce):


properties.TotalForce.value.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(
        mapper=(
            'get_output_data',
            ['.@'],
            dict(path='observables.forces.total'),
        )
    )
)

### SUBSECTIONS

properties.TotalForce.contributions.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(
        mapper=(
            'get_contributions',
            ['.@'],
            dict(path='observables.forces', exclude=['total']),
        )
    )
)


## class Temperature(properties.Temperature):


properties.Temperature.value.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(
        mapper=(
            'get_output_data',
            ['.@'],
            dict(path='observables.temperatures'),
        )
    )
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


TrajectoryOutputs.m_def.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(mapper=('get_output_steps', ['observables']))
)

TrajectoryOutputs.step.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(mapper='.step')
)

TrajectoryOutputs.time.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(mapper='.time')
)

### SUBSECTIONS

TrajectoryOutputs.total_energies.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(mapper='.@')
)

TrajectoryOutputs.total_forces.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(mapper='.@')
)

TrajectoryOutputs.temperatures.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(mapper='.@')
)

TrajectoryOutputs.custom_outputs.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(
        mapper=(
            'get_custom_outputs',
            ['.@'],
            dict(
                path='observables',
                exclude=STANDARD_H5MD_OBSERVABLES,
            ),
        )
    )
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
    x_h5md_version.m_annotations.setdefault('mapping', {})['hdf5'] = MapperAnnotation(
        mapper='h5md."@version"',
    )

    x_h5md_author = SubSection(sub_section=Author.m_def)

    x_h5md_author.m_annotations.setdefault('mapping', {})['hdf5'] = MapperAnnotation(
        mapper='h5md.author'
    )

    x_h5md_creator = SubSection(sub_section=general.Program.m_def)

    x_h5md_creator.m_annotations.setdefault('mapping', {})['hdf5'] = MapperAnnotation(
        mapper='h5md.creator'
    )

    model_system = SubSection(sub_section=ModelSystem.m_def, repeats=True)


Simulation.m_def.m_annotations.setdefault('mapping', {})['hdf5'] = MapperAnnotation(
    mapper='@'
)

Simulation.program.m_annotations.setdefault('mapping', {})['hdf5'] = MapperAnnotation(
    mapper='h5md.program'
)

Simulation.model_system.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(mapper=('get_traj_data', ['particles.all']))
)


# WORKFLOW --> archive.workflow2
h5md_path_md = 'parameters.workflow.molecular_dynamics'

# WORKFLOW.METHOD --> archive.workflow2.method


## class ThermostatParameters(molecular_dynamics.ThermostatParameters):


h5md_path_thermostat = f'{h5md_path_md}.thermostat_parameters'

molecular_dynamics.ThermostatParameters.thermostat_type.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(
    mapper=(
        'map_value',
        [h5md_path_thermostat],
        dict(key='thermostat_type', enum_spec='lower'),
    )
)

molecular_dynamics.ThermostatParameters.reference_temperature.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(
    mapper=(
        'map_value',
        [h5md_path_thermostat],
        dict(key='reference_temperature'),
    )
)

molecular_dynamics.ThermostatParameters.coupling_constant.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(
    mapper=(
        'map_value',
        [h5md_path_thermostat],
        dict(key='coupling_constant'),
    )
)

molecular_dynamics.ThermostatParameters.effective_mass.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(
    mapper=(
        'map_value',
        [h5md_path_thermostat],
        dict(key='effective_mass'),
    )
)

molecular_dynamics.ThermostatParameters.temperature_profile.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(
    mapper=(
        'map_value',
        [h5md_path_thermostat],
        dict(key='temperature_profile', enum_spec='lower'),
    )
)

molecular_dynamics.ThermostatParameters.reference_temperature_start.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(
    mapper=(
        'map_value',
        [h5md_path_thermostat],
        dict(key='reference_temperature_start'),
    )
)

molecular_dynamics.ThermostatParameters.reference_temperature_end.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(
    mapper=(
        'map_value',
        [h5md_path_thermostat],
        dict(key='reference_temperature_end'),
    )
)

molecular_dynamics.ThermostatParameters.temperature_update_frequency.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(
    mapper=(
        'map_value',
        [h5md_path_thermostat],
        dict(key='temperature_update_frequency'),
    )
)

molecular_dynamics.ThermostatParameters.temperature_update_delta.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(
    mapper=(
        'map_value',
        [h5md_path_thermostat],
        dict(key='temperature_update_delta'),
    )
)

molecular_dynamics.ThermostatParameters.temperature_update_factor.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(
    mapper=(
        'map_value',
        [h5md_path_thermostat],
        dict(key='temperature_update_factor'),
    )
)

molecular_dynamics.ThermostatParameters.step_start.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(
    mapper=(
        'map_value',
        [h5md_path_thermostat],
        dict(key='step_start'),
    )
)

molecular_dynamics.ThermostatParameters.step_end.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(
    mapper=(
        'map_value',
        [h5md_path_thermostat],
        dict(key='step_end'),
    )
)


## class BarostatParameters(molecular_dynamics.BarostatParameters):


h5md_path_barostat = f'{h5md_path_md}.barostat_parameters'

molecular_dynamics.BarostatParameters.barostat_type.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(
    mapper=(
        'map_value',
        [h5md_path_barostat],
        dict(key='barostat_type', enum_spec='lower'),
    )
)

molecular_dynamics.BarostatParameters.coupling_type.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(
    mapper=(
        'map_value',
        [h5md_path_barostat],
        dict(key='coupling_type', enum_spec='lower'),
    )
)

molecular_dynamics.BarostatParameters.reference_pressure.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(
    mapper=(
        'map_value',
        [h5md_path_barostat],
        dict(key='reference_pressure'),
    )
)

molecular_dynamics.BarostatParameters.coupling_constant.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(
    mapper=(
        'map_value',
        [h5md_path_barostat],
        dict(key='coupling_constant'),
    )
)

molecular_dynamics.BarostatParameters.compressibility.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(
    mapper=(
        'map_value',
        [h5md_path_barostat],
        dict(key='compressibility'),
    )
)

molecular_dynamics.BarostatParameters.pressure_profile.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(
    mapper=(
        'map_value',
        [h5md_path_barostat],
        dict(key='pressure_profile', enum_spec='lower'),
    )
)

molecular_dynamics.BarostatParameters.reference_pressure_start.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(
    mapper=(
        'map_value',
        [h5md_path_barostat],
        dict(key='reference_pressure_start'),
    )
)

molecular_dynamics.BarostatParameters.reference_pressure_end.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(
    mapper=(
        'map_value',
        [h5md_path_barostat],
        dict(key='reference_pressure_end'),
    )
)

molecular_dynamics.BarostatParameters.pressure_update_frequency.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(
    mapper=(
        'map_value',
        [h5md_path_barostat],
        dict(key='pressure_update_frequency'),
    )
)

molecular_dynamics.BarostatParameters.pressure_update_delta.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(
    mapper=(
        'map_value',
        [h5md_path_barostat],
        dict(key='pressure_update_delta'),
    )
)

molecular_dynamics.BarostatParameters.pressure_update_factor.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(
    mapper=(
        'map_value',
        [h5md_path_barostat],
        dict(key='pressure_update_factor'),
    )
)

molecular_dynamics.BarostatParameters.step_start.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(
    mapper=(
        'map_value',
        [h5md_path_barostat],
        dict(key='step_start'),
    )
)

molecular_dynamics.BarostatParameters.step_end.m_annotations.setdefault('mapping', {})[
    'hdf5'
] = MapperAnnotation(
    mapper=(
        'map_value',
        [h5md_path_barostat],
        dict(key='step_end'),
    )
)


## class ShearParameters(molecular_dynamics.ShearParameters):


h5md_path_shear = f'{h5md_path_md}.shear_parameters'

molecular_dynamics.ShearParameters.shear_type.m_annotations.setdefault('mapping', {})[
    'hdf5'
] = MapperAnnotation(
    mapper=(
        'map_value',
        [h5md_path_shear],
        dict(key='shear_type', enum_spec='lower'),
    )
)

molecular_dynamics.ShearParameters.shear_rate.m_annotations.setdefault('mapping', {})[
    'hdf5'
] = MapperAnnotation(
    mapper=(
        'map_value',
        [h5md_path_shear],
        dict(key='shear_rate'),
    )
)

molecular_dynamics.ShearParameters.step_start.m_annotations.setdefault('mapping', {})[
    'hdf5'
] = MapperAnnotation(
    mapper=(
        'map_value',
        [h5md_path_shear],
        dict(key='step_start'),
    )
)

molecular_dynamics.ShearParameters.step_end.m_annotations.setdefault('mapping', {})[
    'hdf5'
] = MapperAnnotation(
    mapper=(
        'map_value',
        [h5md_path_shear],
        dict(key='step_end'),
    )
)


## class FreeEnergyCalculationParameters(
##     molecular_dynamics.FreeEnergyCalculationParameters
## ):


h5md_path_FEC = f'{h5md_path_md}.free_energy_calculation_parameters'

# TODO Change this to fec_type in the schema
molecular_dynamics.FreeEnergyCalculationParameters.calc_type.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(
    mapper=(
        'map_value',
        [h5md_path_FEC],
        dict(key='type', enum_spec='lower'),
    )
)

molecular_dynamics.FreeEnergyCalculationParameters.current_lambda_index.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(
    mapper=(
        'map_value',
        [h5md_path_FEC],
        dict(key='lambda_index'),
    )
)

# molecular_dynamics.FreeEnergyCalculationParameters.atom_indices.m_annotations.\
#     setdefault(
#     'mapping', {}
# )['hdf5'] = MapperAnnotation(
#     mapper=(
#         'map_value',
#         [h5md_path_FEC],
#         dict(key='atom_indices'),
#     )
# )

# molecular_dynamics.FreeEnergyCalculationParameters.initial_state_vdw.m_annotations.\
#     setdefault(
#     'mapping', {}
# )['hdf5'] = MapperAnnotation(
#     mapper=(
#         'map_value',
#         [h5md_path_FEC],
#         dict(key='initial_state_vdw'),
#     )
# )

# molecular_dynamics.FreeEnergyCalculationParameters.final_state_vdw.m_annotations.\
#     setdefault(
#     'mapping', {}
# )['hdf5'] = MapperAnnotation(
#     mapper=(
#         'map_value',
#         [h5md_path_FEC],
#         dict(key='final_state_vdw'),
#     )
# )

# molecular_dynamics.FreeEnergyCalculationParameters.initial_state_coloumb.\
#     m_annotations.setdefault(
#     'mapping', {}
# )['hdf5'] = MapperAnnotation(
#     mapper=(
#         'map_value',
#         [h5md_path_FEC],
#         dict(key='initial_state_coloumb'),
#     )
# )

# molecular_dynamics.FreeEnergyCalculationParameters.final_state_coloumb.\
#     m_annotations.setdefault(
#     'mapping', {}
# )['hdf5'] = MapperAnnotation(
#     mapper=(
#         'map_value',
#         [h5md_path_FEC],
#         dict(key='final_state_coloumb'),
#     )
# )

# molecular_dynamics.FreeEnergyCalculationParameters.initial_state_bonded.\
#     m_annotations.setdefault(
#     'mapping', {}
# )['hdf5'] = MapperAnnotation(
#     mapper=(
#         'map_value',
#         [h5md_path_FEC],
#         dict(key='initial_state_bonded'),
#     )
# )

# molecular_dynamics.FreeEnergyCalculationParameters.final_state_bonded.\
#     m_annotations.setdefault(
#     'mapping', {}
# )['hdf5'] = MapperAnnotation(
#     mapper=(
#         'map_value',
#         [h5md_path_FEC],
#         dict(key='final_state_bonded'),
#     )
# )

### SUBSECTIONS

molecular_dynamics.FreeEnergyCalculationParameters.lambdas.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(mapper='@')


## class Lambdas(molecular_dynamics.Lambdas):
# TODO add lambdas to test data

# ? Not sure about where this info is going in h5md
h5md_path_lambdas = f'{h5md_path_FEC}.lambdas'

molecular_dynamics.Lambdas.interaction_type.m_annotations.setdefault('mapping', {})[
    'hdf5'
] = MapperAnnotation(
    mapper=(
        'map_value',
        [h5md_path_lambdas],
        dict(key='type', enum_spec='lower'),
    )
)

molecular_dynamics.Lambdas.values.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(
        mapper=(
            'map_value',
            [h5md_path_lambdas],
            dict(key='value'),
        )
    )
)


## class MolecularDynamicsMethod(molecular_dynamics.MolecularDynamicsMethod):


molecular_dynamics.MolecularDynamicsMethod.thermodynamic_ensemble.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(
    mapper=(
        'map_value',
        [h5md_path_md],
        dict(key='thermodynamic_ensemble', enum_spec='upper'),
    )
)

molecular_dynamics.MolecularDynamicsMethod.integrator_type.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(
    mapper=(
        'map_value',
        [h5md_path_md],
        dict(key='integrator_type', enum_spec='lower'),
    )
)

molecular_dynamics.MolecularDynamicsMethod.integration_timestep.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(
    mapper=(
        'map_value',
        [h5md_path_md],
        dict(key='integration_timestep'),
    )
)

molecular_dynamics.MolecularDynamicsMethod.n_steps.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(
    mapper=(
        'map_value',
        [h5md_path_md],
        dict(key='n_steps'),
    )
)

molecular_dynamics.MolecularDynamicsMethod.coordinate_save_frequency.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(
    mapper=(
        'map_value',
        [h5md_path_md],
        dict(key='coordinate_save_frequency'),
    )
)

molecular_dynamics.MolecularDynamicsMethod.velocity_save_frequency.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(
    mapper=(
        'map_value',
        [h5md_path_md],
        dict(key='velocity_save_frequency'),
    )
)

molecular_dynamics.MolecularDynamicsMethod.force_save_frequency.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(
    mapper=(
        'map_value',
        [h5md_path_md],
        dict(key='force_save_frequency'),
    )
)

molecular_dynamics.MolecularDynamicsMethod.thermodynamics_save_frequency.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(
    mapper=(
        'map_value',
        [h5md_path_md],
        dict(key='thermodynamics_save_frequency'),
    )
)

### SUBSECTIONS

molecular_dynamics.MolecularDynamicsMethod.thermostat_parameters.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(mapper='@')

molecular_dynamics.MolecularDynamicsMethod.barostat_parameters.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(mapper='@')

molecular_dynamics.MolecularDynamicsMethod.shear_parameters.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(mapper='@')

molecular_dynamics.MolecularDynamicsMethod.free_energy_calculation_parameters.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(mapper='@')

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


class MolecularDynamics(molecular_dynamics.MolecularDynamics):
    results = SubSection(sub_section=MolecularDynamicsResults.m_def)


# MolecularDynamicsOutputs.m_def.m_annotations.setdefault('mapping', {})['hdf5'] = (
#     MapperAnnotation(mapper=('get_output_data', ['observables']))
# )

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
EnsembleProperty.label.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(mapper='.label')
)

EnsembleProperty.bins_magnitude.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(mapper='.bins')
)

EnsembleProperty.bins_unit.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(mapper='.bins_unit')
)

EnsembleProperty.value_magnitude.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(mapper='.value_magnitude')
)

EnsembleProperty.value_unit.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(mapper='.value_unit')
)

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
CorrelationFunction.label.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(mapper='.label')
)

CorrelationFunction.times.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(mapper='.times')
)

CorrelationFunction.value_magnitude.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(mapper='.value_magnitude')
)

CorrelationFunction.value_unit.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(mapper='.value_unit')
)

CorrelationFunction.direction.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(mapper='.direction')
)

CorrelationFunction.n_times.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(mapper='.n_times')
)

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

# ? Needed? It just points to the trajectory properties? I guess it collects data here?
MolecularDynamicsResults.radii_of_gyration.m_annotations.setdefault('mapping', {})[
    'hdf5'
] = MapperAnnotation(mapper='.@')

# ! multi-ensemble property!
MolecularDynamicsResults.free_energy_calculations.m_annotations.setdefault(
    'mapping', {}
)['hdf5'] = MapperAnnotation(mapper='.@')


# Use our custom MolecularDynamics class with custom results
MolecularDynamics.m_def.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(mapper='@')
)

MolecularDynamics.method.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(mapper='@')
)

MolecularDynamics.results.m_annotations.setdefault('mapping', {})['hdf5'] = (
    MapperAnnotation(mapper='@')
)

# ? Needed?
molecular_dynamics.MolecularDynamics.outputs.m_annotations.setdefault('mapping', {})[
    'hdf5'
] = MapperAnnotation(mapper='@')
# MolecularDynamics.outputs.m_annotations.setdefault('mapping', {})['hdf5'] = (
#     MapperAnnotation(mapper=('get_output_data', ['observables']))
# )

try:
    m_package.__init_metainfo__()
except Exception:
    pass

# TODO Check parameters for enums and add enum_spec
