from nomad.metainfo import SchemaPackage
from nomad_simulations.schema_packages import (
    general,
    model_system,
    outputs,
)
from nomad_simulations.schema_packages.workflow import (
    geometry_optimization,
    molecular_dynamics,
)

from nomad_simulation_parsers.schema_packages.utils import add_mapping_annotation

m_package = SchemaPackage()


LOG_KEY = 'log'
TPR_KEY = 'tpr'
EDR_KEY = 'edr'


class Program(general.Program):
    add_mapping_annotation(
        general.Program.version, LOG_KEY, ('get_version', ['.version'])
    )


class AtomsState(model_system.AtomsState):
    add_mapping_annotation(model_system.AtomsState.label, TPR_KEY, '.@')


class AtomicCell(model_system.Representation):
    """
    Map the representation quantities used by GROMACS to the unified
    Representation fields so annotations and conversion remain correct.
    """

    add_mapping_annotation(
        model_system.Representation.lattice_vectors, TPR_KEY, '.lattice_vectors'
    )
    add_mapping_annotation(
        model_system.Representation.periodic_boundary_conditions, LOG_KEY, '.pbc'
    )


class ModelSystem(model_system.ModelSystem):
    add_mapping_annotation(model_system.ModelSystem.velocities, TPR_KEY, '.velocities')
    add_mapping_annotation(model_system.ModelSystem.positions, TPR_KEY, '.positions')
    add_mapping_annotation(model_system.AtomsState.m_def, TPR_KEY, '.labels')
    add_mapping_annotation(model_system.Representation.m_def, LOG_KEY, '.@')
    add_mapping_annotation(model_system.Representation.m_def, TPR_KEY, '.@')


class TotalEnergy(outputs.TotalEnergy):
    add_mapping_annotation(outputs.TotalEnergy.name, LOG_KEY, '.label')
    add_mapping_annotation(outputs.TotalEnergy.name, EDR_KEY, '.label')
    add_mapping_annotation(outputs.TotalEnergy.value, LOG_KEY, '.value')
    add_mapping_annotation(outputs.TotalEnergy.value, EDR_KEY, '.value')
    add_mapping_annotation(outputs.TotalEnergy.contributions, LOG_KEY, '.contributions')
    add_mapping_annotation(outputs.TotalEnergy.contributions, EDR_KEY, '.contributions')


class TotalForce(outputs.TotalForce):
    add_mapping_annotation(outputs.TotalForce.value, TPR_KEY, '.@')


class Outpus(outputs.Outputs):
    add_mapping_annotation(outputs.Outputs.total_energies, LOG_KEY, '.energy')
    add_mapping_annotation(outputs.Outputs.total_energies, EDR_KEY, '.energy')
    add_mapping_annotation(outputs.Outputs.model_system_ref, LOG_KEY, '.system_ref')
    add_mapping_annotation(outputs.Outputs.model_system_ref, EDR_KEY, '.system_ref')
    add_mapping_annotation(outputs.Outputs.total_forces, TPR_KEY, '.forces')


class Simulation(general.Simulation):
    add_mapping_annotation(general.Simulation.program, LOG_KEY, '.header')
    add_mapping_annotation(
        general.Simulation.model_system, LOG_KEY, ('get_configurations', [])
    )
    add_mapping_annotation(
        general.Simulation.model_system, TPR_KEY, ('get_configurations', [])
    )
    add_mapping_annotation(general.Simulation.outputs, LOG_KEY, ('get_outputs', []))
    add_mapping_annotation(general.Simulation.outputs, TPR_KEY, ('get_outputs', []))
    add_mapping_annotation(general.Simulation.outputs, EDR_KEY, ('get_outputs', ['.@']))


add_mapping_annotation(general.Simulation.m_def, LOG_KEY, '@')
add_mapping_annotation(general.Simulation.m_def, TPR_KEY, '@')
add_mapping_annotation(general.Simulation.m_def, EDR_KEY, '@')


class GeometryOptimizationModel(geometry_optimization.GeometryOptimization):
    add_mapping_annotation(
        geometry_optimization.GeometryOptimizationModel.optimization_method,
        LOG_KEY,
        ('get_integrator_type', ['.input_parameters.integrator']),
    )
    add_mapping_annotation(
        geometry_optimization.GeometryOptimizationModel.n_steps_maximum,
        LOG_KEY,
        '.input_parameters.nsteps',
    )
    add_mapping_annotation(
        geometry_optimization.GeometryOptimizationModel.convergence_tolerance_force_maximum,
        LOG_KEY,
        '.input_parameters.emtol',
        unit='kilojoule/avogadro_number/nanometer',
    )


class GeometryOptimizationResults(geometry_optimization.GeometryOptimizationResults):
    add_mapping_annotation(
        geometry_optimization.GeometryOptimizationResults.energies,
        LOG_KEY,
        ('get_energies', []),
    )
    add_mapping_annotation(
        geometry_optimization.GeometryOptimizationResults.energies,
        EDR_KEY,
        ('get_energies', []),
    )
    add_mapping_annotation(
        geometry_optimization.GeometryOptimizationResults.final_force_maximum,
        LOG_KEY,
        'maximum_force',
        unit='kilojoule/avogadro_number/nanometer',
    )


class GeometryOptimization(geometry_optimization.GeometryOptimization):
    add_mapping_annotation(
        geometry_optimization.GeometryOptimizationModel.m_def, LOG_KEY, '.@'
    )
    add_mapping_annotation(
        geometry_optimization.GeometryOptimizationResults.m_def, LOG_KEY, '.@'
    )
    add_mapping_annotation(
        geometry_optimization.GeometryOptimizationResults.m_def, EDR_KEY, '.@'
    )


class MolecularDynamicsModel(molecular_dynamics.MolecularDynamicsModel):
    add_mapping_annotation(
        molecular_dynamics.MolecularDynamicsModel.integrator_type,
        LOG_KEY,
        ('get_integrator_type', ['.input_parameters.integrator']),
    )
    add_mapping_annotation(
        molecular_dynamics.MolecularDynamicsModel.integration_timestep,
        LOG_KEY,
        '.input_parameters.dt',
        unit='picosecond',
    )


class MolecularDynamicsResults(molecular_dynamics.MolecularDynamicsResults):
    # parse from xvg
    pass


class MolecularDynamics(molecular_dynamics.MolecularDynamics):
    add_mapping_annotation(
        molecular_dynamics.MolecularDynamicsModel.m_def, LOG_KEY, '.@'
    )


# Workflow
add_mapping_annotation(geometry_optimization.GeometryOptimization.m_def, LOG_KEY, '@')
add_mapping_annotation(geometry_optimization.GeometryOptimization.m_def, EDR_KEY, '@')
add_mapping_annotation(molecular_dynamics.MolecularDynamics.m_def, LOG_KEY, '@')


try:
    m_package.__init_metainfo__()
except Exception:
    pass
