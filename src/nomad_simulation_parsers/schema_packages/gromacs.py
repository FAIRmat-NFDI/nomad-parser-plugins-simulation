from nomad.metainfo import SchemaPackage
from nomad_simulations.schema_packages import (
    general,
    model_system,
    outputs,
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


class AtomicCell(model_system.AtomicCell):
    add_mapping_annotation(
        model_system.AtomicCell.lattice_vectors, TPR_KEY, '.lattice_vectors'
    )
    add_mapping_annotation(
        model_system.AtomicCell.periodic_boundary_conditions, LOG_KEY, '.pbc'
    )


class ModelSystem(model_system.ModelSystem):
    add_mapping_annotation(model_system.ModelSystem.velocities, TPR_KEY, '.velocities')
    add_mapping_annotation(model_system.ModelSystem.positions, TPR_KEY, '.positions')
    add_mapping_annotation(model_system.AtomsState.m_def, TPR_KEY, '.labels')
    add_mapping_annotation(model_system.AtomicCell.m_def, LOG_KEY, '.@')
    add_mapping_annotation(model_system.AtomicCell.m_def, TPR_KEY, '.@')


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


try:
    m_package.__init_metainfo__()
except Exception:
    pass
