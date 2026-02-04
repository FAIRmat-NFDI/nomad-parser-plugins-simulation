from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass

from nomad.metainfo import SchemaPackage
from nomad_simulations.schema_packages import (
    atoms_state,
    general,
    model_system,
)

from nomad_simulation_parsers.schema_packages.utils import add_mapping_annotation

OUT_KEY = 'out'

m_package = SchemaPackage()

add_mapping_annotation(general.Simulation.m_def, OUT_KEY, '@')


class Simulation(general.Simulation):
    add_mapping_annotation(
        general.Simulation.program, OUT_KEY, ('get_program_data', ['.@'])
    )
    add_mapping_annotation(
        general.Simulation.model_system, OUT_KEY, ('get_systems', ['.@'])
    )
    # add_mapping_annotation(
    #     general.Simulation.outputs, OUT_KEY, ('get_outputs', ['.@'])
    # )


class Program(general.Program):
    add_mapping_annotation(general.Program.name, OUT_KEY, '.name')
    add_mapping_annotation(general.Program.version, OUT_KEY, '.version')


class ModelSystem(model_system.ModelSystem):
    add_mapping_annotation(model_system.ModelSystem.positions, OUT_KEY, '.positions')
    add_mapping_annotation(
        model_system.ModelSystem.total_charge, OUT_KEY, '.total_charge'
    )
    add_mapping_annotation(model_system.ModelSystem.total_spin, OUT_KEY, '.total_spin')
    add_mapping_annotation(atoms_state.AtomsState.m_def, OUT_KEY, '.particle_states')


class AtomsState(atoms_state.AtomsState):
    add_mapping_annotation(
        atoms_state.AtomsState.chemical_symbol, OUT_KEY, '.chemical_symbol'
    )


# class Outputs(outputs.Outputs):
#     add_mapping_annotation(outputs.Outputs.total_energies, OUT_KEY, '.total_energies')
#     add_mapping_annotation(outputs.Outputs.total_forces, OUT_KEY, '.total_forces')


# class TotalEnergy(properties.energies.TotalEnergy):
#     add_mapping_annotation(
#         properties.energies.TotalEnergy.value,
#         OUT_KEY,
#         '.value || .energy_total || .@',
#     )


# class TotalForce(properties.forces.TotalForce):
#     add_mapping_annotation(
#         properties.forces.TotalForce.value, OUT_KEY, '.value || .forces || .@'
#     )

try:
    m_package.__init_metainfo__()
except Exception:
    pass
