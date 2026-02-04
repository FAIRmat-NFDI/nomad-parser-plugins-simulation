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


############# Simulation + Program ###################

add_mapping_annotation(general.Simulation.m_def, OUT_KEY, '@')

add_mapping_annotation(
    general.Simulation.program, OUT_KEY, ('get_program_data', ['.@'])
)

add_mapping_annotation(general.Program.name, OUT_KEY, '.name')
add_mapping_annotation(general.Program.version, OUT_KEY, '.version')

############# Atoms / ModelSystem ###################

add_mapping_annotation(
    general.Simulation.model_system,
    OUT_KEY,
    ('get_atoms', ['.@']),
)

add_mapping_annotation(
    model_system.ModelSystem.positions,
    OUT_KEY,
    '.positions',
)

add_mapping_annotation(
    model_system.ModelSystem.particle_states,
    OUT_KEY,
    '.particle_states',
)

add_mapping_annotation(
    atoms_state.AtomsState.m_def,
    OUT_KEY,
    '.particle_states',
)

add_mapping_annotation(
    atoms_state.AtomsState.atomic_number,
    OUT_KEY,
    '.atomic_number',
)

try:
    m_package.__init_metainfo__()
except Exception:
    pass
