from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass

from nomad.metainfo import SchemaPackage
from nomad_simulations.schema_packages import (
    atoms_state,
    general,
    model_method,
    model_system,
)

from nomad_simulation_parsers.schema_packages.utils import add_mapping_annotation

m_package = SchemaPackage()

OUT_KEY = 'out'


############# Simulation + Program ###################

add_mapping_annotation(general.Simulation.m_def, OUT_KEY, '@')

add_mapping_annotation(general.Simulation.program, OUT_KEY, '.@')

add_mapping_annotation(general.Program.version, OUT_KEY, '.program_version')

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
    atoms_state.AtomsState.m_def,
    OUT_KEY,
    '.particle_states',
)

add_mapping_annotation(
    atoms_state.AtomsState.chemical_symbol,
    OUT_KEY,
    '.chemical_symbol',
)

############# DFT ###################

add_mapping_annotation(
    model_method.DFT.m_def,
    OUT_KEY,
    ('get_dft', ['.@']),
)

add_mapping_annotation(
    model_method.DFT.xc,
    OUT_KEY,
    '.xc',
)

add_mapping_annotation(
    model_method.XCFunctional.global_exact_exchange,
    OUT_KEY,
    '.global_exact_exchange',
)


try:
    m_package.__init_metainfo__()
except Exception:
    pass
