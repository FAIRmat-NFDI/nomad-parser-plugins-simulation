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
add_mapping_annotation(
    model_method.XCFunctional.functional_key,
    OUT_KEY,
    '.functional_key',
)

############# Multireference (CAS) ###################

add_mapping_annotation(
    model_method.MultireferenceSCF.m_def,
    OUT_KEY,
    ('get_multireference_methods', ['.@']),
)
add_mapping_annotation(model_method.MultireferenceSCF.type, OUT_KEY, '.type')
add_mapping_annotation(
    model_method.MultireferenceSCF.reference_type, OUT_KEY, '.reference_type'
)
add_mapping_annotation(
    model_method.MultireferenceSCF.n_state_groups, OUT_KEY, '.n_state_groups'
)
add_mapping_annotation(
    model_method.MultireferenceSCF.state_multiplicities,
    OUT_KEY,
    '.state_multiplicities',
)
add_mapping_annotation(
    model_method.MultireferenceSCF.n_roots_per_multiplicity,
    OUT_KEY,
    '.n_roots_per_multiplicity',
)
add_mapping_annotation(
    model_method.MultireferenceSCF.state_weights, OUT_KEY, '.state_weights'
)
add_mapping_annotation(
    model_method.MultireferenceSCF.active_space, OUT_KEY, '.active_space'
)
add_mapping_annotation(
    model_method.ActiveSpace.n_active_electrons,
    OUT_KEY,
    '.active_space.n_active_electrons',
)
add_mapping_annotation(
    model_method.ActiveSpace.n_active_orbitals,
    OUT_KEY,
    '.active_space.n_active_orbitals',
)
add_mapping_annotation(
    model_method.ActiveSpace.orbital_space_type,
    OUT_KEY,
    '.active_space.orbital_space_type',
)


try:
    m_package.__init_metainfo__()
except Exception:
    pass
