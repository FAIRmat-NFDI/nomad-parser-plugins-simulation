from nomad.metainfo import SchemaPackage
from nomad_simulations.schema_packages import (
    atoms_state,
    general,
    model_system,
    outputs,
)
from nomad_simulations.schema_packages.properties import molecular_orbitals

from nomad_simulation_parsers.schema_packages.utils import add_mapping_annotation

m_package = SchemaPackage()

OUT_KEY = 'out'

add_mapping_annotation(general.Simulation.m_def, OUT_KEY, '@')

add_mapping_annotation(general.Simulation.program, OUT_KEY, '.@')
add_mapping_annotation(general.Program.name, OUT_KEY, ('get_program_name', ['.@']))
add_mapping_annotation(general.Program.version, OUT_KEY, '.program_version')

add_mapping_annotation(general.Simulation.model_system, OUT_KEY, ('get_atoms', ['.@']))
add_mapping_annotation(
    model_system.ModelSystem.is_representative,
    OUT_KEY,
    '.is_representative',
)
add_mapping_annotation(model_system.ModelSystem.positions, OUT_KEY, '.positions')
add_mapping_annotation(model_system.ModelSystem.total_charge, OUT_KEY, '.total_charge')
add_mapping_annotation(model_system.ModelSystem.total_spin, OUT_KEY, '.total_spin')
add_mapping_annotation(atoms_state.AtomsState.m_def, OUT_KEY, '.particle_states')
add_mapping_annotation(
    atoms_state.AtomsState.chemical_symbol,
    OUT_KEY,
    '.chemical_symbol',
)

add_mapping_annotation(general.Simulation.outputs, OUT_KEY, ('get_outputs', ['.@']))
add_mapping_annotation(outputs.Outputs.model_system_ref, OUT_KEY, '.model_system_ref')
add_mapping_annotation(
    molecular_orbitals.MolecularOrbitals.m_def,
    OUT_KEY,
    '.electronic_eigenvalues',
)
add_mapping_annotation(molecular_orbitals.MolecularOrbitals.n_mo, OUT_KEY, '.n_mo')
add_mapping_annotation(molecular_orbitals.MolecularOrbitals.n_ao, OUT_KEY, '.n_ao')
add_mapping_annotation(
    molecular_orbitals.MolecularOrbitals.mo_occupations,
    OUT_KEY,
    '.mo_occupations',
)
add_mapping_annotation(
    molecular_orbitals.MolecularOrbitals.mo_energies,
    OUT_KEY,
    '.mo_energies',
)
add_mapping_annotation(
    molecular_orbitals.MolecularOrbitals.mo_type,
    OUT_KEY,
    '.mo_type',
)

try:
    m_package.__init_metainfo__()
except Exception:
    pass
