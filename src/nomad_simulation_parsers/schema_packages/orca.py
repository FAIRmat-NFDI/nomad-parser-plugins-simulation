from nomad.metainfo import SchemaPackage
from nomad_simulations.schema_packages import (
    atoms_state,
    basis_set,
    general,
    model_method,
    model_system,
    numerical_settings,
    outputs,
)
from nomad_simulations.schema_packages.properties import molecular_orbitals

from nomad_simulation_parsers.schema_packages.utils import add_mapping_annotation

m_package = SchemaPackage()

OUT_KEY = 'orca_out'

add_mapping_annotation(general.Simulation.m_def, OUT_KEY, '@')

add_mapping_annotation(general.Simulation.program, OUT_KEY, '.@')
add_mapping_annotation(general.Program.version, OUT_KEY, '.program_version')

add_mapping_annotation(general.Simulation.model_system, OUT_KEY, ('get_atoms', ['.@']))
add_mapping_annotation(
    model_system.ModelSystem.is_representative,
    OUT_KEY,
    '.is_representative',
)
add_mapping_annotation(model_system.ModelSystem.positions, OUT_KEY, '.positions')
add_mapping_annotation(model_system.ModelSystem.total_charge, OUT_KEY, '.total_charge')
add_mapping_annotation(
    model_system.ModelSystem.total_spin_multiplicity,
    OUT_KEY,
    '.total_spin_multiplicity',
)
add_mapping_annotation(atoms_state.AtomsState.m_def, OUT_KEY, '.particle_states')
add_mapping_annotation(
    atoms_state.AtomsState.chemical_symbol,
    OUT_KEY,
    '.chemical_symbol',
)

############# DFT ###################

add_mapping_annotation(
    model_method.DFT.m_def,
    OUT_KEY,
    ('get_dft_methods', ['.single_point.self_consistent.scf_settings']),
    update_mode='append',
)

add_mapping_annotation(
    model_method.DFT.xc,
    OUT_KEY,
    '.xc',
)
add_mapping_annotation(
    model_method.DFT.reference_form,
    OUT_KEY,
    '.reference_form',
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

add_mapping_annotation(
    model_method.HF.m_def,
    OUT_KEY,
    ('get_hf_methods', ['.single_point.ci']),
    update_mode='append',
)
add_mapping_annotation(model_method.HF.reference_form, OUT_KEY, '.reference_form')

############# Multireference (CAS) ###################

add_mapping_annotation(
    model_method.MultireferenceSCF.m_def,
    OUT_KEY,
    ('get_multireference_scf_methods', ['.single_point.casscf', '.input_file']),
    update_mode='append',
)
add_mapping_annotation(model_method.MultireferenceSCF.type, OUT_KEY, '.type')
add_mapping_annotation(
    model_method.MultireferenceSCF.state_treatment, OUT_KEY, '.state_treatment'
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
    '.orbital_space_type',
)

add_mapping_annotation(
    model_method.MultireferencePT.m_def,
    OUT_KEY,
    ('get_multireference_pt_methods', ['.single_point.casscf', '.input_file']),
    update_mode='append',
)
add_mapping_annotation(model_method.MultireferencePT.type, OUT_KEY, '.type')
add_mapping_annotation(model_method.MultireferencePT.name, OUT_KEY, '.name')
add_mapping_annotation(model_method.ActiveSpace.m_def, OUT_KEY, '.active_space')
add_mapping_annotation(
    model_method.ActiveSpace.n_active_electrons, OUT_KEY, '.n_active_electrons'
)
add_mapping_annotation(
    model_method.ActiveSpace.n_active_orbitals, OUT_KEY, '.n_active_orbitals'
)

add_mapping_annotation(
    model_method.MultireferenceCI.m_def,
    OUT_KEY,
    ('get_multireference_ci_methods', ['.single_point.casscf', '.input_file']),
    update_mode='append',
)
add_mapping_annotation(model_method.MultireferenceCI.type, OUT_KEY, '.type')
add_mapping_annotation(
    model_method.MultireferenceCI.state_treatment, OUT_KEY, '.state_treatment'
)
add_mapping_annotation(
    model_method.MultireferenceCI.n_state_groups, OUT_KEY, '.n_state_groups'
)
add_mapping_annotation(
    model_method.MultireferenceCI.state_multiplicities,
    OUT_KEY,
    '.state_multiplicities',
)
add_mapping_annotation(
    model_method.MultireferenceCI.n_roots_per_multiplicity,
    OUT_KEY,
    '.n_roots_per_multiplicity',
)
add_mapping_annotation(
    model_method.MultireferenceCI.state_weights, OUT_KEY, '.state_weights'
)
add_mapping_annotation(
    model_method.MultireferenceCI.active_space, OUT_KEY, '.active_space'
)

############# Local CC ###################

add_mapping_annotation(
    model_method.OrbitalLocalization.m_def,
    OUT_KEY,
    ('get_orbital_localization_methods', ['.single_point.loc']),
)
add_mapping_annotation(model_method.OrbitalLocalization.method, OUT_KEY, '.method')
add_mapping_annotation(
    model_method.OrbitalLocalization.n_localized_orbitals,
    OUT_KEY,
    '.n_localized_orbitals',
)

add_mapping_annotation(
    model_method.CC.m_def,
    OUT_KEY,
    ('get_coupled_cluster_methods', ['.single_point.ci', '.input_file']),
)
add_mapping_annotation(model_method.CC.type, OUT_KEY, '.type')
add_mapping_annotation(model_method.CC.excitation_order, OUT_KEY, '.excitation_order')
add_mapping_annotation(
    model_method.CC.perturbative_correction, OUT_KEY, '.perturbative_correction'
)
add_mapping_annotation(
    model_method.CC.perturbative_correction_order,
    OUT_KEY,
    '.perturbative_correction_order',
)
add_mapping_annotation(
    model_method.CC.explicit_correlation, OUT_KEY, '.explicit_correlation'
)
add_mapping_annotation(model_method.CC.local_correlation, OUT_KEY, '.local_correlation')

add_mapping_annotation(model_method.LocalCorrelation.type, OUT_KEY, '.type')
add_mapping_annotation(model_method.LocalCorrelation.spaces, OUT_KEY, '.spaces')

add_mapping_annotation(
    model_method.LocalCorrelationSpace.space_kind, OUT_KEY, '.space_kind'
)
add_mapping_annotation(
    model_method.LocalCorrelationSpace.occupied_tuple_kind,
    OUT_KEY,
    '.occupied_tuple_kind',
)
add_mapping_annotation(
    model_method.LocalCorrelationSpace.virtual_space_type,
    OUT_KEY,
    '.virtual_space_type',
)
add_mapping_annotation(
    model_method.LocalCorrelationSpace.excitation_order,
    OUT_KEY,
    '.excitation_order',
)

add_mapping_annotation(
    model_method.ModelMethod.numerical_settings,
    OUT_KEY,
    '.numerical_settings',
    m_def=numerical_settings.LocalCorrelationSettings.m_def,
)
add_mapping_annotation(
    numerical_settings.LocalCorrelationSettings.screening_thresholds,
    OUT_KEY,
    '.screening_thresholds',
)
add_mapping_annotation(
    numerical_settings.LocalCorrelationThreshold.name, OUT_KEY, '.name'
)
add_mapping_annotation(
    numerical_settings.LocalCorrelationThreshold.value, OUT_KEY, '.value'
)
add_mapping_annotation(
    numerical_settings.LocalCorrelationThreshold.applies_to,
    OUT_KEY,
    '.applies_to',
)

add_mapping_annotation(model_method.PerturbationMethod.type, OUT_KEY, '.type')
add_mapping_annotation(model_method.PerturbationMethod.order, OUT_KEY, '.order')
add_mapping_annotation(
    model_method.PerturbationMethod.m_def,
    OUT_KEY,
    ('get_perturbation_methods', ['.single_point.ci', '.single_point.mp2', '.input_file']),
    update_mode='append',
)

add_mapping_annotation(
    model_method.LocalCorrelation.m_def, OUT_KEY, '.local_correlation'
)
add_mapping_annotation(model_method.LocalCorrelation.type, OUT_KEY, '.type')

add_mapping_annotation(basis_set.AtomCenteredBasisSet.basis_set, OUT_KEY, '.basis_set')
add_mapping_annotation(basis_set.AtomCenteredBasisSet.type, OUT_KEY, '.type')
add_mapping_annotation(basis_set.AtomCenteredBasisSet.role, OUT_KEY, '.role')
add_mapping_annotation(
    basis_set.AtomCenteredBasisSet.n_total_basis_functions,
    OUT_KEY,
    '.n_total_basis_functions',
)
add_mapping_annotation(
    basis_set.AtomCenteredBasisSet.m_def, OUT_KEY, ('get_basis_set_components', ['.@'])
)
add_mapping_annotation(
    basis_set.BasisSetContainer.m_def, OUT_KEY, '.@', update_mode='append'
)
add_mapping_annotation(
    basis_set.BasisSetContainer.native_tier, OUT_KEY, '.main_basis_set'
)


add_mapping_annotation(general.Simulation.outputs, OUT_KEY, ('get_outputs', ['.@']))
add_mapping_annotation(outputs.Outputs.model_system_ref, OUT_KEY, '.model_system_ref')
add_mapping_annotation(
    molecular_orbitals.MolecularOrbitals.m_def,
    OUT_KEY,
    '.molecular_orbitals',
)
add_mapping_annotation(molecular_orbitals.MolecularOrbitals.n_mo, OUT_KEY, '.n_mo')
add_mapping_annotation(molecular_orbitals.MolecularOrbitals.n_ao, OUT_KEY, '.n_ao')
add_mapping_annotation(
    molecular_orbitals.MolecularOrbitals.occupations,
    OUT_KEY,
    '.occupations',
)
# `value` overrides the `PhysicalProperty.value` placeholder; its class-attribute
# access can resolve to `None`, so read the quantity from the `m_def` registry.
add_mapping_annotation(
    molecular_orbitals.MolecularOrbitals.m_def.all_quantities['value'],
    OUT_KEY,
    '.energies',
)
add_mapping_annotation(
    molecular_orbitals.MolecularOrbitals.coefficients,
    OUT_KEY,
    '.coefficients',
)
add_mapping_annotation(
    molecular_orbitals.MolecularOrbitals.kind,
    OUT_KEY,
    '.kind',
)

add_mapping_annotation(model_method.RelativityModel.level, OUT_KEY, '.level')
add_mapping_annotation(
    model_method.RelativityModel.approximation, OUT_KEY, '.approximation'
)
add_mapping_annotation(model_method.RelativityModel.dkh_order, OUT_KEY, '.dkh_order')
add_mapping_annotation(
    model_method.ModelMethod.contributions,
    OUT_KEY,
    ('get_relativity_models', ['@']),
    cache=True,
    m_def=model_method.RelativityModel.m_def,
)

try:
    m_package.__init_metainfo__()
except Exception:
    pass
