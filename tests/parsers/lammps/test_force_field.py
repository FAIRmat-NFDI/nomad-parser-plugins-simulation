from types import SimpleNamespace
from unittest.mock import Mock

import numpy as np
import pytest
from nomad.client import normalize_all
from nomad.datamodel import EntryArchive, EntryMetadata
from nomad.units import ureg
from nomad.utils import get_logger
from nomad_simulations.schema_packages.atoms_state import ParticleState
from nomad_simulations.schema_packages.force_field import (
    ForceCalculations,
    ForceField,
    HarmonicAngle,
    HarmonicBond,
    ParticleParametersContainer,
    Potential,
)

from nomad_simulation_parsers.parsers.lammps.file_parsers import get_unit
from nomad_simulation_parsers.parsers.lammps.force_field import (
    attach_particle_parameters,
    build_force_field,
)
from nomad_simulation_parsers.parsers.lammps.parser import LammpsParser


def _reader(values, mainfile=None):
    return SimpleNamespace(mainfile=mainfile, get=values.get)


@pytest.mark.integration
def test_log_only_force_field_settings_are_retained_without_a_system(tmp_path):
    mainfile = tmp_path / 'log.lammps'
    mainfile.write_text(
        'LAMMPS (test)\nunits metal\npair_style lj/cut 8.0\npair_coeff 1 1 0.25 2.5\n'
    )
    archive = EntryArchive()
    LammpsParser().parse(str(mainfile), archive, get_logger(__name__))

    assert archive.data.model_system == []
    assert len(archive.data.model_method) == 1
    method = archive.data.model_method[0]
    assert isinstance(method, ForceField)
    assert method.numerical_settings[0].vdw_cutoff.to('angstrom').magnitude == (
        pytest.approx(8)
    )
    assert len(method.contributions) == 1
    parameters = {entry.name: entry for entry in method.contributions[0].parameters}
    assert float(parameters['epsilon'].value) == pytest.approx(
        (0.25 * ureg.eV).to('joule').magnitude
    )
    assert float(parameters['sigma'].value) == pytest.approx(2.5e-10)

    restored = EntryArchive.m_from_dict(archive.m_to_dict())
    restored.metadata = EntryMetadata()
    normalize_all(restored, logger=get_logger(__name__))
    assert isinstance(restored.data.model_method[0], ForceField)
    assert restored.data.model_method[0].numerical_settings[0].vdw_cutoff.to(
        'angstrom'
    ).magnitude == pytest.approx(8)
    assert restored.data.model_method[0].contributions[0].name == ('lj/cut types 1-1')


@pytest.mark.integration
def test_hexane_force_field_preserves_distinct_types_and_interactions(hexane_archive):
    simulation = hexane_archive.data
    method = simulation.model_method[0]
    assert isinstance(method, ForceField)
    settings = method.numerical_settings[0]
    assert isinstance(settings, ForceCalculations)
    assert settings.vdw_cutoff.to('angstrom').magnitude == pytest.approx(12)
    assert settings.coulomb_cutoff.to('angstrom').magnitude == pytest.approx(12)
    assert settings.coulomb_type == 'cutoff'
    # ``neigh_modify every 10 ... check yes`` is not a guaranteed rebuild rate.
    assert settings.neighbor_update_frequency is None

    bonds = [term for term in method.contributions if isinstance(term, HarmonicBond)]
    angles = [term for term in method.contributions if isinstance(term, HarmonicAngle)]
    dihedrals = [
        term for term in method.contributions if term.functional_form == 'OPLS'
    ]
    assert [term.parameters[0].value for term in bonds] == ['1', '2']
    assert [term.parameters[0].value for term in angles] == ['1', '2', '3']
    assert [len(term.particle_indices) for term in bonds] == [198, 468]
    assert [len(term.particle_indices) for term in angles] == [180, 828, 288]
    assert [len(term.particle_indices) for term in dihedrals] == [162, 756, 864]
    assert sum(len(term.particle_indices) for term in bonds) == 666
    assert sum(len(term.particle_indices) for term in angles) == 1296
    assert sum(len(term.particle_indices) for term in dihedrals) == 1782
    assert all(
        np.all(term.particle_indices >= 0)
        and np.all(term.particle_indices < simulation.model_system[0].n_particles)
        for term in [*bonds, *angles, *dihedrals]
    )
    # The schema uses 1/2 k Δ² whereas LAMMPS includes the 1/2 in its K.
    real_units = get_unit('real')
    bond_constant = (
        bonds[0].force_constant * real_units['distance'] ** 2 / real_units['energy']
    )
    angle_constant = angles[0].force_constant * ureg.radian**2 / real_units['energy']
    assert bond_constant.to('').magnitude == pytest.approx(536)
    assert angle_constant.to('').magnitude == pytest.approx(116.7)
    assert angles[0].equilibrium_value.to('degree').magnitude == pytest.approx(112.7)
    pairs = [
        term
        for term in method.contributions
        if type(term) is Potential and term.type == 'nonbonded'
    ]
    assert len(pairs) == 2
    assert pairs[0].n_interactions is None
    assert {parameter.name for parameter in pairs[0].parameters} == {
        'type_i',
        'type_j',
        'epsilon',
        'sigma',
    }
    parameters = next(
        setting
        for setting in method.numerical_settings
        if isinstance(setting, ParticleParametersContainer)
    )
    assert [len(entry.species_scope) for entry in parameters.particle_parameters] == [
        36,
        180,
        468,
    ]
    assert [entry.particle_type for entry in parameters.particle_parameters] == [
        '1',
        '1',
        '2',
    ]
    assert [
        entry.partial_charge.to('elementary_charge').magnitude
        for entry in parameters.particle_parameters
    ] == pytest.approx([-0.18, -0.12, 0.06])
    assert parameters.particle_parameters[0].effective_mass.to('amu').magnitude == (
        pytest.approx(12.011)
    )

    restored = EntryArchive.m_from_dict(hexane_archive.m_to_dict())
    restored_bonds = [
        term
        for term in restored.data.model_method[0].contributions
        if isinstance(term, HarmonicBond)
    ]
    assert restored_bonds[0].particle_indices.shape == (198, 2)
    restored_parameters = next(
        setting
        for setting in restored.data.model_method[0].numerical_settings
        if isinstance(setting, ParticleParametersContainer)
    )
    assert len(restored_parameters.particle_parameters[0].species_scope) == 36
    assert (
        restored_parameters.particle_parameters[0]
        .species_scope[0]
        .m_path()
        .startswith('/data/model_system/0/particle_states/')
    )


@pytest.mark.integration
@pytest.mark.parametrize('fixture_name', ['methane_dcd_archive', 'methane_xyz_archive'])
def test_data_backed_md_analysis_trajectories_retain_force_field_scopes(
    request, fixture_name
):
    archive = request.getfixturevalue(fixture_name)
    method = archive.data.model_method[0]
    assert isinstance(method, ForceField)
    bonds = [term for term in method.contributions if isinstance(term, HarmonicBond)]
    angles = [term for term in method.contributions if isinstance(term, HarmonicAngle)]
    assert len(bonds) == len(angles) == 1
    assert bonds[0].particle_indices.shape == (256, 2)
    # The DATA Atoms rows are unsorted; the reader's atom-ID order is 1..320.
    assert bonds[0].particle_indices[0].tolist() == [265, 266]
    assert angles[0].particle_indices.shape == (384, 3)
    particles = next(
        setting
        for setting in method.numerical_settings
        if isinstance(setting, ParticleParametersContainer)
    )
    assert [len(entry.species_scope) for entry in particles.particle_parameters] == [
        64,
        256,
    ]
    assert [
        entry.partial_charge.to('elementary_charge').magnitude
        for entry in particles.particle_parameters
    ] == pytest.approx([-0.24, 0.06])


@pytest.mark.integration
def test_normalization_preserves_typed_force_field_and_interaction_rows(
    methyl_naphthalene_archive,
):
    normalized = EntryArchive.m_from_dict(methyl_naphthalene_archive.m_to_dict())
    normalized.metadata = EntryMetadata()
    normalize_all(normalized, logger=get_logger(__name__))

    method = normalized.data.model_method[0]
    assert isinstance(method, ForceField)
    settings = next(
        item
        for item in method.numerical_settings
        if isinstance(item, ForceCalculations)
    )
    assert settings.vdw_cutoff.to('angstrom').magnitude == pytest.approx(12)
    assert settings.coulomb_type == 'particle_particle_particle_mesh'
    bonds = [term for term in method.contributions if isinstance(term, HarmonicBond)]
    assert len(bonds) == 6
    assert sum(term.n_interactions for term in bonds) == 1188
    assert bonds[0].particle_indices.shape == (270, 2)
    particle_parameters = next(
        item
        for item in method.numerical_settings
        if isinstance(item, ParticleParametersContainer)
    )
    assert (
        sum(len(item.species_scope) for item in particle_parameters.particle_parameters)
        == normalized.data.model_system[0].n_particles
    )


@pytest.mark.unit
def test_force_field_maps_atom_ids_and_log_overrides_data_coefficients():
    log = _reader(
        {
            'bond_style': ['harmonic'],
            'bond_coeff': [[2, 250.0, 1.2]],
            'neigh_modify': [['every', 5, 'check', 'no']],
        }
    )
    data = _reader(
        {
            'Bond Coeffs': [(None, np.asarray([[1, 100.0, 1.0], [2, 200.0, 1.1]]))],
            'Bonds': [(None, np.asarray([[1, 2, 30, 10], [2, 1, 10, 20]]))],
        },
        mainfile='data.test',
    )
    method = build_force_field(log, data, [10, 20, 30], get_unit('real'))
    assert isinstance(method, ForceField)
    bonds = [term for term in method.contributions if isinstance(term, HarmonicBond)]
    assert bonds[0].particle_indices.tolist() == [[0, 1]]
    assert bonds[1].particle_indices.tolist() == [[2, 0]]
    real_units = get_unit('real')
    overridden_constant = (
        bonds[1].force_constant * real_units['distance'] ** 2 / real_units['energy']
    )
    assert overridden_constant.to('').magnitude == pytest.approx(500)
    assert method.numerical_settings[0].neighbor_update_frequency == 5


@pytest.mark.unit
def test_force_field_does_not_guess_interactions_or_hybrid_cutoffs():
    log = _reader(
        {
            'bond_style': ['harmonic'],
            'bond_coeff': [[1, 100.0, 1.0]],
            'pair_style': [['hybrid', 'lj/cut', 12.0, 'coul/cut', 10.0]],
        }
    )
    data = _reader(
        {'Bonds': [(None, np.asarray([[1, 1, 1, 999]]))]},
        mainfile='data.test',
    )
    logger = Mock()
    method = build_force_field(log, data, [1, 2], get_unit('real'), logger)
    assert isinstance(method, ForceField)
    assert method.contributions[0].particle_indices is None
    assert method.numerical_settings == []
    logger.warning.assert_called_once()
    assert 'hybrid' in logger.warning.call_args.args[0]


@pytest.mark.unit
def test_explicit_pair_coefficients_and_cutoffs_keep_source_types():
    log = _reader(
        {
            'pair_style': [['lj/cut/coul/long', 9.0, 11.0]],
            'pair_coeff': [[1, 2, 0.5, 3.0, 8.0]],
            'kspace_style': [['pppm', 0.0001]],
            'neighbor': [[1.0, 'bin']],
        }
    )
    data = _reader(
        {'Pair Coeffs': [(None, np.asarray([[1, 0.1, 2.0], [2, 0.2, 2.5]]))]},
        mainfile='data.test',
    )
    method = build_force_field(log, data, None, get_unit('real'))
    assert isinstance(method, ForceField)
    settings = method.numerical_settings[0]
    assert settings.vdw_cutoff.to('angstrom').magnitude == pytest.approx(9)
    assert settings.coulomb_cutoff.to('angstrom').magnitude == pytest.approx(11)
    assert settings.neighbor_update_cutoff.to('angstrom').magnitude == pytest.approx(12)
    assert settings.coulomb_type == 'particle_particle_particle_mesh'
    pairs = method.contributions
    assert len(pairs) == 3
    assert [pair.name for pair in pairs] == [
        'lj/cut/coul/long types 1-1',
        'lj/cut/coul/long types 1-2',
        'lj/cut/coul/long types 2-2',
    ]
    middle = {parameter.name: parameter for parameter in pairs[1].parameters}
    assert float(middle['epsilon'].value) == pytest.approx(
        (0.5 * get_unit('real')['energy']).to('joule').magnitude
    )
    assert float(middle['vdw_cutoff'].value) == pytest.approx(8e-10)
    assert 'coulomb_cutoff' not in middle


@pytest.mark.unit
def test_neighbor_cutoff_includes_explicit_pair_override():
    log = _reader(
        {
            'pair_style': [['lj/cut/coul/cut', 9.0, 8.0]],
            'pair_coeff': [[1, 1, 0.5, 3.0, 12.0, 11.0]],
            'neighbor': [[1.0, 'bin']],
        }
    )
    method = build_force_field(log, _reader({}), None, get_unit('real'))
    assert isinstance(method, ForceField)
    settings = method.numerical_settings[0]
    assert settings.vdw_cutoff.to('angstrom').magnitude == pytest.approx(9)
    assert settings.coulomb_cutoff.to('angstrom').magnitude == pytest.approx(8)
    assert settings.neighbor_update_cutoff.to('angstrom').magnitude == (
        pytest.approx(13)
    )
    parameters = {item.name: item for item in method.contributions[0].parameters}
    assert float(parameters['vdw_cutoff'].value) == pytest.approx(1.2e-9)
    assert float(parameters['coulomb_cutoff'].value) == pytest.approx(1.1e-9)


@pytest.mark.unit
def test_opls_coefficients_remain_separate_and_source_scaled():
    log = _reader(
        {
            'dihedral_style': ['opls'],
            'dihedral_coeff': [[1, 1.74, -0.157, 0.279, 0.0]],
        }
    )
    data = _reader(
        {'Dihedrals': [(None, np.asarray([[1, 1, 40, 10, 20, 30]]))]},
        mainfile='data.test',
    )
    method = build_force_field(log, data, [10, 20, 30, 40], get_unit('real'))
    assert isinstance(method, ForceField)
    potential = method.contributions[0]
    assert potential.type == 'dihedral'
    assert potential.particle_indices.tolist() == [[3, 0, 1, 2]]
    parameters = {entry.name: entry for entry in potential.parameters}
    assert set(parameters) == {'lammps_type', 'K1', 'K2', 'K3', 'K4'}
    assert float(parameters['K1'].value) == pytest.approx(
        (1.74 * get_unit('real')['energy']).to('joule').magnitude
    )
    assert float(parameters['K2'].value) < 0


@pytest.mark.unit
def test_unresolved_pair_coefficients_are_diagnosed_not_invented():
    log = _reader(
        {
            'pair_style': [['lj/cut', 10.0]],
            'pair_coeff': [['*', '*', 0.5, 3.0]],
        }
    )
    logger = Mock()
    method = build_force_field(log, _reader({}), None, get_unit('real'), logger)
    assert isinstance(method, ForceField)
    assert method.contributions == []
    assert method.numerical_settings[0].vdw_cutoff is not None
    logger.warning.assert_called_once()
    assert 'pair coefficients' in logger.warning.call_args.args[0]


@pytest.mark.unit
def test_atomic_masses_require_exact_atom_id_match():
    log = _reader({'atom_style': ['atomic']})
    data = _reader(
        {
            'Atoms': [(None, np.asarray([[10, 1, 0, 0, 0], [20, 2, 1, 0, 0]]))],
            'Masses': [(None, np.asarray([[1, 12.0], [2, 1.0]]))],
        },
        mainfile='data.test',
    )
    states = [ParticleState(label='C'), ParticleState(label='H')]
    invalid = ForceField()
    attach_particle_parameters(invalid, log, data, [10, 30], states, get_unit('real'))
    assert invalid.numerical_settings == []

    valid = ForceField()
    attach_particle_parameters(valid, log, data, [20, 10], states, get_unit('real'))
    container = valid.numerical_settings[0]
    assert isinstance(container, ParticleParametersContainer)
    assert [entry.particle_type for entry in container.particle_parameters] == [
        '1',
        '2',
    ]
    assert container.particle_parameters[0].species_scope[0] is states[1]
    assert container.particle_parameters[1].species_scope[0] is states[0]
    assert all(entry.partial_charge is None for entry in container.particle_parameters)


@pytest.mark.unit
def test_neighbor_rebuild_controls_accumulate_without_inventing_a_rate():
    log = _reader(
        {
            'neigh_modify': [
                ['every', 5],
                ['check', 'no'],
                ['once', 'yes'],
            ]
        }
    )
    method = build_force_field(log, _reader({}), None, get_unit('real'))
    assert method is None

    log = _reader({'neigh_modify': [['every', 5], ['check', 'no']]})
    method = build_force_field(log, _reader({}), None, get_unit('real'))
    assert isinstance(method, ForceField)
    assert method.numerical_settings[0].neighbor_update_frequency == 5
