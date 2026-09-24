import numpy as np
from nomad.client import normalize_all
from nomad.datamodel import EntryArchive, EntryMetadata
from nomad.utils import get_logger

from tests.parsers.common import (
    SimulationParserPipelineTestSuite,
    approx,
    assert_approx,
    assert_identity_populated_once,
)

LOGGER = get_logger(__name__)


def normalize_archive(archive):
    normalized = EntryArchive.m_from_dict(archive.m_to_dict())
    normalized.metadata = EntryMetadata()
    normalize_all(normalized, logger=LOGGER)
    return normalized


class TestLammpsPipeline(SimulationParserPipelineTestSuite):
    archive_fixture = 'methyl_naphthalene_archive'
    expected_program_name = 'LAMMPS'

    def test_methyl_naphthalene_normalized_values(self, methyl_naphthalene_archive):
        normalized = normalize_archive(methyl_naphthalene_archive)
        systems = normalized.data.model_system

        assert len(systems) == 4
        assert np.shape(systems[0].positions) == (1134, 3)
        assert systems[0].n_particles == 1134
        assert systems[0].particle_states[100].chemical_symbol == 'H'
        assert systems[0].particle_states[100].label == 'H'
        assert_identity_populated_once(normalized)
        assert systems[2].positions[567][1].to('angstrom').magnitude == approx(
            -5.88475
        )
        assert systems[3].lattice_vectors[2][2].to(
            'angstrom'
        ).magnitude == approx(21.468)
        assert systems[3].periodic_boundary_conditions == [True, True, True]
        assert_approx(systems[0].bond_list[200], np.array([189, 192]))
        assert systems[0].dimensionality == 3
        assert systems[0].is_molecule() is False

    def test_methyl_naphthalene_molecular_hierarchy(self, methyl_naphthalene_archive):
        normalized = normalize_archive(methyl_naphthalene_archive)
        groups = normalized.data.model_system[0].sub_systems

        assert len(groups) == 1
        assert groups[0].particle_states == []
        assert groups[0].name == 'group_0'
        assert groups[0].branch_label == 'molecule_group'
        assert groups[0].composition_formula == '0(54)'
        assert groups[0].particle_indices[13] == 13
        assert groups[0].is_molecule() is False

        monomers = groups[0].sub_systems
        assert len(monomers) == 54
        assert monomers[0].name == '0'
        assert monomers[0].branch_label == 'molecule'
        assert monomers[0].composition_formula == 'C(11)H(10)'
        assert monomers[0].particle_indices[20] == 20
        assert monomers[0].is_molecule() is True


class TestLammpsXyzPipeline(SimulationParserPipelineTestSuite):
    archive_fixture = 'xyz_archive'
    expected_program_name = 'LAMMPS'

    def test_xyz_normalized_velocities(self, xyz_archive):
        normalized = normalize_archive(xyz_archive)
        systems = normalized.data.model_system

        assert np.shape(systems[0].velocities) == (500, 3)
        assert systems[100].velocities[250][2].to(
            'angstrom/ps'
        ).magnitude == approx(0.0256726)
