import numpy as np
from nomad.client import normalize_all
from nomad.datamodel import EntryArchive, EntryMetadata
from nomad.utils import get_logger

from nomad_simulation_parsers.parsers.lammps.parser import LammpsArchiveWriter
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

        assert normalized.workflow2.m_def.name == 'MolecularDynamics'
        assert normalized.workflow2.method.thermodynamic_ensemble == 'NPT'
        assert (
            normalized.workflow2.method.barostat_parameters[0].coupling_type
            == 'isotropic'
        )
        assert len(systems) == 4
        assert np.shape(systems[0].positions) == (1134, 3)
        assert systems[0].n_particles == 1134
        assert systems[0].particle_states[100].chemical_symbol == 'H'
        assert systems[0].particle_states[100].label == 'H'
        assert_identity_populated_once(normalized)
        assert systems[2].positions[567][1].to('angstrom').magnitude == approx(-5.88475)
        assert systems[3].lattice_vectors[2][2].to('angstrom').magnitude == approx(
            21.468
        )
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
        assert {monomer.name for monomer in monomers} == {'0'}
        assert len({tuple(m.particle_indices) for m in monomers}) == 54
        assert len(normalized.workflow2.results.radial_distribution_functions) == 1
        molecule_topology = next(
            entry
            for entry in normalized.results.material.topology
            if entry.label == '0'
        )
        assert len(molecule_topology.indices) == 54

    def test_hexane_molecule_names_identify_types(self, hexane_archive):
        groups = hexane_archive.data.model_system[0].sub_systems

        assert len(groups) == 2
        assert sorted(len(group.sub_systems) for group in groups) == [18, 18]
        assert {
            group.name: {molecule.name for molecule in group.sub_systems}
            for group in groups
        } == {'group_0': {'0'}, 'group_1': {'1'}}

        normalized = normalize_archive(hexane_archive)
        results = normalized.workflow2.results
        assert len(results.radial_distribution_functions) == 12
        assert len(results.mean_squared_displacements) == 2
        assert len(results.diffusion_constants) == 2
        topology = {
            entry.label: entry
            for entry in normalized.results.material.topology
            if entry.label in {'0', '1'}
        }
        assert len(topology['0'].indices) == 18
        assert len(topology['1'].indices) == 18

    def test_nested_molecule_names_stay_distinct(self):
        particle_arrays = {
            'molnums': np.repeat([0, 1], 6),
            'resids': np.repeat([0, 1, 2, 3], 3),
            'resnames': np.array(['X'] * 12),
        }
        writer = LammpsArchiveWriter()

        molecules = [
            writer._create_molecule(molecule, 'P', particle_arrays)
            for molecule in (0, 1)
        ]

        assert [molecule.name for molecule in molecules] == ['0', '1']
        assert all(len(molecule.sub_systems) == 1 for molecule in molecules)
        assert all(
            len(molecule.sub_systems[0].sub_systems) == 2 for molecule in molecules
        )

    def test_normalization_preserves_md_method(self, methyl_naphthalene_archive):
        normalized = normalize_archive(methyl_naphthalene_archive)
        assert normalized.workflow2.method.thermodynamic_ensemble == 'NPT'
        assert normalized.workflow2.method.n_steps == 30000
        assert normalized.workflow2.results.n_steps == len(normalized.data.model_system)
        assert len(normalized.workflow2.results.trajectory) == len(
            normalized.data.model_system
        )
