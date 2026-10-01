import numpy as np
import pytest
from nomad.utils import get_logger

from tests.parsers.common import (
    SimulationParserTestSuite,
    approx,
    assert_identity_populated_once,
)

LOGGER = get_logger(__name__)


class LammpsParserIntegrationSuite(SimulationParserTestSuite):
    expected_program_name = 'LAMMPS'
    required_simulation_sections = ('model_method', 'model_system')

    def test_representative_system_is_complete(self, archive):
        representative = archive.data.model_system[-1]
        assert representative.positions is not None


@pytest.mark.integration
class TestLammpsMethylNaphthaleneArchive(LammpsParserIntegrationSuite):
    archive_fixture = 'methyl_naphthalene_archive'

    def test_archive_contains_trajectory_and_topology(self, archive):
        simulation = archive.data

        assert len(simulation.model_system) == 4
        assert simulation.model_system[0].positions.shape == (1134, 3)
        assert simulation.model_system[0].bond_list.shape[1] == 2
        assert simulation.model_system[0].particle_states[100].chemical_symbol == 'H'
        assert_identity_populated_once(archive)

    def test_archive_has_valid_particle_indices(self, archive):
        system = archive.data.model_system[0]
        assert system.sub_systems
        for group in system.sub_systems:
            assert np.all(group.particle_indices >= 0)
            assert np.all(group.particle_indices < system.n_particles)


@pytest.mark.integration
class TestLammpsHexaneArchive(LammpsParserIntegrationSuite):
    archive_fixture = 'hexane_archive'

    def test_archive_contains_nvt_structure(self, archive):
        simulation = archive.data
        assert len(simulation.model_system) == 201
        assert simulation.model_method

        system = simulation.model_system[5]
        assert system.lattice_vectors[1][1].to('nanometer').magnitude == approx(2.24235)
        assert system.periodic_boundary_conditions == [True, True, True]
        assert simulation.model_system[80].particle_states == []
        assert simulation.model_system[0].particle_states[91].chemical_symbol == 'H'
        assert simulation.model_system[0].particle_states[94].chemical_symbol == 'C'
        assert simulation.model_system[0].bond_list[200, 0] == 194
        assert_identity_populated_once(archive)


@pytest.mark.integration
class TestLammpsXyzArchive(LammpsParserIntegrationSuite):
    archive_fixture = 'xyz_archive'

    def test_archive_contains_velocities(self, archive):
        system = archive.data.model_system[0]
        assert system.positions.shape == (500, 3)
        assert system.velocities.shape == (500, 3)
        assert archive.data.model_system[100].velocities[250][2].to(
            'angstrom/ps'
        ).magnitude == approx(0.0256726)
        assert archive.data.model_system[1].positions[452][2].to(
            'meter'
        ).magnitude == approx(5.99898)
        assert archive.data.model_system[2].velocities[457][-2].to(
            'meter/second'
        ).magnitude == approx(-0.928553)


@pytest.mark.integration
class TestLammpsPolymerMeltArchive(LammpsParserIntegrationSuite):
    archive_fixture = 'polymer_melt_minimization_archive'

    def test_archive_contains_minimization_trajectory(self, archive):
        simulation = archive.data
        assert len(simulation.model_system) == 159
        assert simulation.model_method

        topology = simulation.model_system[0]
        assert topology.positions.shape == (7200, 3)
        assert topology.n_particles == 7200
        lattice_length = topology.lattice_vectors[0][0].to('nanometer').magnitude
        assert lattice_length == approx(4.80158)
        assert topology.periodic_boundary_conditions == [True, True, True]
        assert topology.bond_list.shape == (7100, 2)
        assert topology.sub_systems
        assert simulation.model_system[1].particle_states == []
        assert simulation.model_system[0].positions[468][0].to(
            'angstrom'
        ).magnitude == approx(4.2341122)
        assert_identity_populated_once(archive)


@pytest.mark.integration
class TestLammpsDcdArchive(LammpsParserIntegrationSuite):
    archive_fixture = 'methane_dcd_archive'

    def test_reads_dcd_trajectory(self, archive):
        simulation = archive.data
        assert len(simulation.model_system) == 201

        topology = simulation.model_system[0]
        assert topology.positions.shape == (320, 3)
        assert topology.n_particles == 320
        assert len(topology.particle_states) == 320
        assert simulation.model_system[56].particle_states == []
        position = simulation.model_system[56].positions[7][0].to('meter').magnitude
        assert position == approx(-1.0751316e-10)
        assert_identity_populated_once(archive)


@pytest.mark.integration
class TestLammpsMethaneXyzArchive(LammpsParserIntegrationSuite):
    archive_fixture = 'methane_xyz_archive'

    def test_reads_xyz_trajectory(self, archive):
        simulation = archive.data
        assert len(simulation.model_system) == 201

        topology = simulation.model_system[0]
        assert topology.positions.shape == (320, 3)
        assert topology.n_particles == 320
        assert len(topology.particle_states) == 320
        assert simulation.model_system[13].particle_states == []
        position = simulation.model_system[13].positions[7][0].to('meter').magnitude
        assert position == approx(-8.00436e-10)
        assert_identity_populated_once(archive)
