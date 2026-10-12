import numpy as np
import pytest

from tests.parsers.common import SimulationParserTestSuite


class PhonopyParserIntegrationSuite(SimulationParserTestSuite):
    expected_program_name = 'Phonopy'
    required_simulation_sections = ('model_system',)


@pytest.mark.large_fixture
class TestVaspPhonopyArchive(PhonopyParserIntegrationSuite):
    archive_fixture = 'vasp_phonopy_archive'

    @pytest.mark.integration
    def test_identity_populated_once(self, archive):
        # The unit cell and supercell are distinct structures, not trajectory
        # frames, so both sections legitimately carry particle identities.
        assert all(system.particle_states for system in archive.data.model_system)

    @pytest.mark.integration
    def test_archive_contract(self, archive):
        simulation = archive.data

        assert len(simulation.model_system) == 2
        assert [
            state.chemical_symbol
            for state in simulation.model_system[0].particle_states
        ][:3] == ['C', 'C', 'C']
        assert simulation.model_system[0].positions.shape == (48, 3)
        assert simulation.model_system[1].positions.shape == (384, 3)
        np.testing.assert_array_equal(
            simulation.model_system[1].representations[0].supercell_matrix,
            np.diag([2, 2, 2]),
        )


class TestNoncanonicalHexagonalArchive(PhonopyParserIntegrationSuite):
    archive_fixture = 'cp2k_hexagonal_noncanonical_phonopy_archive'

    @pytest.mark.integration
    def test_identity_populated_once(self, archive):
        # Unit cell and supercell are distinct structures, not trajectory
        # frames, so both sections legitimately carry particle identities.
        assert all(system.particle_states for system in archive.data.model_system)

    @pytest.mark.integration
    def test_archive_contract(self, archive):
        simulation = archive.data

        assert len(simulation.model_system) == 2
        assert simulation.model_system[0].positions.shape == (288, 3)
        assert simulation.model_system[1].positions.shape == (288, 3)
        np.testing.assert_array_equal(
            simulation.model_system[1].representations[0].supercell_matrix,
            np.eye(3, dtype=int),
        )

    @pytest.mark.integration
    def test_mapping_populates_method_and_outputs(self, archive):
        simulation = archive.data
        assert len(simulation.model_method) == 1
        assert simulation.model_method[0].name == 'harmonic lattice dynamics'
        assert len(simulation.outputs) == 1
        outputs = simulation.outputs[0]
        assert len(outputs.force_constants) == 1
        assert len(outputs.phonon_band_structures) == 9
        assert len(outputs.phonon_dos) == 1
        assert len(outputs.vibrational_free_energies) == 1
        assert outputs.model_system_ref is simulation.model_system[1]
