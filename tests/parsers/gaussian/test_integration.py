import pytest

from tests.parsers.common import SimulationParserTestSuite, approx, assert_approx


class GaussianParserIntegrationSuite(SimulationParserTestSuite):
    expected_program_name = 'Gaussian'


class TestAlSCFArchive(GaussianParserIntegrationSuite):
    archive_fixture = 'al_scf_archive'

    @pytest.mark.integration
    def test_archive_contract(self, archive):
        simulation = archive.data

        assert simulation.program.version.startswith('09 EM64L-G09RevB.01')
        assert len(simulation.model_method) == 1
        assert simulation.model_method[0].m_def.name == 'DFT'
        assert simulation.model_method[0].xc.functional_key == 'B3LYP'
        assert len(simulation.model_system) == 1
        system = simulation.model_system[0]
        assert system.positions.shape == (1, 3)
        assert_approx(system.positions.to('angstrom').magnitude, [[0.0, 0.0, 0.0]])
        assert system.lattice_vectors is None
        assert system.particle_states[0].atomic_number == 13
        assert system.particle_states[0].chemical_symbol == 'Al'
        assert len(simulation.outputs) == 1
        assert simulation.outputs[0].total_energies[0].value.to(
            'hartree'
        ).magnitude == approx(-242.389700911934)


class TestAlMultistepArchive(GaussianParserIntegrationSuite):
    archive_fixture = 'al_multistep_archive'

    @pytest.mark.integration
    def test_archive_contract(self, archive):
        simulation = archive.data

        assert len(simulation.model_method) == 2
        assert all(
            method.xc.functional_key == 'PBEPBE' for method in simulation.model_method
        )
        assert len(simulation.model_system) == 7
        assert len(simulation.outputs) == 6
        system = simulation.model_system[0]
        assert system.positions.shape == (2, 3)
        assert_approx(
            system.positions.to('angstrom').magnitude,
            [[0.0, 0.0, 1.375], [0.0, 0.0, -1.375]],
        )
        assert system.lattice_vectors is None
        assert [
            (state.atomic_number, state.chemical_symbol)
            for state in system.particle_states
        ] == [(13, 'Al'), (13, 'Al')]
        assert all(not frame.particle_states for frame in simulation.model_system[1:])
        assert simulation.outputs[4].total_energies[0].value.to(
            'hartree'
        ).magnitude == approx(-484.282281400868)
        assert simulation.outputs[4].total_forces[0].value.shape == (2, 3)


class TestCHOFrequencyArchive(GaussianParserIntegrationSuite):
    archive_fixture = 'cho_frequency_archive'

    @pytest.mark.integration
    def test_archive_contract(self, archive):
        simulation = archive.data

        assert len(simulation.model_method) == 2
        assert all(
            method.xc.functional_key == 'SVWN' for method in simulation.model_method
        )
        assert len(simulation.model_system) == 25
        assert len(simulation.outputs) == 24
        system = simulation.model_system[0]
        assert system.positions.shape == (13, 3)
        assert_approx(
            system.positions.to('angstrom').magnitude[[0, -1]],
            [[0.672383, -2.128358, 0.069305], [1.595194, 1.205311, 0.059439]],
        )
        assert system.lattice_vectors is None
        assert [
            (state.atomic_number, state.chemical_symbol)
            for state in system.particle_states[:4]
        ] == [
            (1, 'H'),
            (7, 'N'),
            (8, 'O'),
            (8, 'O'),
        ]
        assert all(not frame.particle_states for frame in simulation.model_system[1:])
        assert simulation.outputs[-1].total_energies[0].value.to(
            'hartree'
        ).magnitude == approx(-322.091639810912)


class TestNOMPArchive(GaussianParserIntegrationSuite):
    archive_fixture = 'no_mp_archive'
    required_simulation_sections = ('model_system', 'outputs')

    @pytest.mark.integration
    def test_archive_contract(self, archive):
        simulation = archive.data

        assert simulation.program.version.startswith('09 EM64L-G09RevD.01')
        assert len(simulation.model_method) == 0
        assert len(simulation.model_system) == 1
        system = simulation.model_system[0]
        assert [
            (state.atomic_number, state.chemical_symbol)
            for state in system.particle_states
        ] == [(7, 'N'), (7, 'N'), (8, 'O'), (8, 'O')]
        assert system.positions.shape == (4, 3)
        assert_approx(
            system.positions.to('angstrom').magnitude,
            [
                [0.0, 1.1782, 0.624158],
                [0.0, -1.1782, 0.624158],
                [0.0, 1.263496, -0.546138],
                [0.0, -1.263496, -0.546138],
            ],
        )
        assert system.lattice_vectors is None
        assert len(simulation.outputs) == 1
        assert simulation.outputs[0].total_energies[0].value.to(
            'hartree'
        ).magnitude == approx(-258.103194703929)
