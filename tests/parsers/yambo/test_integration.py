import pytest

from tests.parsers.common import SimulationParserTestSuite, approx, assert_approx


class TestYamboHBNArchive(SimulationParserTestSuite):
    archive_fixture = 'hbn_archive'
    expected_program_name = 'YAMBO'
    required_simulation_sections = ('model_method', 'model_system', 'outputs')

    @pytest.mark.integration
    def test_exposes_band_outputs(self, archive):
        assert archive.data.outputs
        assert archive.data.outputs[0].electronic_band_structures
        assert archive.data.outputs[0].electronic_eigenvalues

    @pytest.mark.integration
    def test_maps_structure_and_band_assertions(self, archive):
        system = archive.data.model_system[0]
        eigenvalues = archive.data.outputs[0].electronic_eigenvalues[0]
        gap = archive.data.outputs[0].electronic_band_gaps[0]
        kpoints = archive.data.model_method[0].numerical_settings[0].k_mesh[0]

        assert archive.data.program.version == '5.0.4 Revision 19598'
        assert system.particle_states[3].chemical_symbol == 'N'
        assert_approx(
            system.lattice_vectors.to('bohr').magnitude,
            [
                [4.716000080108643, -2.3580000400543213, 0.0],
                [0.0, 4.0841755867004395, 0.0],
                [0.0, 0.0, 12.176712036132812],
            ],
        )
        assert_approx(kpoints.all_points[6][1], -0.49999997)
        assert gap.value.to('eV').magnitude == approx(3.878048)
        assert eigenvalues.value.shape == (14, 100)
        assert eigenvalues.value[6][7].to('eV').magnitude == approx(0.17268175)


class TestYamboLiFArchive(SimulationParserTestSuite):
    archive_fixture = 'lif_archive'
    expected_program_name = 'YAMBO'
    required_simulation_sections = ('model_method', 'outputs')

    @pytest.mark.integration
    def test_maps_qp_eigenvalue_archive(self, archive):
        eigenvalues = archive.data.outputs[0].electronic_eigenvalues[0]

        assert archive.data.program.version == '4.4.0 Revision 148'
        assert eigenvalues.value.shape == (10, 40)


class TestYamboGaSbArchive(SimulationParserTestSuite):
    archive_fixture = 'gasb_archive'
    expected_program_name = 'YAMBO'
    required_simulation_sections = ('model_method', 'outputs')

    @pytest.mark.integration
    def test_maps_gw_eigenvalue_archive(self, archive):
        eigenvalues = archive.data.outputs[0].electronic_eigenvalues[0]

        assert archive.data.program.version == '4.0.0 Revision 4245'
        assert eigenvalues.value.shape == (29, 40)
        assert eigenvalues.value[18][5].to('eV').magnitude == approx(-2.61243)


class TestYamboAluminumArchive(SimulationParserTestSuite):
    archive_fixture = 'aluminum_archive'
    expected_program_name = 'YAMBO'
    required_simulation_sections = ('model_method', 'outputs')

    @pytest.mark.integration
    def test_maps_lifetime_eigenvalue_archive(self, archive):
        output = archive.data.outputs[0]
        eigenvalues = output.electronic_eigenvalues[0]

        assert archive.data.program.version == '4.4.0 Revision 148'
        assert output.electronic_band_structures[0].highest_occupied.to(
            'eV'
        ).magnitude == approx(8.363973)
        assert eigenvalues.value.shape == (29, 16)
        assert eigenvalues.value[16][9].to('eV').magnitude == approx(26.26274)


class TestYamboCH4Archive(SimulationParserTestSuite):
    archive_fixture = 'ch4_archive'
    expected_program_name = 'YAMBO'
    required_simulation_sections = ('model_method', 'model_system', 'outputs')

    @pytest.mark.integration
    def test_maps_setup_structure(self, archive):
        system = archive.data.model_system[0]

        assert archive.data.program.version == '5.0.0 Revision 18656'
        assert len(system.positions) == 5
        assert system.particle_states[0].chemical_symbol == 'H'
        assert system.particle_states[4].chemical_symbol == 'C'
        assert_approx(
            system.positions.to('bohr').magnitude[0],
            [1.2755651, 1.2755651, 1.2755651],
            rtol=1e-6,
        )
        assert_approx(
            system.lattice_vectors.to('bohr').magnitude,
            [
                [9.44863, 0.0, 0.0],
                [0.0, 9.44863, 0.0],
                [0.0, 0.0, 9.44863],
            ],
        )
