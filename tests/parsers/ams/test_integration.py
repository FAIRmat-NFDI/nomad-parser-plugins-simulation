import pytest

from tests.parsers.common import SimulationParserTestSuite, WorkflowTestSuite


def approx(value, abs=0, rel=1e-6):
    return pytest.approx(value, abs=abs, rel=rel)


class AMSParserIntegrationSuite(SimulationParserTestSuite, WorkflowTestSuite):
    expected_program_name = 'AMS'


class TestSCFArchive(AMSParserIntegrationSuite):
    archive_fixture = 'scf_archive'
    workflow_name = 'SinglePoint'

    @pytest.mark.integration
    def test_archive_contract(self, archive):
        errors, warnings = archive.m_validate()
        assert errors == []
        assert warnings == []
        simulation = archive.data
        assert simulation.program.version == '70630 2018-11-24'
        assert len(simulation.model_system) == 1
        assert len(simulation.model_system[0].particle_states) == 11
        assert [
            state.chemical_symbol
            for state in simulation.model_system[0].particle_states
        ] == ['C'] * 6 + ['H'] * 5
        assert simulation.model_system[0].positions.shape == (11, 3)
        assert simulation.model_system[0].positions.to('bohr').magnitude[
            2, 1
        ] == approx(4.36723014)
        output = simulation.outputs[0]
        assert output.total_energies[0].value.to('hartree').magnitude == approx(
            -2.55556188
        )
        assert simulation.model_method[0].xc.functional_key == 'TPSS'
        assert len(output.scf_steps.delta_energies_total) == 20
        assert output.scf_steps.delta_energies_total[-1].to(
            'hartree'
        ).magnitude == approx(5.64e-7)
        assert output.electronic_band_gaps[0].value.to('hartree').magnitude == approx(
            0.097, rel=1e-3
        )
        assert len(output.electronic_eigenvalues) == 2
        assert output.electronic_eigenvalues[0].value.to('hartree').magnitude[
            0, 3
        ] == approx(-0.57762)
        assert output.electronic_eigenvalues[1].value.to('hartree').magnitude[
            1, 16
        ] == approx(-0.05324)
        assert output.electronic_eigenvalues[1].occupation[0, 14] == approx(0.0)


class TestBandArchive(AMSParserIntegrationSuite):
    archive_fixture = 'band_archive'
    workflow_name = 'GeometryOptimization'
    required_simulation_sections = ('model_system', 'outputs')

    @pytest.mark.integration
    def test_archive_contract(self, archive):
        simulation = archive.data
        assert len(simulation.model_system) == 1
        assert [
            state.chemical_symbol
            for state in simulation.model_system[0].particle_states
        ] == ['Cs', 'Cl']
        assert simulation.model_system[0].positions.shape == (2, 3)
        output = simulation.outputs[0]
        assert output.total_energies[0].value.to('hartree').magnitude == approx(
            -0.23505514075497302
        )
        assert output.total_forces[0].value.shape == (2, 3)
        assert len(output.electronic_dos) == 1
        assert output.electronic_band_gaps[0].value.to('hartree').magnitude == approx(
            0.19174663326676872
        )
        assert output.electronic_eigenvalues[0].occupation[0, 3] == approx(
            1.9999999999999984
        )


class TestDOSArchive(AMSParserIntegrationSuite):
    archive_fixture = 'dos_archive'
    workflow_name = 'SinglePoint'

    @pytest.mark.integration
    def test_archive_contract(self, archive):
        simulation = archive.data
        assert [
            state.chemical_symbol
            for state in simulation.model_system[0].particle_states
        ] == ['Ni', 'O']
        assert len(simulation.outputs) == 1
        output = simulation.outputs[0]
        assert len(output.electronic_dos) == 2
        dos_up, dos_down = output.electronic_dos
        assert dos_up.spin_channel == 0
        assert dos_down.spin_channel == 1
        assert dos_up.energies.points.to('hartree').magnitude[78] == approx(-0.01989053)
        assert dos_down.value.shape == (158,)
        assert dos_down.value.to('1 / hartree').magnitude[19] == approx(
            0.0020337946278972557
        )
        assert len(output.scf_steps.delta_energies_total) == 21
        assert output.total_energies[0].value.to('hartree').magnitude == approx(
            -0.4088499099999999
        )


class TestRestrictedDOSArchive(AMSParserIntegrationSuite):
    archive_fixture = 'restricted_dos_archive'
    workflow_name = 'SinglePoint'

    @pytest.mark.integration
    def test_archive_contract(self, archive):
        simulation = archive.data
        assert [
            state.chemical_symbol
            for state in simulation.model_system[0].particle_states
        ] == ['Ni', 'O']
        assert len(simulation.outputs) == 1
        output = simulation.outputs[0]
        assert len(output.electronic_dos) == 1
        dos = output.electronic_dos[0]
        assert dos.spin_channel == 0
        assert dos.energies.points.shape == (154,)
        assert dos.energies.points.to('hartree').magnitude[78] == approx(
            -0.004229074982928244
        )
        assert len(output.scf_steps.delta_energies_total) == 17
        assert output.total_energies[0].value.to('hartree').magnitude == approx(
            -0.40559589
        )


class TestGeometryOptimizationArchive(AMSParserIntegrationSuite):
    archive_fixture = 'geometry_optimization_archive'
    workflow_name = 'GeometryOptimization'

    @pytest.mark.integration
    def test_archive_contract(self, archive):
        simulation = archive.data
        assert len(simulation.model_system) == 21
        assert len(simulation.outputs) == 22
        assert all(output.model_system_ref is not None for output in simulation.outputs)
        assert all(output.electronic_band_gaps for output in simulation.outputs)
        assert simulation.outputs[6].scf_steps.delta_energies_total.shape == (11,)
        assert simulation.outputs[4].electronic_eigenvalues[0].value.to(
            'hartree'
        ).magnitude[0, 15] == approx(-0.06523)
        assert simulation.outputs[17].total_energies[0].value.to(
            'hartree'
        ).magnitude == approx(-2.53916048)
        position = simulation.model_system[3].positions.to('bohr').magnitude
        assert position[10, 1] == approx(10.99639586)


class TestLargeGeometryOptimizationArchive(AMSParserIntegrationSuite):
    archive_fixture = 'large_geometry_optimization_archive'
    workflow_name = 'GeometryOptimization'
    required_simulation_sections = ('model_system', 'outputs')

    @pytest.mark.integration
    def test_archive_contract(self, archive):
        simulation = archive.data
        assert len(simulation.model_system) == 14
        assert len(simulation.outputs) == 14
        output = simulation.outputs[-1]
        assert output.total_energies[0].value.to('hartree').magnitude == approx(
            -974.97593754
        )
        assert output.electronic_eigenvalues[0].value.to('hartree').magnitude[
            0, 11
        ] == approx(-0.76666)
        assert output.electronic_eigenvalues[0].occupation[0, 775] == approx(2.0)
        assert output.total_forces[0].value.to('hartree / bohr').magnitude[
            163, 0
        ] == approx(-1.7468629999999995e-4)
        position = simulation.model_system[4].positions.to('bohr').magnitude
        assert position[10, 2] == approx(46.48260740098912)


class TestADFArchive(AMSParserIntegrationSuite):
    archive_fixture = 'adf_archive'
    workflow_name = 'SinglePoint'
    required_simulation_sections = ('model_system',)

    @pytest.mark.integration
    def test_representative_system_is_complete(self, archive):
        assert archive.data.model_system

    @pytest.mark.integration
    def test_archive_contract(self, archive):
        assert archive.data.program.name == 'AMS'
        assert archive.data.program.version == '98925 2021-11-25'
        assert len(archive.data.model_system) == 1
        assert [
            state.chemical_symbol
            for state in archive.data.model_system[0].particle_states
        ] == ['O', 'O']
        position = archive.data.model_system[0].positions.to('bohr').magnitude
        assert position[1, 0] == approx(-3.1867335262123166)
        assert len(archive.data.model_method) == 1
        assert archive.data.model_method[0].xc.functional_key == 'VWN'
        assert archive.data.model_system[0].periodic_boundary_conditions == [
            False,
            False,
            False,
        ]
        assert len(archive.data.outputs) == 1
        assert archive.data.outputs[0].total_energies[0].value.to(
            'hartree'
        ).magnitude == approx(-0.25256548879865254)
