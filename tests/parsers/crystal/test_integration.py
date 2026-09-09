import pytest

from tests.parsers.common import (
    SimulationParserTestSuite,
    WorkflowTestSuite,
)


def approx(value, abs=0, rel=1e-6):
    return pytest.approx(value, abs=abs, rel=rel)


class CrystalParserIntegrationSuite(SimulationParserTestSuite, WorkflowTestSuite):
    expected_program_name = 'Crystal'
    is_periodic = True

    @pytest.mark.integration
    def test_basic(self, archive):
        errors, warnings = archive.m_validate()
        assert errors == []
        assert warnings == []

        simulation = archive.data
        assert simulation.program.version is not None

        for system in simulation.model_system:
            assert system.positions is not None
            assert system.particle_states
            assert all(
                state.chemical_symbol for state in system.particle_states
            )
            if not self.is_periodic:
                continue
            assert system.lattice_vectors is not None
            assert system.lattice_vectors.shape == (3, 3)
            assert len(system.periodic_boundary_conditions) == 3
            assert all(
                isinstance(periodic, bool)
                for periodic in system.periodic_boundary_conditions
            )

class CrystalParserSimulationIntegrationSuite(SimulationParserTestSuite):
    expected_program_name = 'Crystal'


class TestSinglePointArchive(CrystalParserIntegrationSuite):
    archive_fixture = 'single_point_archive'
    workflow_name = 'SinglePoint'

    @pytest.mark.integration
    def test_archive_contract(self, archive):
        simulation = archive.data
        assert simulation.model_method
        assert simulation.model_method[0].m_def.name == 'DFT'

        assert len(simulation.model_system) == 1
        system = simulation.model_system[0]
        assert system.positions.shape == (2, 3)
        assert [state.chemical_symbol for state in system.particle_states] == ['Si', 'Si']
        assert {
            component.canonical_label
            for component in simulation.model_method[0].xc.components
        } == {'LDA_C_PZ', 'GGA_X_B88'}
        assert len(simulation.outputs) == 1
        output = simulation.outputs[0]
        assert len(output.scf_steps.energies_total) == 8
        assert output.scf_steps.energies_total[-1].to('hartree').magnitude == approx(
            -573.300583798
        )
        assert output.scf_steps.delta_energies_total[-1].to('hartree').magnitude == approx(
            5.73e-8
        )
        assert output.total_energies[0].value.to('hartree').magnitude == approx(
            -573.30058382967
        )


class TestGeometryOptimizationArchive(CrystalParserIntegrationSuite):
    archive_fixture = 'geometry_optimization_archive'
    workflow_name = 'GeometryOptimization'

    @pytest.mark.integration
    def test_geometry_optimization(self, archive):
        workflow = archive.workflow2
        assert workflow.m_def.name == 'GeometryOptimization'
        assert workflow.method.convergence_targets
        assert workflow.method.convergence_targets[0].threshold is not None

    @pytest.mark.integration
    def test_forces(self, archive):
        forces = archive.data.outputs[0].total_forces
        assert forces
        assert forces[0].value.shape == (4, 3)

    @pytest.mark.integration
    def test_archive_contract(self, archive):
        simulation = archive.data
        assert simulation.model_method
        assert simulation.model_method[0].m_def.name == 'DFT'
        assert len(simulation.model_system) > 1
        assert len(simulation.outputs) > 1
        threshold = archive.workflow2.method.convergence_targets[0].threshold
        assert threshold.to('hartree').magnitude == approx(1e-7)
        assert simulation.model_system[0].positions.shape == (4, 3)


class TestBandArchive(CrystalParserIntegrationSuite):
    archive_fixture = 'band_archive'
    workflow_name = 'SinglePoint'
    required_simulation_sections = ('model_system', 'outputs')

    @pytest.mark.integration
    def test_band_structures_and_dos(self, archive):
        output = archive.data.outputs[0]
        assert output.electronic_band_structures
        for band_structure in output.electronic_band_structures:
            assert band_structure.value is not None
            assert band_structure.value.ndim == 2
        assert output.electronic_dos
        for dos in output.electronic_dos:
            assert dos.energies is not None
            assert dos.energies.points is not None
            assert dos.value is not None
            assert dos.energies.points.shape == dos.value.shape

    @pytest.mark.integration
    def test_archive_contract(self, archive):
        simulation = archive.data
        assert len(simulation.outputs) == 1
        output = simulation.outputs[0]
        assert len(output.electronic_band_structures) == 4
        assert output.electronic_band_structures[0].value.shape == (20, 18)
        assert output.electronic_band_structures[0].value.to('hartree').magnitude[
            0, 14
        ] == approx(0.128591)
        assert output.electronic_dos[0].value.shape == (302,)
        assert output.electronic_dos[0].value.to('1 / hartree').magnitude[9] == approx(
            295.046
        )


class TestDOSArchive(CrystalParserIntegrationSuite):
    archive_fixture = 'dos_archive'
    workflow_name = 'SinglePoint'
    required_simulation_sections = ('model_system', 'outputs')

    @pytest.mark.integration
    def test_band_structures_and_dos(self, archive):
        output = archive.data.outputs[0]
        assert output.electronic_band_structures
        for band_structure in output.electronic_band_structures:
            assert band_structure.value is not None
            assert band_structure.value.ndim == 2
        assert output.electronic_dos
        for dos in output.electronic_dos:
            assert dos.energies is not None
            assert dos.energies.points is not None
            assert dos.value is not None
            assert dos.energies.points.shape == dos.value.shape

    @pytest.mark.integration
    def test_archive_contract(self, archive):
        simulation = archive.data
        assert len(simulation.outputs) == 1
        output = simulation.outputs[0]
        assert len(output.electronic_dos) == 1
        dos = output.electronic_dos[0]
        assert dos.energies.points.shape == (302,)
        assert dos.energies.points.to('hartree').magnitude[0] == approx(-3.9685800e-01)
        assert dos.energies.points.to('hartree').magnitude[78] == approx(-0.2320713)
        assert dos.value.shape == (302,)
        assert dos.value.to('1 / hartree').magnitude[9] == approx(295.046)


class TestMoleculeArchive(CrystalParserIntegrationSuite):
    archive_fixture = 'molecule_archive'
    workflow_name = 'GeometryOptimization'
    required_simulation_sections = ('outputs',)
    is_periodic = False

    @pytest.mark.integration
    def test_geometry_optimization(self, archive):
        workflow = archive.workflow2
        assert workflow.m_def.name == 'GeometryOptimization'
        assert workflow.method.convergence_targets
        assert workflow.method.convergence_targets[0].threshold is not None

    @pytest.mark.integration
    def test_archive_contract(self, archive):
        simulation = archive.data
        assert len(simulation.outputs) == 5


class TestBandNoF25Archive(CrystalParserIntegrationSuite):
    archive_fixture = 'band_no_f25_archive'
    workflow_name = 'SinglePoint'
    required_simulation_sections = ('model_system', 'outputs')


class TestBandNoF25VariantArchive(CrystalParserIntegrationSuite):
    archive_fixture = 'band_no_f25_variant_archive'
    workflow_name = 'SinglePoint'
    required_simulation_sections = ('model_system',)


class TestConstraintsArchive(CrystalParserSimulationIntegrationSuite):
    archive_fixture = 'constraints_archive'
    required_simulation_sections = ('model_system', 'outputs')


class TestDisplacementArchive(CrystalParserSimulationIntegrationSuite):
    archive_fixture = 'displacement_archive'
    required_simulation_sections = ('model_system', 'outputs')


class TestGhostsArchive(CrystalParserSimulationIntegrationSuite):
    archive_fixture = 'ghosts_archive'
    required_simulation_sections = ('model_system', 'outputs')


class TestNATArchive(CrystalParserIntegrationSuite):
    archive_fixture = 'nat_archive'
    workflow_name = 'SinglePoint'


class TestSubstitutionArchive(CrystalParserSimulationIntegrationSuite):
    archive_fixture = 'substitution_archive'


class TestNanotubeGeometryOptimizationArchive(CrystalParserIntegrationSuite):
    archive_fixture = 'nanotube_geometry_optimization_archive'
    workflow_name = 'GeometryOptimization'


class TestNanotubeSCFArchive(CrystalParserIntegrationSuite):
    archive_fixture = 'nanotube_scf_archive'
    workflow_name = 'SinglePoint'


class TestSinglePointForcesArchive(CrystalParserIntegrationSuite):
    archive_fixture = 'single_point_forces_archive'
    workflow_name = 'SinglePoint'

    @pytest.mark.integration
    def test_forces(self, archive):
        assert archive.data.outputs[0].total_forces


class TestSinglePointHFArchive(CrystalParserIntegrationSuite):
    archive_fixture = 'single_point_hf_archive'
    workflow_name = 'SinglePoint'
    required_simulation_sections = ('model_system', 'outputs')


class TestSurfaceArchive(CrystalParserSimulationIntegrationSuite):
    archive_fixture = 'surface_archive'
    required_simulation_sections = ('model_system', 'outputs')


class TestPBE0Archive(CrystalParserIntegrationSuite):
    archive_fixture = 'pbe0_archive'
    workflow_name = 'SinglePoint'
    required_simulation_sections = ('model_system',)


class TestPBEArchive(CrystalParserIntegrationSuite):
    archive_fixture = 'pbe_archive'
    workflow_name = 'SinglePoint'


class TestPW91HybridArchive(CrystalParserIntegrationSuite):
    archive_fixture = 'pw91_hybrid_archive'
    workflow_name = 'SinglePoint'
    required_simulation_sections = ('model_system', 'outputs')


class TestWC1LYPArchive(CrystalParserIntegrationSuite):
    archive_fixture = 'wc1lyp_archive'
    workflow_name = 'SinglePoint'
    required_simulation_sections = ('model_system', 'outputs')
