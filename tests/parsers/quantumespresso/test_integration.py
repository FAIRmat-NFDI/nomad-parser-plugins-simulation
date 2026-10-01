import numpy as np
import pytest

from tests.parsers.common import (
    SimulationParserTestSuite,
    WorkflowTestSuite,
    approx,
    assert_approx,
    assert_identity_populated_once,
)


class QuantumEspressoIntegrationSuite(SimulationParserTestSuite, WorkflowTestSuite):
    expected_program_name = 'Quantum Espresso'
    identity_policy = 'once'
    require_lattice_vectors = True
    require_periodic_boundary_conditions = True

    @pytest.mark.integration
    def test_identity_populated_once(self, archive):
        systems = archive.data.model_system
        if not any(system.particle_states for system in systems):
            pytest.skip('fixture does not contain particle-state identities')
        if self.identity_policy == 'every':
            assert all(system.particle_states for system in systems)
        else:
            assert_identity_populated_once(archive)


class TestPWSCFTextArchive(QuantumEspressoIntegrationSuite):
    archive_fixture = 'pwscf_archive'
    workflow_name = 'GeometryOptimization'

    @pytest.mark.integration
    def test_archive_contract(self, archive):
        simulation = archive.data

        assert simulation.program.version == '7.3'
        assert len(simulation.model_system) == 6
        assert len(simulation.outputs) == 6
        assert [
            state.chemical_symbol
            for state in simulation.model_system[0].particle_states
        ] == [
            'Ti',
            'Ti',
            'O',
            'O',
            'O',
            'O',
        ]
        assert simulation.model_method[0].xc.functional_key == 'PBE'
        assert_approx(
            simulation.model_system[0].lattice_vectors.to('angstrom').magnitude,
            np.diag([4.59373445, 4.59373445, 2.95812152]),
        )

        outputs = simulation.outputs
        assert len(outputs[0].scf_steps.energies_total) == 12
        assert len(outputs[1].scf_steps.energies_total) == 6
        assert len(outputs[5].scf_steps.energies_total) == 14
        assert len(outputs[2].scf_steps.delta_energies_total) == 6
        assert len(outputs[3].scf_steps.delta_energies_total) == 5
        assert len(outputs[4].scf_steps.delta_energies_total) == 5
        assert len(outputs[0].total_forces) == 1
        assert len(outputs[5].total_forces) == 1
        assert outputs[0].electronic_eigenvalues[0].occupation is not None
        assert outputs[0].electronic_band_structures[0].highest_occupied is not None

    @pytest.mark.integration
    def test_workflow_convergence_contract(self, archive):
        method = archive.workflow2.method

        assert method.convergence_targets in [None, []]
        targets = method.single_point_convergence_targets
        assert targets is not None
        assert len(targets) == 1
        assert targets[0].m_def.name == 'EnergyConvergenceTarget'
        assert targets[0].threshold_type == 'absolute'
        assert targets[0].threshold.to('rydberg').magnitude == approx(1.0e-8)


class TestPWSCFXMLArchive(QuantumEspressoIntegrationSuite):
    archive_fixture = 'pwscf_xml_archive'
    workflow_name = 'GeometryOptimization'
    required_simulation_sections = ('model_method', 'outputs')

    @pytest.mark.integration
    def test_workflow_and_scf_contract(self, archive):
        targets = archive.workflow2.method.convergence_targets
        assert len(targets) == 2
        targets_by_name = {target.m_def.name: target for target in targets}
        assert targets_by_name['ForceConvergenceTarget'].threshold_type == 'maximum'
        assert targets_by_name['ForceConvergenceTarget'].threshold.to(
            'rydberg / bohr'
        ).magnitude == approx(2.5e-4)
        assert targets_by_name['EnergyConvergenceTarget'].threshold_type == 'absolute'
        assert targets_by_name['EnergyConvergenceTarget'].threshold.to(
            'rydberg'
        ).magnitude == approx(5.0e-5)

        assert len(archive.data.outputs) == 2
        for output in archive.data.outputs:
            assert output.scf_steps is not None
            assert len(output.scf_steps.energies_total) == 5
            assert len(output.scf_steps.delta_energies_total) == 5


class TestDOSArchive(QuantumEspressoIntegrationSuite):
    archive_fixture = 'dos_archive'
    workflow_name = 'SinglePoint'

    @pytest.mark.integration
    def test_dos_contract(self, archive):
        simulation = archive.data
        assert simulation.outputs
        assert all(output.model_system_ref is not None for output in simulation.outputs)

        output = archive.data.outputs[-1]
        assert output.electronic_dos
        dos = output.electronic_dos[-1]
        assert dos.value is not None
        assert dos.energies is not None
        assert dos.energies.points is not None
        assert dos.energies_origin is not None


class TestEPWArchive(QuantumEspressoIntegrationSuite):
    archive_fixture = 'epw_archive'
    required_simulation_sections = ()
    workflow_name = 'SinglePoint'


class TestPhononArchive(QuantumEspressoIntegrationSuite):
    archive_fixture = 'phonon_archive'
    identity_policy = 'every'
    required_simulation_sections = ('model_system',)
    workflow_name = 'SinglePoint'


class TestXSpectraArchive(QuantumEspressoIntegrationSuite):
    archive_fixture = 'xspectra_archive'
    required_simulation_sections = ()
    workflow_name = 'SinglePoint'


class TestGIPAWNMRTextArchive(QuantumEspressoIntegrationSuite):
    archive_fixture = 'gipaw_nmr_text_archive'
    required_simulation_sections = ()
    workflow_name = 'SinglePoint'


class TestGIPAWNMRXMLArchive(QuantumEspressoIntegrationSuite):
    archive_fixture = 'gipaw_nmr_xml_archive'
    required_simulation_sections = ()
    workflow_name = 'SinglePoint'


class TestGIPAWEFGTextArchive(QuantumEspressoIntegrationSuite):
    archive_fixture = 'gipaw_efg_text_archive'
    required_simulation_sections = ()
    workflow_name = 'SinglePoint'


class TestGIPAWEFGXMLArchive(QuantumEspressoIntegrationSuite):
    archive_fixture = 'gipaw_efg_xml_archive'
    required_simulation_sections = ()
    workflow_name = 'SinglePoint'


class TestGIPAWHyperfineTextArchive(QuantumEspressoIntegrationSuite):
    archive_fixture = 'gipaw_hyperfine_text_archive'
    required_simulation_sections = ()
    workflow_name = 'SinglePoint'


class TestGIPAWHyperfineXMLArchive(QuantumEspressoIntegrationSuite):
    archive_fixture = 'gipaw_hyperfine_xml_archive'
    required_simulation_sections = ()
    workflow_name = 'SinglePoint'


class TestGIPAWDeltaGTextArchive(QuantumEspressoIntegrationSuite):
    archive_fixture = 'gipaw_delta_g_text_archive'
    required_simulation_sections = ()
    workflow_name = 'SinglePoint'


class TestGIPAWDeltaGXMLArchive(QuantumEspressoIntegrationSuite):
    archive_fixture = 'gipaw_delta_g_xml_archive'
    required_simulation_sections = ()
    workflow_name = 'SinglePoint'
