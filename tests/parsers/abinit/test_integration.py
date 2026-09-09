import numpy as np
import pytest
from nomad.datamodel import EntryArchive

from tests.parsers.common import SimulationParserTestSuite, WorkflowTestSuite


def approx(value, abs=0, rel=1e-6):
    return pytest.approx(value, abs=abs, rel=rel)


class AbinitParserIntegrationSuite(SimulationParserTestSuite, WorkflowTestSuite):
    expected_program_name = 'ABINIT'


class TestFeSinglePointArchive(AbinitParserIntegrationSuite):
    archive_fixture = 'fe_archive'
    workflow_name = 'SinglePoint'

    @pytest.mark.integration
    def test_archive_contract(self, archive):  # noqa: PLR0915
        errors, warnings = archive.m_validate()
        assert errors == []
        assert warnings == []

        simulation = archive.data
        assert simulation.program.name == 'ABINIT'
        assert simulation.program.version == '7.8.2'
        assert len(simulation.model_system) == 1
        assert len(simulation.model_method) > 0
        system = simulation.model_system[0]
        assert [state.chemical_symbol for state in system.particle_states] == ['Fe']
        assert system.positions.shape == (1, 3)
        assert system.lattice_vectors.shape == (3, 3)
        assert list(system.periodic_boundary_conditions) == [True, True, True]

        assert len(simulation.outputs) == 2
        first_steps = simulation.outputs[0].scf_steps
        second_steps = simulation.outputs[1].scf_steps
        assert len(first_steps.energies_total) == 16
        assert len(second_steps.energies_total) == 30
        assert first_steps.energies_total[0].to('hartree').magnitude == approx(
            -23.45196
        )
        assert first_steps.energies_total[-1].to('hartree').magnitude == approx(
            -24.6617073
        )
        assert first_steps.delta_energies_total[-1].to('hartree').magnitude == approx(
            1.243e-13
        )
        assert second_steps.delta_energies_total[-1].to('hartree').magnitude == approx(
            0.0
        )

        band_structure = simulation.outputs[0].electronic_band_structures[0]
        assert band_structure.value.shape == (3, 8)
        assert band_structure.value.to('hartree').magnitude[0, 0] == approx(-0.47539)
        dos = simulation.outputs[1].electronic_dos[0]
        assert dos.value.shape == (1601,)
        assert dos.energies.points.shape == (1601,)
        assert dos.energies.points.to('hartree').magnitude[543] == approx(-0.257)
        assert dos.value.to('1 / hartree').magnitude[151] == approx(0.001316)

        assert archive.workflow2.m_def.name == 'SinglePoint'
        assert archive.workflow2.method is None

    @pytest.mark.integration
    def test_archive_serialization_round_trip(self, archive):
        restored = EntryArchive.m_from_dict(archive.m_to_dict())
        assert restored.data.program.name == 'ABINIT'
        assert len(restored.data.model_system) == 1
        assert len(restored.data.outputs) == 2
        assert restored.data.outputs[0].total_energies[0].value.magnitude == approx(
            archive.data.outputs[0].total_energies[0].value.magnitude
        )


class TestH2GeometryOptimizationArchive(AbinitParserIntegrationSuite):
    archive_fixture = 'h2_archive'
    workflow_name = 'GeometryOptimization'

    @pytest.mark.integration
    def test_archive_contract(self, archive):
        workflow = archive.workflow2
        assert workflow.m_def.name == 'GeometryOptimization'
        assert workflow.method.optimization_method == 'bfgs'
        targets = {
            target.m_def.name: target for target in workflow.method.convergence_targets
        }
        assert set(targets) == {'EnergyConvergenceTarget', 'ForceConvergenceTarget'}
        assert targets['EnergyConvergenceTarget'].threshold.to(
            'hartree'
        ).magnitude == approx(0.0)
        assert targets['ForceConvergenceTarget'].threshold.to(
            'hartree / bohr'
        ).magnitude == approx(5.0e-4)
        assert archive.data.outputs[3].total_energies[0].value.to(
            'hartree'
        ).magnitude == approx(-1.13305855)
        assert len(archive.data.outputs[3].scf_steps.energies_total) == 5
        assert archive.data.outputs[3].scf_steps.energies_total[-1].to(
            'hartree'
        ).magnitude == approx(-1.1330585487294)


class TestSiSinglePointArchive(AbinitParserIntegrationSuite):
    archive_fixture = 'si_archive'
    workflow_name = 'SinglePoint'

    @pytest.mark.integration
    def test_archive_contract(self, archive):
        errors, warnings = archive.m_validate()
        assert errors == []
        assert warnings == []
        simulation = archive.data
        assert simulation.program.version == '7.8.2'
        assert simulation.datetime.timestamp() == approx(1467132480.0)
        assert len(simulation.model_system) == 1
        system = simulation.model_system[0]
        assert [state.chemical_symbol for state in system.particle_states] == [
            'Si',
            'Si',
        ]
        assert system.positions.shape == (2, 3)
        assert system.lattice_vectors.shape == (3, 3)
        assert list(system.periodic_boundary_conditions) == [True, True, True]
        assert system.positions[1][1].to('bohr').magnitude == approx(2.545)
        assert system.lattice_vectors[2][0].to('bohr').magnitude == approx(5.09)
        assert (
            simulation.model_method[0].xc.components[0].canonical_label
            == 'LDA_XC_TETER93'
        )

        output = simulation.outputs[0]
        assert np.max(
            output.total_forces[0].value.to('hartree / bohr').magnitude
        ) == approx(0.0)
        assert output.total_energies[0].value.to('hartree').magnitude == approx(
            -8.866223895975928
        )
        assert output.electronic_band_structures[0].value.shape == (2, 5)
        assert archive.workflow2.method is not None
        assert len(archive.workflow2.method.convergence_targets) == 1


class TestZrO2GWArchive(AbinitParserIntegrationSuite):
    archive_fixture = 'zro2_gw_archive'
    workflow_name = 'SinglePoint'

    @pytest.mark.integration
    def test_archive_contract(self, archive):
        simulation = archive.data
        assert simulation.program.name == 'ABINIT'
        assert simulation.program.version == '9.4.0'
        assert len(simulation.model_system) == 1
        assert [
            state.chemical_symbol
            for state in simulation.model_system[0].particle_states
        ] == ['Zr', 'O', 'O']
        assert len(simulation.outputs) == 1
        assert simulation.outputs[0].electronic_band_structures[0].value.shape == (
            1,
            512,
        )
