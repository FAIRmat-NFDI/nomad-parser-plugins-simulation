import numpy as np
import pytest
from nomad.datamodel import EntryArchive
from nomad.utils import get_logger

from tests.parsers.common import (
    SimulationParserTestSuite,
    WorkflowTestSuite,
)

LOGGER = get_logger(__name__)


def approx(value, abs=0, rel=1e-6):
    return pytest.approx(value, abs=abs, rel=rel)


class ExcitingParserIntegrationSuite(SimulationParserTestSuite, WorkflowTestSuite):
    """Exciting-specific configuration of the shared integration suite."""

    expected_program_name = 'exciting'
    require_lattice_vectors = True
    require_periodic_boundary_conditions = True


class TestCMinimalArchive(ExcitingParserIntegrationSuite):
    archive_fixture = 'c_minimal_archive'
    workflow_name = 'SinglePoint'

    @pytest.mark.integration
    def test_archive_contract(self, archive):
        simulation = archive.data

        assert len(simulation.model_system) == 1
        assert [
            state.chemical_symbol
            for state in simulation.model_system[0].particle_states
        ] == ['C', 'C']
        assert simulation.model_system[0].positions.shape == (2, 3)
        assert simulation.model_system[0].lattice_vectors.shape == (3, 3)
        assert (
            list(simulation.model_system[0].periodic_boundary_conditions) == [True] * 3
        )

        method = archive.workflow2.method
        assert method is not None
        targets = method.convergence_targets
        assert [target.m_def.name for target in targets] == [
            'PotentialConvergenceTarget',
            'EnergyConvergenceTarget',
            'DensityConvergenceTarget',
        ]
        assert targets[0].threshold_type == 'rms'
        assert targets[1].threshold_type == 'absolute'
        assert targets[2].threshold_type == 'absolute'
        assert targets[0].threshold.to('hartree').magnitude == approx(
            0.100000e-05, rel=1e-12
        )

        output = simulation.outputs[0]
        output.scf_steps.delta_energies_total[9].to('hartree').magnitude == 0.991184e-07
        output.scf_steps.delta_potential_rms[9].to('hartree').magnitude == 0.324873e-08
        output.scf_steps.delta_charge_abs[9].to('coulomb').magnitude == 0.235501e-08
        assert output.total_energies[0].value.to('hartree').magnitude == approx(
            -75.88903685
        )
        assert simulation.outputs[0].electronic_band_gaps[0].value is not None

    @pytest.mark.integration
    def test_electronic_outputs_mapping(self, archive):
        output = archive.data.outputs[0]

        if output.electronic_dos:
            dos = output.electronic_dos[0]
            assert dos.value is not None
            assert dos.energies is not None
            assert dos.energies.points is not None

        if output.electronic_band_structures:
            band_structure = output.electronic_band_structures[0]
            assert band_structure.value is not None

        if output.electronic_band_gaps:
            assert output.electronic_band_gaps[0].value is not None

    @pytest.mark.integration
    def test_archive_serialization_round_trip(self, archive):
        super().test_archive_serialization_round_trip(archive)
        restored = EntryArchive.m_from_dict(archive.m_to_dict())

        assert restored.workflow2.m_def.name == 'SinglePoint'
        np.testing.assert_allclose(
            restored.data.outputs[0].scf_steps.energies_total.to('hartree').magnitude,
            archive.data.outputs[0].scf_steps.energies_total.to('hartree').magnitude,
        )


class TestCGroundStateArchive(ExcitingParserIntegrationSuite):
    archive_fixture = 'c_gs_archive'
    workflow_name = 'GeometryOptimization'

    @pytest.mark.integration
    def test_archive_contract(self, archive):
        simulation = archive.data

        assert simulation.program.name == 'exciting'
        model_system = simulation.model_system
        assert len(model_system) == 2
        assert [state.chemical_symbol for state in model_system[0].particle_states] == [
            'C',
            'C',
        ]
        assert model_system[0].lattice_vectors[0][0].to('bohr').magnitude == approx(
            3.2559441762
        )
        assert model_system[0].positions[1][0].to('bohr').magnitude == approx(
            1.627972089188372
        )

        outputs = simulation.outputs
        assert len(outputs) == 2
        assert outputs[0].total_energies[0].value.to('hartree').magnitude == approx(
            -75.89058007
        )
        assert np.mean(
            outputs[0].total_forces[0].value.to('hartree / bohr').magnitude
        ) == approx(0.0)

        assert outputs[0].scf_steps is not None
        assert len(outputs[0].scf_steps.energies_total) == 12
        assert outputs[0].scf_steps.energies_total[-1].to(
            'hartree'
        ).magnitude == approx(-75.89058007)
        assert outputs[0].total_energies[0].value.to('hartree').magnitude == approx(
            -75.89058007
        )


class TestGaOSodiumGeometryOptimization(ExcitingParserIntegrationSuite):
    archive_fixture = 'ga_o_sodium_archive'
    workflow_name = 'GeometryOptimization'

    @pytest.mark.integration
    def test_archive_contract(self, archive):
        simulation = archive.data
        workflow = archive.workflow2
        assert len(simulation.model_system) == 4
        assert len(simulation.outputs) == 4

        targets = workflow.method.convergence_targets
        assert [target.m_def.name for target in targets] == ['ForceConvergenceTarget']
        assert targets[0].threshold_type == 'maximum'
        assert targets[0].threshold.to('hartree / bohr').magnitude == approx(0.05)

        output = simulation.outputs[0]
        assert output.scf_steps.delta_energies_total[23].to(
            'hartree'
        ).magnitude == approx(0.337153e-05)
        assert output.scf_steps.delta_potential_rms[23].to(
            'hartree'
        ).magnitude == approx(0.489106e-07)
        assert output.scf_steps.delta_charge_abs[23].to('coulomb').magnitude == approx(
            4.16073e-08
        )
        assert output.scf_steps.delta_force_abs[23].to(
            'hartree / bohr'
        ).magnitude == approx(0.127181e-05)


class TestGaOStructureOptimization(ExcitingParserIntegrationSuite):
    archive_fixture = 'ga_o_strucopt_archive'
    workflow_name = 'GeometryOptimization'

    @pytest.mark.integration
    def test_archive_contract(self, archive):
        simulation = archive.data
        workflow = archive.workflow2

        assert len(simulation.model_system) == 15
        assert len(simulation.outputs) == 15
        assert simulation.model_system[0].positions.shape == (10, 3)
        assert simulation.model_system[0].lattice_vectors.shape == (3, 3)
        assert simulation.model_system[3].positions[1][1].to(
            'bohr'
        ).magnitude == approx(5.814610146853088)
        assert simulation.model_system[10].positions[-1][0].to(
            'bohr'
        ).magnitude == approx(0.6938259404131845)
        assert simulation.model_system[1].lattice_vectors[2][1].magnitude == approx(
            simulation.model_system[13].lattice_vectors[2][1].magnitude
        )

        assert simulation.outputs[0].scf_steps is not None
        assert len(simulation.outputs[0].scf_steps.energies_total) == 19
        assert simulation.outputs[-1].total_energies[0].value.to(
            'hartree'
        ).magnitude == approx(-8221.02233684)
        assert workflow.method.convergence_targets
        assert workflow.method.convergence_targets[0].threshold is not None


class TestCeODOS(ExcitingParserIntegrationSuite):
    archive_fixture = 'ce_o_dos_archive'
    workflow_name = 'SinglePoint'

    @pytest.mark.integration
    def test_archive_contract(self, archive):
        simulation = archive.data

        assert len(simulation.model_system) == 1
        assert len(simulation.outputs) == 1

        output = simulation.outputs[0]
        assert len(output.electronic_dos) == 2

        dos_up, dos_down = output.electronic_dos
        assert dos_up.value.shape == (500,)
        assert dos_down.value.shape == (500,)
        assert dos_up.energies.points.shape == (500,)
        assert dos_down.energies.points.shape == (500,)
        assert dos_up.value[126].to('1 / hartree').magnitude == approx(20.83182629)
        assert dos_down.value[136].to('1 / hartree').magnitude == approx(2.109103733)
        assert dos_up.value[220].to('1 / hartree').magnitude == approx(62.06860954)
        assert dos_down.value[78].to('1 / hartree').magnitude == approx(47.70198869)


class TestPbIHybridsArchive(ExcitingParserIntegrationSuite):
    archive_fixture = 'pb_i_hybrids_archive'
    workflow_name = 'SinglePoint'

    @pytest.mark.integration
    def test_method_mapping(self, archive):
        method = archive.data.model_method[0]
        assert method.xc is not None
        assert any(
            component.canonical_label == 'HYB_GGA_XC_HSE03'
            for component in method.xc.components
        )
        assert archive.data.outputs[0].total_energies[0].value.to(
            'hartree'
        ).magnitude == approx(-35198.96925817)
