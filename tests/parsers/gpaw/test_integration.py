import pytest
from nomad.datamodel import EntryArchive

from tests.parsers.common import (
    SimulationParserTestSuite,
    WorkflowTestSuite,
    approx,
    assert_approx,
)


class GPAWParserIntegrationSuite(SimulationParserTestSuite, WorkflowTestSuite):
    """GPAW-specific configuration of the shared integration suite."""

    expected_program_name = 'GPAW'
    require_lattice_vectors = True
    require_periodic_boundary_conditions = True


class TestFe2Archive(GPAWParserIntegrationSuite):
    archive_fixture = 'fe2_archive'
    workflow_name = 'SinglePoint'

    @pytest.mark.integration
    def test_archive_contract(self, archive):
        simulation = archive.data

        assert len(simulation.model_system) == 1
        system = simulation.model_system[0]
        assert [state.chemical_symbol for state in system.particle_states] == [
            'Fe',
            'Fe',
        ]
        assert system.positions.shape == (2, 3)
        assert system.lattice_vectors.shape == (3, 3)
        assert list(system.periodic_boundary_conditions) == [True] * 3

        output = simulation.outputs[0]
        assert output.total_energies[0].value.to('eV').magnitude == approx(
            -7.301259879298866
        )
        assert len(output.electronic_eigenvalues) == 2

        target = archive.workflow2.method.convergence_targets[0]
        assert target.m_def.name == 'EnergyConvergenceTarget'
        assert target.threshold_type == 'absolute'
        assert target.threshold.to('eV').magnitude == approx(4.0930753554401515e-08)

    @pytest.mark.integration
    def test_electronic_outputs_mapping(self, archive):
        output = archive.data.outputs[0]

        for eigenvalues in output.electronic_eigenvalues:
            assert eigenvalues.value is not None
            assert eigenvalues.occupation is not None
            assert eigenvalues.highest_occupied is not None

        if output.electronic_band_structures:
            assert output.electronic_band_structures[0].value is not None
        if output.electronic_band_gaps:
            assert output.electronic_band_gaps[0].value is not None

    @pytest.mark.integration
    def test_archive_serialization_round_trip(self, archive):
        super().test_archive_serialization_round_trip(archive)
        restored = EntryArchive.m_from_dict(archive.m_to_dict())

        assert_approx(
            restored.data.outputs[0].total_energies[0].value.to('eV').magnitude,
            archive.data.outputs[0].total_energies[0].value.to('eV').magnitude,
        )


class TestH2Archive(GPAWParserIntegrationSuite):
    archive_fixture = 'h2_archive'
    workflow_name = 'SinglePoint'

    @pytest.mark.integration
    def test_archive_contract(self, archive):
        simulation = archive.data
        system = simulation.model_system[0]
        output = simulation.outputs[0]
        eigenvalues = output.electronic_eigenvalues[0]

        assert [state.chemical_symbol for state in system.particle_states] == ['H', 'H']
        assert list(system.periodic_boundary_conditions) == [False] * 3
        assert_approx(
            system.lattice_vectors.to('angstrom').magnitude,
            [
                [4.0, 0.0, 0.0],
                [0.0, 4.0, 0.0],
                [0.0, 0.0, 4.73716558],
            ],
        )
        assert_approx(
            system.positions.to('angstrom').magnitude,
            [
                [2.0, 2.0, 2.73716576],
                [2.0, 2.0, 2.0],
            ],
        )
        assert output.total_energies[0].value.to('eV').magnitude == approx(
            -9.45057938402476
        )
        assert eigenvalues.value.shape == (1, 2)
        assert_approx(
            eigenvalues.value.to('eV').magnitude,
            [[-9.70134724, 3.46080199]],
            rtol=1e-7,
        )
        assert_approx(eigenvalues.occupation, [[2.0, 0.0]])
        assert output.electronic_band_gaps


class TestSiPWArchive(GPAWParserIntegrationSuite):
    archive_fixture = 'si_pw_archive'
    workflow_name = 'SinglePoint'

    @pytest.mark.integration
    def test_archive_contract(self, archive):
        simulation = archive.data
        system = simulation.model_system[0]
        output = simulation.outputs[0]

        assert simulation.program.version == '1.1.1b1'
        assert [state.chemical_symbol for state in system.particle_states] == [
            'Si',
            'Si',
        ]
        assert list(system.periodic_boundary_conditions) == [True] * 3
        assert system.lattice_vectors[0][1].to('angstrom').magnitude == approx(2.715)
        assert system.positions[1][0].to('angstrom').magnitude == approx(1.3575)
        assert output.electronic_eigenvalues[0].value.shape == (10, 8)
        assert output.electronic_eigenvalues[0].value[7][4].to(
            'eV'
        ).magnitude == approx(5.93437008109612)
        assert output.electronic_eigenvalues[0].occupation[0][0] == approx(1.0)
        assert output.electronic_eigenvalues[0].highest_occupied is not None
        assert output.total_energies[0].contributions[1].value.to(
            'eV'
        ).magnitude == approx(-13.707842612674868)


class TestSiLCAOArchive(GPAWParserIntegrationSuite):
    archive_fixture = 'si_lcao_archive'
    workflow_name = 'SinglePoint'

    @pytest.mark.integration
    def test_archive_contract(self, archive):
        simulation = archive.data
        system = simulation.model_system[0]
        output = simulation.outputs[0]

        assert simulation.program.version == '1.1.1b1'
        assert [state.chemical_symbol for state in system.particle_states] == [
            'Si',
            'Si',
        ]
        assert system.positions.shape == (2, 3)
        assert system.lattice_vectors.shape == (3, 3)
        assert output.electronic_eigenvalues[0].value.shape == (10, 8)
        assert output.total_energies


class TestHSpinPolarizedArchive(GPAWParserIntegrationSuite):
    archive_fixture = 'hspinpol_archive'
    workflow_name = 'SinglePoint'

    @pytest.mark.integration
    def test_spin_resolved_eigenvalues(self, archive):
        eigenvalues = archive.data.outputs[0].electronic_eigenvalues

        assert len(eigenvalues) == 2
        assert eigenvalues[0].spin_channel == 0
        assert eigenvalues[1].spin_channel == 1
        assert eigenvalues[0].value.shape == (1, 1)
        assert eigenvalues[1].value.shape == (1, 1)
        assert eigenvalues[0].value[0][0].to('eV').magnitude == approx(
            -7.551182149856757
        )
