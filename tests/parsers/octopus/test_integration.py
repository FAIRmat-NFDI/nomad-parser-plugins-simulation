import numpy as np
import pytest

from tests.parsers.common import SimulationParserTestSuite, approx, assert_approx


class OctopusParserIntegrationSuite(SimulationParserTestSuite):
    expected_program_name = 'Octopus'
    require_lattice_vectors = True
    require_periodic_boundary_conditions = True

    @pytest.mark.integration
    def test_model_system_contract(self, archive):
        simulation = archive.data

        assert len(simulation.model_system) >= 1
        representative = next(
            (
                system
                for system in simulation.model_system
                if getattr(system, 'is_representative', False)
            ),
            simulation.model_system[0],
        )

        assert representative.positions is not None
        assert representative.positions.shape[1] == 3
        assert representative.lattice_vectors is not None
        assert representative.lattice_vectors.shape == (3, 3)
        assert list(representative.periodic_boundary_conditions) == [True] * 3
        assert representative.particle_states
        assert all(
            state.chemical_symbol is not None
            for state in representative.particle_states
        )

    @pytest.mark.integration
    def test_outputs_contract(self, archive):
        outputs = archive.data.outputs

        assert outputs
        output = outputs[-1]
        assert output.total_energies
        assert output.total_energies[0].value is not None
        assert output.total_forces
        assert output.total_forces[0].value is not None

        if output.electronic_eigenvalues:
            eigenvalues = output.electronic_eigenvalues[0]
            assert eigenvalues.occupation is not None

        if output.electronic_band_structures:
            assert output.electronic_band_structures[0].value is not None

        if output.electronic_band_gaps:
            assert output.electronic_band_gaps[0].value is not None


class TestSiScfArchive(OctopusParserIntegrationSuite):
    archive_fixture = 'si_scf_archive'

    @pytest.mark.integration
    def test_archive_contract(self, archive):
        simulation = archive.data

        assert simulation.program.version == 'wolfi'

        method = simulation.model_method[0]
        assert method.m_def.name == 'DFT'
        assert [component.canonical_label for component in method.xc.components] == [
            'LDA_X',
            'LDA_C_PZ_MOD',
        ]
        system = simulation.model_system[0]
        assert [state.chemical_symbol for state in system.particle_states] == [
            'Si',
            'Si',
            'Si',
            'Si',
        ]
        assert_approx(
            system.lattice_vectors.to('angstrom').magnitude,
            [
                [3.83958972, 0.0, 0.0],
                [0.0, 3.83958972, 0.0],
                [0.0, 0.0, 5.43000008],
            ],
        )
        assert system.positions[1][0].to('angstrom').magnitude == approx(1.91979671)

        output = simulation.outputs[0]
        assert output.total_energies[0].value.to('eV').magnitude == approx(
            -431.67878767
        )
        assert_approx(
            output.total_forces[0].value.to('eV / angstrom').magnitude,
            0.0,
        )

        eigenvalues = output.electronic_eigenvalues[0]
        assert eigenvalues.occupation.shape == (18, 8)
        assert eigenvalues.occupation[16][1] == approx(2.0)

        band_structure = output.electronic_band_structures[0]
        assert band_structure.value.shape == (18, 8)
        assert band_structure.value[4][6].to('eV').magnitude == approx(3.28702661)


class TestFeSpinpolArchive(OctopusParserIntegrationSuite):
    archive_fixture = 'fe_spinpol_archive'

    @pytest.mark.integration
    def test_archive_contract(self, archive):
        simulation = archive.data
        assert simulation.program.version == 'mimus'
        assert [
            state.chemical_symbol
            for state in simulation.model_system[0].particle_states
        ] == ['Fe', 'Fe']
        assert_approx(
            simulation.model_system[0].lattice_vectors.to('angstrom').magnitude,
            np.eye(3) * 2.87,
        )

        output = archive.data.outputs[-1]
        assert output.total_energies
        assert output.total_forces
        assert output.total_energies[0].value.to('eV').magnitude == approx(
            -6437.29578945
        )
