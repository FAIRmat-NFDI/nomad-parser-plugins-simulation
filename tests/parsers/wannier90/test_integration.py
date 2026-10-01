import pytest

from tests.parsers.common import (
    SimulationParserTestSuite,
    WorkflowTestSuite,
    approx,
    assert_approx,
)


class Wannier90IntegrationSuite(SimulationParserTestSuite, WorkflowTestSuite):
    expected_program_name = 'Wannier90'
    workflow_name = 'SinglePoint'
    required_simulation_sections = ('model_method', 'outputs')


class TestLcoMlwfArchive(Wannier90IntegrationSuite):
    archive_fixture = 'lco_mlwf_archive'

    @pytest.mark.integration
    def test_archive_contract(self, archive):
        simulation = archive.data
        assert len(simulation.model_system) == 1
        system = simulation.model_system[0]
        assert_approx(
            system.positions.to('angstrom').magnitude,
            [
                [0.0, 0.0, 4.77028],
                [1.90914, 1.90914, 1.83281],
                [0.0, 0.0, 0.0],
                [0.0, 1.90914, 0.0],
                [1.90914, 0.0, 0.0],
                [0.0, 0.0, 2.45222],
                [1.90914, 1.90914, 4.15087],
            ],
        )
        assert_approx(
            system.lattice_vectors.to('angstrom').magnitude,
            [
                [-1.909145, 1.909145, 6.603098],
                [1.909145, -1.909145, 6.603098],
                [1.909145, 1.909145, -6.603098],
            ],
        )
        assert len(simulation.model_method) == 1
        assert len(simulation.outputs) == 1
        method = simulation.model_method[0]
        assert method.is_maximally_localized is True
        assert_approx(method.energy_window_outer.to('eV').magnitude, [10.0, 16.0])

    @pytest.mark.integration
    def test_electronic_outputs_and_workflow(self, archive):
        output = archive.data.outputs[0]
        assert output.electronic_band_structures
        bands = output.electronic_band_structures[0]
        assert bands.value is not None
        # First and last values from lco_band.dat (energies as written by
        # Wannier90; the current mapping stores them with a joule unit).
        assert_approx(bands.value.magnitude[:3, 0], [10.837560, 10.838521, 10.841402])
        assert_approx(bands.value.magnitude[-1, 0], 10.837618)
        assert bands.highest_occupied.to('eV').magnitude == approx(11.375)
        assert output.crystal_field_splittings
        assert output.crystal_field_splittings[0].n_orbitals == 1
        assert archive.workflow2.m_def.name == 'SinglePoint'
