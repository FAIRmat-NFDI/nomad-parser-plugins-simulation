import MDAnalysis
import numpy as np
import pytest
from nomad.units import ureg
from packaging.version import Version

from nomad_simulation_parsers.parsers.gromacs.parser import ENERGY_UNIT
from tests.parsers.common import (
    SimulationParserTestSuite,
    WorkflowTestSuite,
    assert_approx,
    assert_identity_populated_once,
)


class GromacsParserIntegrationSuite(SimulationParserTestSuite, WorkflowTestSuite):
    expected_program_name = 'GROMACS'
    required_simulation_sections = ('model_method', 'model_system', 'outputs')


@pytest.mark.integration
class TestGromacsWaterArchive(GromacsParserIntegrationSuite):
    archive_fixture = 'water_archive'
    workflow_name = 'MolecularDynamics'

    def test_archive_contains_trajectory_topology_and_bonds(self, archive):
        simulation = archive.data

        assert simulation.program.name == 'GROMACS'
        assert simulation.model_system
        assert simulation.outputs
        assert_identity_populated_once(archive)

        system = simulation.model_system[0]
        assert system.positions is not None
        assert system.bond_list is not None
        assert system.bond_list.shape[1] == 2
        assert system.bond_list.shape[0] == 432
        assert system.sub_systems

    def test_archive_has_valid_molecular_hierarchy(self, archive):
        system = archive.data.model_system[0]
        molecule_group = system.sub_systems[0]

        assert molecule_group.branch_label == 'molecule_group'
        assert molecule_group.sub_systems
        molecule = molecule_group.sub_systems[0]
        assert molecule.branch_label == 'molecule'
        assert molecule.composition_formula
        assert molecule.particle_indices is not None

        n_particles = system.n_particles

        def check_indices(subsystem, parent_indices=None):
            indices = subsystem.particle_indices
            assert indices is not None
            assert len(indices) > 0
            assert np.all(indices >= 0)
            assert np.all(indices < n_particles)
            if parent_indices is not None:
                assert np.all(np.isin(indices, parent_indices))
            for child in subsystem.sub_systems or []:
                check_indices(child, indices)

        for group in system.sub_systems:
            check_indices(group)


@pytest.mark.integration
class TestGromacsProteinSmallArchive(GromacsParserIntegrationSuite):
    archive_fixture = 'protein_small_archive'
    workflow_name = 'GeometryOptimization'

    def test_protein_small_archive_contains_polymer_hierarchy(self, archive):
        system = archive.data.model_system[0]

        assert system.sub_systems
        multi_residue_group = next(
            (
                group
                for group in system.sub_systems
                if group.sub_systems and group.sub_systems[0].sub_systems
            ),
            None,
        )
        assert multi_residue_group is not None
        assert multi_residue_group.branch_label == 'molecule_group'

        molecule = multi_residue_group.sub_systems[0]
        assert molecule.branch_label == 'molecule'
        assert molecule.sub_systems

        monomer_group = molecule.sub_systems[0]
        assert monomer_group.branch_label == 'monomer_group'
        assert 'group_' in monomer_group.name
        assert monomer_group.sub_systems

        monomer = monomer_group.sub_systems[0]
        assert monomer.branch_label == 'monomer'
        assert monomer.composition_formula
        assert monomer.particle_indices is not None
        assert len(monomer.particle_indices) > 0

        for group in system.sub_systems:
            assert group.branch_label == 'molecule_group'
            for molecule in group.sub_systems or []:
                assert molecule.branch_label == 'molecule'
                for monomer_group in molecule.sub_systems or []:
                    assert monomer_group.branch_label == 'monomer_group'
                    for monomer in monomer_group.sub_systems or []:
                        assert monomer.branch_label == 'monomer'


class GromacsIntegratorArchiveSuite(GromacsParserIntegrationSuite):
    required_simulation_sections = ('model_method',)
    workflow_name = 'MolecularDynamics'


@pytest.mark.integration
class TestGromacsIntegratorSdArchive(GromacsIntegratorArchiveSuite):
    archive_fixture = 'integrator_sd_archive'

    def test_integrator_settings(self, archive):
        method = archive.workflow2.method
        assert method.thermodynamic_ensemble == 'NVE'
        assert method.integrator_type == 'langevin_goga'
        assert method.thermostat_parameters


@pytest.mark.integration
class TestGromacsIntegratorMdvvArchive(GromacsIntegratorArchiveSuite):
    archive_fixture = 'integrator_mdvv_archive'

    def test_integrator_settings(self, archive):
        method = archive.workflow2.method
        assert method.thermodynamic_ensemble == 'NVE'
        assert method.integrator_type == 'velocity_verlet'


@pytest.mark.integration
class TestGromacsIntegratorBdArchive(GromacsIntegratorArchiveSuite):
    archive_fixture = 'integrator_bd_archive'

    def test_integrator_settings(self, archive):
        method = archive.workflow2.method
        assert method.thermodynamic_ensemble == 'NVE'
        assert method.integrator_type == 'brownian'


@pytest.mark.integration
class TestGromacsIntegratorVRescaleArchive(GromacsIntegratorArchiveSuite):
    archive_fixture = 'integrator_vrescale_archive'

    def test_integrator_settings(self, archive):
        method = archive.workflow2.method
        assert method.thermodynamic_ensemble == 'NVT'
        assert method.integrator_type == 'leap_frog'
        assert method.thermostat_parameters[0].thermostat_type == 'velocity_rescaling'


@pytest.mark.integration
class TestGromacsIntegratorNosehooverArchive(GromacsIntegratorArchiveSuite):
    archive_fixture = 'integrator_nosehoover_archive'

    def test_integrator_settings(self, archive):
        method = archive.workflow2.method
        assert method.thermodynamic_ensemble == 'NPT'
        assert method.integrator_type == 'leap_frog'
        assert method.thermostat_parameters[0].thermostat_type == 'nose_hoover'
        assert method.barostat_parameters[0].barostat_type == 'parrinello_rahman'


@pytest.mark.integration
class TestGromacsPolymerMeltArchive(GromacsParserIntegrationSuite):
    archive_fixture = 'polymer_melt_archive'
    required_simulation_sections = ('model_method', 'outputs')
    workflow_name = 'GeometryOptimization'

    def test_minimization_workflow_and_outputs(self, archive):
        workflow = archive.workflow2

        assert workflow.method.optimization_method == 'steepest_descent'
        assert workflow.method.n_steps_maximum == 5000
        assert len(workflow.results.energies) == 11
        assert_approx(
            workflow.results.energies[2].to(ENERGY_UNIT).magnitude,
            49650.90234375001,
        )
        assert_approx(
            workflow.results.final_force_maximum.to(
                ENERGY_UNIT / ureg.nanometer
            ).magnitude,
            676.0214199999999,
        )
        assert len(archive.data.outputs) == 11

    @pytest.mark.skipif(
        Version(MDAnalysis.__version__) > Version('2.9'),
        reason='Incompatible polymer_melt TPR file for MDAnalysis',
    )
    def test_minimization_sub_systems(self, archive):
        system = archive.data.model_system[0]

        assert len(system.sub_systems) == 1
        molecule_group = system.sub_systems[0]
        assert molecule_group.branch_label == 'molecule_group'
        assert molecule_group.name == 'group_S1P1'
        assert len(molecule_group.sub_systems) == 100
        assert molecule_group.sub_systems[52].branch_label == 'molecule'

        molecule = molecule_group.sub_systems[52]
        assert molecule.name == 'S1P1'
        assert molecule.sub_systems
        monomer_group = molecule.sub_systems[0]
        assert monomer_group.branch_label == 'monomer_group'
        assert monomer_group.name == 'group_ETHOX'
        assert monomer_group.sub_systems
        assert monomer_group.sub_systems[7].branch_label == 'monomer'


@pytest.mark.integration
class TestGromacsFeTestArchive(GromacsParserIntegrationSuite):
    archive_fixture = 'fe_test_archive'
    workflow_name = 'MolecularDynamics'

    def test_md_verbose_archive_content(self, archive):
        simulation = archive.data
        method = archive.workflow2.method
        outputs = simulation.outputs

        assert simulation.program.name == 'GROMACS'
        assert method.thermodynamic_ensemble == 'NPT'
        assert method.integrator_type == 'leap_frog'
        assert_approx(method.integration_timestep, 5e-16, rtol=1e-12)
        assert method.n_steps == 20
        assert method.coordinate_save_frequency == 20
        assert method.thermodynamics_save_frequency == 5
        assert method.thermostat_parameters[0].thermostat_type == 'berendsen'
        assert method.barostat_parameters[0].barostat_type == 'berendsen'
        assert method.barostat_parameters[0].coupling_type == 'isotropic'

        assert len(outputs) == 5
        assert_approx(outputs[3].temperatures[0].value.magnitude, 291.80401611328125)
        assert_approx(
            outputs[2].total_energies[0].value.to(ENERGY_UNIT).magnitude,
            -19663.486328125,
        )
        assert_approx(
            outputs[0].total_energies[0].value.to(ENERGY_UNIT).magnitude,
            -19699.189453125,
        )
        contributions = outputs[0].total_energies[0].contributions
        assert len(contributions) == 3
        assert_approx(
            contributions[0].value.to(ENERGY_UNIT).magnitude,
            -23393.28125,
        )
        assert_approx(
            contributions[1].value.to(ENERGY_UNIT).magnitude,
            3694.0922851562505,
        )
        assert_approx(
            contributions[2].value.to(ENERGY_UNIT).magnitude,
            0.9067516922950745,
        )
        assert outputs[0].total_forces[0].value.shape == (1516, 3)
        assert_approx(
            outputs[0]
            .total_forces[0]
            .value[5][2]
            .to(ENERGY_UNIT / ureg.nanometer)
            .magnitude,
            -477.73454314849005,
        )

        assert len(simulation.model_system) == 2
        first_system, second_system = simulation.model_system
        assert first_system.positions.shape == (1516, 3)
        assert_approx(
            second_system.positions[800][1].to('nanometer').magnitude,
            2.4740036685955142,
        )
        assert_approx(
            first_system.velocities[500][0].to('nanometer / picosecond').magnitude,
            0.8694772949218749,
        )
        assert_approx(
            first_system.lattice_vectors[2][2].to('nanometer').magnitude,
            2.469330073751052,
        )
        assert first_system.bond_list[200, 0] == 289
        assert simulation.model_method[0].contributions


@pytest.mark.integration
class TestGromacsFreeEnergyArchive(GromacsParserIntegrationSuite):
    archive_fixture = 'fep_archive'
    workflow_name = 'MolecularDynamics'

    @pytest.mark.integration
    def test_representative_system_is_complete(self, archive):
        representative = archive.data.model_system[0]
        assert representative.positions is not None
        assert representative.particle_states

    def test_fep_xvg_fields_are_populated(self, archive):
        workflow = archive.workflow2
        assert workflow is not None
        method = workflow.method
        assert method is not None
        parameters = method.free_energy_calculation_parameters
        assert parameters

        fep = parameters[0]
        assert fep.n_frames == 5001
        assert fep.n_states == 11
        assert fep.times is not None
        assert len(fep.times) == 5001
        assert fep.energy_derivative is not None
        assert fep.energy_differences is not None
        assert fep.energy_differences.shape == (5001, 11)
        assert fep.pv_energy is not None
