import numpy as np
import pytest

from tests.parsers.common import (
    SimulationParserTestSuite,
    WorkflowTestSuite,
    approx,
    assert_approx,
)


class H5MDParserIntegrationSuite(SimulationParserTestSuite, WorkflowTestSuite):
    expected_program_name = 'OpenMM'
    workflow_name = 'MolecularDynamics'
    require_lattice_vectors = True
    require_periodic_boundary_conditions = True

    @pytest.mark.integration
    def test_representative_system_is_complete(self, archive):
        """Validate H5MD topology identity and trajectory cell data separately.

        H5MD stores particle identity on the topology frame only, while the
        representative frame selected by normalization can be a later frame.
        Consequently, the generic single-frame completeness contract does not
        apply to this trajectory representation.
        """
        topology = archive.data.model_system[0]
        assert topology.positions is not None
        assert topology.particle_states
        assert all(
            state.chemical_symbol is not None for state in topology.particle_states
        )

        assert any(
            system.lattice_vectors is not None
            for system in archive.data.model_system
        )
        assert any(
            system.periodic_boundary_conditions is not None
            for system in archive.data.model_system
        )
    required_simulation_sections = ['model_system', 'outputs']


class TestH5MDArchive(H5MDParserIntegrationSuite):
    archive_fixture = 'h5md_archive'

    def test_h5md_archive_contract(self, archive):
        simulation = archive.data

        assert simulation.program.name == 'OpenMM'
        assert simulation.program.version == '-1.-1.-1'
        assert list(simulation.x_h5md_version) == [1, 0]
        assert simulation.x_h5md_author.name == 'Joseph F. Rudzinski'
        assert simulation.x_h5md_author.email == 'joseph.rudzinski@physik.hu-berlin.de'
        assert simulation.x_h5md_creator.name == 'h5py'
        assert simulation.x_h5md_creator.version == '3.6.0'

        systems = simulation.model_system
        assert len(systems) == 5
        assert systems[0].n_particles == 728
        assert systems[0].positions.shape == (728, 3)
        assert systems[0].velocities.shape == (728, 3)
        assert systems[2].positions[80][1].to('angstrom').magnitude == approx(
            28.748762
        )
        assert systems[2].velocities[50][2].to('angstrom/ps').magnitude == approx(400.0)
        assert systems[3].lattice_vectors[2][2].to('angstrom').magnitude == approx(
            68.22318
        )
        assert list(systems[3].periodic_boundary_conditions) == [True] * 3
        assert systems[0].bond_list[200][0] == 198
        assert systems[0].dimensionality == 3
        assert systems[0].is_molecule() is False

        assert len(systems[0].sub_systems) == 4
        atoms_group = systems[0].sub_systems[0]
        assert atoms_group.particle_states == []
        assert atoms_group.name == 'group_1ZNF'
        assert atoms_group.branch_label == 'molecule_group'
        assert atoms_group.composition_formula == '1ZNF(1)'
        assert atoms_group.particle_indices[159] == 159
        assert atoms_group.is_molecule() is True

        protein = atoms_group.sub_systems[0]
        assert len(atoms_group.sub_systems) == 1
        assert protein.name == '1ZNF'
        assert protein.branch_label == 'molecule'
        assert protein.composition_formula == (
            'ACE(1)TYR(1)LYS(3)CYS(2)GLY(1)LEU(2)GLU(2)ARG(3)SER(3)PHE(1)'
            'VAL(2)ALA(1)HIS(2)GLN(1)ASN(1)NH2(1)'
        )
        assert protein.particle_indices[400] == 400
        assert protein.is_molecule() is True

        residue_groups = protein.sub_systems
        assert len(residue_groups) == 16
        assert residue_groups[13].name == 'group_ARG'
        assert residue_groups[13].composition_formula == 'ARG(3)'
        assert residue_groups[14].branch_label == 'monomer_group'
        assert residue_groups[14].particle_indices[2] == 136
        assert residue_groups[14].is_molecule() is False

        residues = residue_groups[13].sub_systems
        assert len(residues) == 3
        assert residues[0].name == 'ARG'
        assert residues[0].branch_label == 'monomer'
        assert residues[0].particle_indices[10] == 120
        assert residues[0].is_molecule() is False
        assert residues[0].composition_formula == (
            'C(1)CA(1)CB(1)CD(1)CG(1)CZ(1)H(1)HA(1)HB2(1)HB3(1)HD2(1)'
            'HD3(1)HE(1)HG2(1)HG3(1)HH11(1)HH12(1)HH21(1)HH22(1)N(1)'
            'NE(1)NH1(1)NH2(1)O(1)'
        )

        outputs = simulation.outputs
        assert len(outputs) == 5
        assert outputs[3].step == 3
        assert outputs[2].time.to('ps').magnitude == approx(2.0)
        assert outputs[2].temperatures[0].value.to('kelvin').magnitude == approx(300.0)
        assert outputs[2].total_energies[0].value.to('kilojoule').magnitude == approx(
            6.0
        )
        total_energy = outputs[2].total_energies[0]
        assert total_energy.contributions[0].name == 'BaseEnergy'
        assert total_energy.contributions[0].contribution_type == 'custom'
        assert total_energy.contributions[0].value.to(
            'kilojoule'
        ).magnitude == approx(3.0)
        assert total_energy.contributions[1].name == 'BaseEnergy'
        assert total_energy.contributions[1].contribution_type == 'kinetic'
        assert total_energy.contributions[1].value.to(
            'kilojoule'
        ).magnitude == approx(2.0)
        assert total_energy.contributions[2].name == 'BaseEnergy'
        assert total_energy.contributions[2].contribution_type == 'potential'
        assert total_energy.contributions[2].value.to(
            'kilojoule'
        ).magnitude == approx(1.0)

        assert outputs[1].total_forces[0].value.shape == (728, 3)
        assert outputs[1].total_forces[0].value[21][2].to('newton').magnitude == approx(
            500.0
        )
        assert outputs[2].total_forces[0].value[11].to('newton').magnitude == approx(
            500.0
        )
        force_contribution = outputs[2].total_forces[0].contributions[0]
        assert force_contribution.name == 'BaseForce'
        assert force_contribution.contribution_type == 'custom'
        assert force_contribution.value[21].to('newton').magnitude == approx(4.0)
        assert outputs[2].custom_outputs[0].m_def.name == 'CustomProperty'
        assert len(outputs[1].custom_outputs) == 1
        assert outputs[1].custom_outputs[0].name == 'custom_thermodynamic_properties'
        assert outputs[1].custom_outputs[0].value == approx(100.0)
        assert outputs[1].custom_outputs[0].unit == 'newton / angstrom ** 2'

    def test_h5md_workflow_contract(self, archive):
        workflow = archive.workflow2
        assert len(workflow.tasks) == 5
        assert workflow.method.integrator_type == 'langevin_leap_frog'
        assert workflow.method.thermodynamic_ensemble == 'NPT'
        assert workflow.method.n_steps == 20000000
        assert workflow.method.coordinate_save_frequency == 10000
        assert workflow.method.integration_timestep.to(
            'picosecond'
        ).magnitude == approx(2e-15)
        assert workflow.method.thermostat_parameters[0].thermostat_type == (
            'langevin_leap_frog'
        )
        thermostat = workflow.method.thermostat_parameters[0]
        assert thermostat.reference_temperature.magnitude == approx(300.0)
        assert thermostat.coupling_constant.to('picosecond').magnitude == approx(1.0)

        barostat = workflow.method.barostat_parameters[0]
        assert barostat.barostat_type == 'berendsen'
        assert barostat.coupling_type == 'isotropic'
        assert_approx(
            barostat.reference_pressure.to('bar').magnitude,
            np.eye(3),
        )
        assert_approx(
            barostat.coupling_constant.to('picosecond').magnitude,
            np.eye(3),
        )
        assert_approx(barostat.compressibility.to('1/bar').magnitude, np.eye(3))

        shear = workflow.method.shear_parameters[0]
        assert shear.shear_type == 'lees_edwards'
        assert_approx(
            shear.shear_rate.to('1 / picosecond').magnitude,
            [[0.0, 0.0, 0.01], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
        )
        assert workflow.method.free_energy_calculation_parameters[0].calc_type == (
            'alchemical'
        )

        results = workflow.results
        rdf = results.radial_distribution_functions
        assert len(rdf) == 3
        assert [item.label for item in rdf] == ['MOL1-MOL1', 'MOL1-MOL2', 'MOL2-MOL2']
        assert len(rdf[0].bins) == 651
        assert len(rdf[0].value) == 651
        assert rdf[0].bins[51].to('nm').magnitude == approx(0.255)
        assert rdf[0].value[51] == approx(0.284764)
        assert rdf[1].bins[51].to('nm').magnitude == approx(0.255)
        assert rdf[1].value[51] == approx(0.284764)

        msd = results.mean_squared_displacements
        assert len(msd) == 2
        assert msd[0].label == 'MOL1'
        assert msd[0].direction == 'xyz'
        assert msd[0].n_times == 51
        assert len(msd[0].times) == 51
        assert len(msd[0].value) == 51
        assert msd[0].times[10].to('ps').magnitude == approx(20.0)
        assert msd[0].value[10].to('nm**2').magnitude == approx(0.679723)
        assert msd[1].label == 'MOL2'
        assert msd[1].direction == 'xyz'
        assert msd[1].n_times == 51

        diffusion = results.diffusion_constants
        assert len(diffusion) == 2
        assert diffusion[0].label == 'MOL1'
        assert diffusion[0].value.to('nm**2/ps').magnitude == approx(1.0)
        assert diffusion[1].label == 'MOL2'
        assert diffusion[1].value.to('nm**2/ps').magnitude == approx(2.0)

        ensemble = results.ensemble_properties
        assert len(ensemble) == 4
        assert ensemble[0].label == 'bond_length_histogram'
        assert len(ensemble[0].bins_magnitude) == 10
        assert len(ensemble[0].value_magnitude) == 9
        assert ensemble[0].bins_magnitude[0] == approx(0.8)
        assert ensemble[0].bins_unit == 'angstrom'
        assert ensemble[0].value_magnitude[1] == approx(0.03076923)
        assert ensemble[1].value_magnitude == approx(-12.7)
        assert ensemble[1].value_unit == 'kilojoule / mole'
        assert ensemble[2].value_magnitude == approx(-5.2)
        assert ensemble[3].value_magnitude == approx(0.0)

        correlations = results.correlation_functions
        assert len(correlations) == 1
        assert correlations[0].label == 'velocity_autocorrelation'
        assert len(correlations[0].times) == 11
        assert len(correlations[0].value_magnitude) == 11
        assert correlations[0].times[1].to('ps').magnitude == approx(0.1)
        assert correlations[0].value_magnitude[0] == approx(1.03528105)
        assert correlations[0].value_unit == 'nanometer ** 2 / picosecond ** 2'

        for task in workflow.tasks:
            assert not any(
                hasattr(inp, 'tasks') and inp.tasks for inp in (task.inputs or [])
            )

    def test_particle_identity_is_stored_on_topology_frame_only(self, archive):
        populated = [
            index
            for index, system in enumerate(archive.data.model_system)
            if system.particle_states
        ]

        assert populated == [0]
        assert archive.data.model_system[0].particle_states[100].chemical_symbol == 'H'
