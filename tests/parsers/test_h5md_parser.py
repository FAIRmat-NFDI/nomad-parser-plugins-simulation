#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD. See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

# tests/parsers/test_h5md_parser.py

import numpy as np
import pytest
from nomad import utils
from nomad.client import normalize_all
from nomad.datamodel import EntryArchive

from nomad_simulation_parsers.parsers.h5md.parser import H5MDParser

logger = utils.get_logger(__name__)


def approx(value, abs=0, rel=1e-6):
    return pytest.approx(value, abs=abs, rel=rel)


@pytest.fixture(scope='module')
def parser():
    return H5MDParser()


def assert_h5md_header(archive: EntryArchive) -> None:
    sec_simulation = archive.data
    assert sec_simulation.program.name == 'OpenMM'
    assert sec_simulation.program.version == '-1.-1.-1'
    assert len(sec_simulation.x_h5md_version) == 2
    assert sec_simulation.x_h5md_version[1] == 0
    assert sec_simulation.x_h5md_author.name == 'Joseph F. Rudzinski'
    assert sec_simulation.x_h5md_author.email == 'joseph.rudzinski@physik.hu-berlin.de'
    assert sec_simulation.x_h5md_creator.name == 'h5py'
    assert sec_simulation.x_h5md_creator.version == '3.6.0'


def assert_systems(archive: EntryArchive) -> None:
    sec_systems = archive.data.model_system
    assert len(sec_systems) == 5
    assert np.shape(sec_systems[0].positions) == (728, 3)
    assert np.shape(sec_systems[0].velocities) == (728, 3)
    assert sec_systems[0].n_particles == 728
    assert sec_systems[0].particle_states[100].chemical_symbol == 'H'
    assert sec_systems[0].particle_states[100].label == 'H'

    assert sec_systems[2].positions[80][1].to('angstrom').magnitude == approx(28.748762)
    assert sec_systems[2].velocities[50][2].to('angstrom/ps').magnitude == approx(400.0)
    assert sec_systems[3].lattice_vectors[2][2].to(
        'angstrom'
    ).magnitude == approx(68.22318)
    assert sec_systems[3].periodic_boundary_conditions == [True, True, True]
    assert sec_systems[0].bond_list[200][0] == 198
    assert sec_systems[0].dimensionality == 3
    assert sec_systems[0].is_molecule() is False


def assert_system_hierarchy(archive: EntryArchive) -> None:
    sec_atoms_group = archive.data.model_system[0].sub_systems
    assert len(sec_atoms_group) == 4
    assert sec_atoms_group[0].particle_states == []
    assert sec_atoms_group[0].name == 'group_1ZNF'
    assert sec_atoms_group[0].branch_label == 'molecule_group'
    assert sec_atoms_group[0].composition_formula == '1ZNF(1)'
    assert sec_atoms_group[0].particle_indices[159] == 159
    assert sec_atoms_group[0].is_molecule() is True

    sec_proteins = sec_atoms_group[0].sub_systems
    assert len(sec_proteins) == 1
    assert sec_proteins[0].name == '1ZNF'
    assert sec_proteins[0].branch_label == 'molecule'
    assert (
        sec_proteins[0].composition_formula
        == 'ACE(1)TYR(1)LYS(3)CYS(2)GLY(1)LEU(2)GLU(2)ARG(3)SER(3)PHE(1)VAL(2)ALA(1)'
        'HIS(2)GLN(1)ASN(1)NH2(1)'
    )
    assert sec_proteins[0].particle_indices[400] == 400
    assert sec_proteins[0].is_molecule() is True

    sec_res_group = sec_proteins[0].sub_systems
    assert len(sec_res_group) == 16
    assert sec_res_group[13].name == 'group_ARG'
    assert sec_res_group[14].branch_label == 'monomer_group'
    assert sec_res_group[13].composition_formula == 'ARG(3)'
    assert sec_res_group[14].particle_indices[2] == 136  # TODO: check explicitly
    assert sec_res_group[14].is_molecule() is False

    sec_res = sec_res_group[13].sub_systems
    assert len(sec_res) == 3
    assert sec_res[0].name == 'ARG'
    assert sec_res[0].branch_label == 'monomer'
    assert (
        sec_res[0].composition_formula
        == 'C(1)CA(1)CB(1)CD(1)CG(1)CZ(1)H(1)HA(1)HB2(1)HB3(1)HD2(1)HD3(1)HE(1)HG2(1)'
        'HG3(1)HH11(1)HH12(1)HH21(1)HH22(1)N(1)NE(1)NH1(1)NH2(1)O(1)'
    )
    assert sec_res[0].particle_indices[10] == 120  # TODO: check explicitly
    assert sec_res[0].is_molecule() is False
    # TODO later:
    # assert sec_res[0].custom_system_attributes[0].name == 'hydrophobicity'
    # assert sec_res[0].custom_system_attributes[0].value == '0.13'
    # assert sec_res[0].custom_system_attributes[0].unit is None


def assert_outputs(archive: EntryArchive) -> None:
    sec_outputs = archive.data.outputs
    assert len(sec_outputs) == 5
    assert sec_outputs[3].step == 3
    assert sec_outputs[2].time.to('ps').magnitude == approx(2.0)

    # Temperature
    assert sec_outputs[2].temperatures[0].value.to('kelvin').magnitude == approx(300.0)

    # Energies
    assert sec_outputs[2].total_energies[0].value.to('kilojoule').magnitude == approx(
        6.0
    )
    assert sec_outputs[2].total_energies[0].contributions[0].name == 'BaseEnergy'
    assert (
        sec_outputs[2].total_energies[0].contributions[0].contribution_type == 'custom'
    )
    assert sec_outputs[2].total_energies[0].contributions[0].value.to(
        'kilojoule'
    ).magnitude == approx(3.0)
    assert sec_outputs[2].total_energies[0].contributions[1].name == 'BaseEnergy'
    assert (
        sec_outputs[2].total_energies[0].contributions[1].contribution_type == 'kinetic'
    )
    assert sec_outputs[2].total_energies[0].contributions[1].value.to(
        'kilojoule'
    ).magnitude == approx(2.0)
    assert sec_outputs[2].total_energies[0].contributions[2].name == 'BaseEnergy'
    assert (
        sec_outputs[2].total_energies[0].contributions[2].contribution_type
        == 'potential'
    )
    assert sec_outputs[2].total_energies[0].contributions[2].value.to(
        'kilojoule'
    ).magnitude == approx(1.0)

    # Forces
    assert np.shape(sec_outputs[1].total_forces[0].value) == (728, 3)
    assert sec_outputs[1].total_forces[0].value[21][2].to('newton').magnitude == approx(
        500.0
    )
    assert sec_outputs[2].total_forces[0].value[11].to('newton').magnitude == approx(
        500.0
    )
    assert sec_outputs[2].total_forces[0].contributions[0].name == 'BaseForce'
    assert sec_outputs[2].total_forces[0].contributions[0].contribution_type == 'custom'
    assert sec_outputs[2].total_forces[0].contributions[0].value[21].to(
        'newton'
    ).magnitude == approx(4.0)

    # Custom outputs
    assert sec_outputs[2].custom_outputs[0].m_def.name == 'CustomProperty'
    assert len(sec_outputs[1].custom_outputs) == 1
    assert sec_outputs[1].custom_outputs[0].name == 'custom_thermodynamic_properties'
    assert sec_outputs[1].custom_outputs[0].value == approx(100.0)
    assert sec_outputs[1].custom_outputs[0].unit == 'newton / angstrom ** 2'


def assert_md_method(sec_workflow) -> None:
    """Test MD method parameters."""
    # MD method
    assert sec_workflow.method.integrator_type == 'langevin_leap_frog'
    assert sec_workflow.method.thermodynamic_ensemble == 'NPT'
    assert sec_workflow.method.integration_timestep.to(
        'picosecond'
    ).magnitude == approx(2e-15)
    assert sec_workflow.method.n_steps == 20000000
    assert sec_workflow.method.coordinate_save_frequency == 10000
    # assert sec_workflow.method.velocity_save_frequency is None
    # assert sec_workflow.method.force_save_frequency is None
    # assert sec_workflow.method.thermodynamics_save_frequency is None


def assert_thermostats_barostats_shear(sec_workflow) -> None:
    """Test thermostat, barostat, and shear parameters."""
    # Thermostat
    sec_thermostat = sec_workflow.method.thermostat_parameters
    assert sec_thermostat[0].thermostat_type == 'langevin_leap_frog'
    assert sec_thermostat[0].reference_temperature.magnitude == approx(300.0)
    assert sec_thermostat[0].coupling_constant.to('picosecond').magnitude == approx(1.0)

    # Barostat
    sec_barostat = sec_workflow.method.barostat_parameters
    assert sec_barostat[0].barostat_type == 'berendsen'
    assert sec_barostat[0].coupling_type == 'isotropic'
    assert np.all(
        sec_barostat[0].reference_pressure.to('bar').magnitude
        == [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    )
    assert np.all(
        sec_barostat[0].coupling_constant.to('picosecond').magnitude
        == [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    )
    assert np.all(
        sec_barostat[0].compressibility.to('1/bar').magnitude
        == [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    )

    # Shear
    sec_shear = sec_workflow.method.shear_parameters
    assert sec_shear[0].shear_type == 'lees_edwards'
    assert np.all(
        sec_shear[0].shear_rate.to('1 / picosecond').magnitude
        == [[0.0, 0.0, 0.01], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
    )

    # Free energy calc parameters
    sec_free_energy = sec_workflow.method.free_energy_calculation_parameters
    assert sec_free_energy[0].calc_type == 'alchemical'
    # TODO: fix lambda index issue
    # assert sec_free_energy[0].current_lambda_index == 7
    # assert sec_free_energy[0].atom_indices.shape == (1,)
    # assert sec_free_energy[0].atom_indices[0] == 0
    # assert sec_free_energy[0].initial_state_vdw is True
    # assert sec_free_energy[0].final_state_vdw is False
    # assert sec_free_energy[0].initial_state_coloumb is False
    # assert sec_free_energy[0].final_state_coloumb is False
    # assert sec_free_energy[0].initial_state_bonded is True
    # assert sec_free_energy[0].final_state_bonded is True


def assert_radial_distribution_functions(sec_workflow_results) -> None:
    """Test radial distribution function results."""
    # MD results - RDF
    assert len(sec_workflow_results.radial_distribution_functions) == 3

    # Check first RDF (MOL1-MOL1)
    rdf_0 = sec_workflow_results.radial_distribution_functions[0]
    assert rdf_0.label == 'MOL1-MOL1'
    # assert rdf_0.error_type == 'ensemble_average'
    assert len(rdf_0.bins) == 651
    assert len(rdf_0.value) == 651
    assert rdf_0.bins[51].to('nm').magnitude == approx(0.255)
    assert rdf_0.value[51] == approx(0.284764)

    # Check second RDF (MOL1-MOL2)
    rdf_1 = sec_workflow_results.radial_distribution_functions[1]
    assert rdf_1.label == 'MOL1-MOL2'
    # assert rdf_1.error_type == 'ensemble_average'
    assert rdf_1.bins[51].to('nm').magnitude == approx(0.255)
    assert rdf_1.value[51] == approx(0.284764)

    # Check third RDF (MOL2-MOL2)
    rdf_2 = sec_workflow_results.radial_distribution_functions[2]
    assert rdf_2.label == 'MOL2-MOL2'
    # assert rdf_2.error_type == 'ensemble_average'


def assert_mean_squared_displacements(sec_workflow_results) -> None:
    """Test mean squared displacement results."""
    # MD results - MSD
    assert len(sec_workflow_results.mean_squared_displacements) == 2

    # Check first MSD (MOL1)
    msd_0 = sec_workflow_results.mean_squared_displacements[0]
    assert msd_0.label == 'MOL1'
    assert msd_0.direction == 'xyz'
    assert msd_0.n_times == 51
    assert len(msd_0.times) == 51
    assert len(msd_0.value) == 51
    assert msd_0.times[10].to('ps').magnitude == approx(20.0)
    assert msd_0.value[10].to('nm**2').magnitude == approx(0.679723)
    # Check errors array exists
    # assert len(msd_0.errors) == 51
    # assert msd_0.errors[10] == approx(0.0)

    # Check second MSD (MOL2)
    msd_1 = sec_workflow_results.mean_squared_displacements[1]
    assert msd_1.label == 'MOL2'
    assert msd_1.direction == 'xyz'
    assert msd_1.n_times == 51


def assert_diffusion_constants(sec_workflow_results) -> None:
    """Test diffusion constant results."""
    # MD results - Diffusion Constants
    assert len(sec_workflow_results.diffusion_constants) == 2

    # Check first diffusion constant (MOL1)
    diff_0 = sec_workflow_results.diffusion_constants[0]
    assert diff_0.label == 'MOL1'
    assert diff_0.value.to('nm**2/ps').magnitude == approx(1.0)

    # Check second diffusion constant (MOL2)
    diff_1 = sec_workflow_results.diffusion_constants[1]
    assert diff_1.label == 'MOL2'
    assert diff_1.value.to('nm**2/ps').magnitude == approx(2.0)


def assert_ensemble_properties(sec_workflow_results) -> None:
    """Test custom ensemble properties."""
    # Test for custom ensemble properties - these should be populated by the parser
    assert hasattr(sec_workflow_results, 'ensemble_properties')

    # Test bond length histogram (ensemble_average)
    ens_props = sec_workflow_results.ensemble_properties
    assert ens_props is not None
    assert len(ens_props) == 4

    # TODO check the application of unit factor for all of the custom obs
    bond_hist = ens_props[0]
    print(bond_hist.__dict__)
    assert bond_hist.label == 'bond_length_histogram'
    assert len(bond_hist.bins_magnitude) == 10
    assert len(bond_hist.value_magnitude) == 9
    assert bond_hist.bins_magnitude[0] == approx(0.8)
    assert bond_hist.bins_unit == 'angstrom'
    assert bond_hist.value_magnitude[1] == approx(0.03076923)

    # Test individual free energy states
    bound_prop = ens_props[1]
    assert bound_prop is not None
    assert bound_prop.value_magnitude == approx(-12.7)
    assert bound_prop.value_unit == 'kilojoule / mole'

    intermediate_prop = ens_props[2]
    assert intermediate_prop is not None
    assert intermediate_prop.value_magnitude == approx(-5.2)
    assert bound_prop.value_unit == 'kilojoule / mole'

    unbound_prop = ens_props[3]
    assert unbound_prop is not None
    assert unbound_prop.value_magnitude == approx(0.0)
    assert bound_prop.value_unit == 'kilojoule / mole'


def assert_correlation_functions(sec_workflow_results) -> None:
    """Test correlation function results."""
    # Test for custom correlation functions
    assert hasattr(sec_workflow_results, 'correlation_functions')

    # Test velocity autocorrelation function
    corr_funcs = sec_workflow_results.correlation_functions
    assert corr_funcs is not None
    assert len(corr_funcs) == 1

    vacf = corr_funcs[0]
    print(vacf.__dict__)
    assert vacf.label == 'velocity_autocorrelation'
    assert len(vacf.times) == 11
    assert len(vacf.value_magnitude) == 11
    assert vacf.times[1].to('ps').magnitude == approx(0.1)
    assert vacf.value_magnitude[0] == approx(1.03528105)
    assert vacf.value_unit == 'nanometer ** 2 / picosecond ** 2'

    # TODO Add Rg tests


def assert_workflow(archive: EntryArchive) -> None:
    """Test workflow data by delegating to specialized assertion functions."""
    sec_workflow = archive.workflow2
    sec_workflow_results = sec_workflow.results

    assert_md_method(sec_workflow)
    assert_thermostats_barostats_shear(sec_workflow)
    assert_radial_distribution_functions(sec_workflow_results)
    assert_mean_squared_displacements(sec_workflow_results)
    assert_diffusion_constants(sec_workflow_results)
    assert_ensemble_properties(sec_workflow_results)
    assert_correlation_functions(sec_workflow_results)

    # MD results
    # sec_workflow_results = sec_workflow.results
    # assert len(sec_workflow_results.ensemble_properties) == 1
    # ensemble_property_0 = sec_workflow_results.ensemble_properties[0]
    # assert ensemble_property_0.label == 'diffusion_constants'
    # assert ensemble_property_0.error_type == 'Pearson_correlation_coefficient'
    # assert len(ensemble_property_0.ensemble_property_values) == 2
    # assert ensemble_property_0.ensemble_property_values[1].label == 'MOL2'
    # assert ensemble_property_0.ensemble_property_values[1].errors == 0.95
    # assert ensemble_property_0.ensemble_property_values[1].value_magnitude == 2.0
    # assert (
    #     ensemble_property_0.ensemble_property_values[1].value_unit
    #     == 'nanometer ** 2 / picosecond'
    # )
    # assert bound_prop is not None
    # assert bound_prop.value_magnitude == approx(-12.7)

    # intermediate_prop = next(
    #     (prop for prop in free_energy_props if 'intermediate' in prop.label), None
    # )
    # assert intermediate_prop is not None
    # assert intermediate_prop.value_magnitude == approx(-5.2)

    # unbound_prop = next(
    #     (prop for prop in free_energy_props if 'unbound' in prop.label), None
    # )
    # assert unbound_prop is not None
    # assert unbound_prop.value_magnitude == approx(0.0)


def test_md(parser):
    archive = EntryArchive()
    parser.parse(
        'tests/data/h5md/test_traj_openmm_reduced-SOL_5frames_07-10-25.h5',
        archive,
        None,
    )
    normalize_all(archive, logger=logger)

    assert_h5md_header(archive)
    assert_systems(archive)
    assert_system_hierarchy(archive)
    assert_outputs(archive)
    assert_workflow(archive)
