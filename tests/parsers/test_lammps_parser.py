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

import pytest
import numpy as np

from nomad.datamodel import EntryArchive
from nomad.utils import get_logger
from nomad_simulation_parsers.parsers.lammps.parser import (
    LammpsParser,
    TrajParser,
    TrajParsers,
)

LOGGER = get_logger(__name__)


@pytest.fixture(scope='module')
def parser():
    return LammpsParser()


# TODO: Extend test to cover all relevant LAMMPS box styles
@pytest.mark.parametrize(
    'description, content, expected_pbc, expected_cell',
    [
        (
            'Orthogonal cell, all dimensions periodic',
            """
ITEM: BOX BOUNDS pp pp pp
-13.4569 13.25
-14.6313 14.1743
-12.4476 12.4476
        """,
            [True, True, True],
            np.diag(
                [
                    13.25 - (-13.4569),
                    14.1743 - (-14.6313),
                    12.4476 - (-12.4476),
                ]
            ),
        ),
        # (
        #     'Description',
        #     """
        # ITEM: BOX BOUNDS
        #         """,
        #     [],
        #     np.array([[], [], []]),
        # ),
    ],
)
def test_pbc_cell_extraction(description, content, expected_pbc, expected_cell):
    parser = TrajParser()
    parser.mainfile = 'dummy'
    parser._file_handler = content.encode('utf-8')
    parser.init_quantities()

    parsers = TrajParsers([parser])
    pbc_cell = parsers.eval('pbc_cell')
    assert len(pbc_cell) == 1, f'{description} - pbc_cell not extracted'

    pbc, cell = pbc_cell[0]
    assert pbc == expected_pbc, f'{description} - wrong PBC'
    assert cell == pytest.approx(expected_cell), f'{description} - wrong cell'


# TODO: re-include output testing alongside migration of parser functionalities
# def test_nvt(parser):
#     archive = EntryArchive()
#     parser.parse(
#         'tests/data/lammps/hexane_cyclohexane/log.hexane_cyclohexane_nvt', archive, None
#     )

#     sec_run = archive.run[0]
#     assert sec_run.program.version == '14 May 2016'

#     sec_workflow = archive.workflow2
#     assert sec_workflow.m_def.name == 'MolecularDynamics'
#     assert sec_workflow.method.thermodynamic_ensemble == 'NVT'
#     assert sec_workflow.method.integrator_type == 'velocity_verlet'
#     assert sec_workflow.method.integration_timestep.magnitude == 2.5e-16
#     assert sec_workflow.method.integration_timestep.units == 'second'
#     assert sec_workflow.method.n_steps == 80000
#     assert sec_workflow.method.coordinate_save_frequency == 400
#     assert sec_workflow.method.thermodynamics_save_frequency == 400
#     assert sec_workflow.method.thermostat_parameters[0].thermostat_type == 'nose_hoover'
#     assert (
#         sec_workflow.method.thermostat_parameters[0].reference_temperature.magnitude
#         == 300.0
#     )
#     assert (
#         sec_workflow.method.thermostat_parameters[0].reference_temperature.units
#         == 'kelvin'
#     )
#     assert (
#         sec_workflow.method.thermostat_parameters[0].coupling_constant.magnitude
#         == 2.5e-14
#     )
#     assert (
#         sec_workflow.method.thermostat_parameters[0].coupling_constant.units == 'second'
#     )

#     sec_method = sec_run.method[0]
#     assert len(sec_method.force_field.model[0].contributions) == 3
#     assert sec_method.force_field.model[0].contributions[1].type == 'bond'
#     assert sec_method.force_field.model[0].contributions[1].n_interactions == 666
#     assert sec_method.force_field.model[0].contributions[1].n_atoms == 2
#     assert sec_method.force_field.model[0].contributions[1].atom_indices[100, 1] == 103
#     assert sec_method.force_field.model[0].contributions[1].atom_labels[350, 0] == '1'
#     assert (
#         sec_method.force_field.force_calculations.coulomb_cutoff.magnitude
#         == 1.2000000000000002e-08
#     )
#     assert sec_method.force_field.force_calculations.coulomb_cutoff.units == 'meter'
#     assert (
#         sec_method.force_field.force_calculations.neighbor_searching.neighbor_update_frequency
#         == 10
#     )

#     sec_system = sec_run.system
#     assert len(sec_system) == 201
#     assert sec_system[5].atoms.lattice_vectors[1][1].magnitude == approx(2.24235e-09)
#     assert False not in sec_system[0].atoms.periodic
#     assert sec_system[80].atoms.labels[91:96] == ['H', 'H', 'H', 'C', 'C']
#     assert sec_system[0].atoms.bond_list[200, 0] == 194

#     sec_scc = sec_run.calculation
#     assert len(sec_scc) == 201
#     assert sec_scc[21].energy.current.value.magnitude == approx(8.86689197e-18)
#     assert sec_scc[56].pressure.magnitude == approx(-77642135.4975)
#     assert sec_scc[103].temperature.magnitude == 291.4591
#     assert sec_scc[11].step == 4400
#     assert len(sec_scc[1].energy.contributions) == 9
#     assert sec_scc[112].energy.contributions[8].kind == 'kspace long range'
#     assert sec_scc[96].energy.contributions[2].value.magnitude == approx(1.19666271e-18)
#     assert sec_scc[47].energy.contributions[4].value.magnitude == approx(1.42166035e-18)
#     assert sec_scc[75].time_physical.magnitude == approx(83.56332225)
#     assert sec_scc[112].time_calculation.magnitude == approx(1.2351 / 400)

#     assert (
#         sec_run.x_lammps_section_control_parameters[0].x_lammps_inout_control_atomstyle
#         == 'full'
#     )


# def test_thermo_format(parser):
#     archive = EntryArchive()
#     parser.parse(
#         'tests/data/lammps/1_methyl_naphthalene/log.1_methyl_naphthalene', archive, None
#     )

#     sec_sccs = archive.run[0].calculation
#     assert len(sec_sccs) == 301
#     assert sec_sccs[98].energy.total.value.magnitude == approx(1.45322428e-17)

#     assert len(archive.run[0].system) == 4


def test_traj_xyz(parser):
    archive = EntryArchive()
    parser.parse(
        'tests/data/lammps/methane_xyz/log.methane_nvt_traj_xyz_thermo_style_custom',
        archive,
        LOGGER,
    )
    sec_systems = archive.data.model_system
    assert len(sec_systems) == 201
    assert sec_systems[13].positions[7][0].magnitude == pytest.approx(-8.00436e-10)


def test_traj_dcd(parser):
    archive = EntryArchive()
    parser.parse(
        'tests/data/lammps/methane_dcd/log.methane_nvt_traj_dcd_thermo_style_custom',
        archive,
        LOGGER,
    )
    # TODO: add assertion for calculation
    # assert len(archive.run[0].calculation) == 201
    sec_systems = archive.data.model_system
    assert np.shape(sec_systems[56].positions) == (320, 3)
    assert (
        len(
            [
                particle_state.label
                for particle_state in sec_systems[107].particle_states
            ]
        )
        == 320
    )


def test_unwrapped_pos(parser):
    archive = EntryArchive()
    parser.parse('tests/data/lammps/1_xyz_files/log.lammps', archive, LOGGER)
    # TODO: add assertion for calculation
    # ? Where has "calculation" moved to now, which section of the parser?
    # assert len(archive.run[0].calculation) == 101
    sec_systems = archive.data.model_system
    assert sec_systems[1].positions[452][2].magnitude == pytest.approx(5.99898)
    assert sec_systems[2].velocities[457][-2].magnitude == pytest.approx(-0.928553)


# ! Positions and velocities are in separate files. MDAnalysis-parser fails to create
# ! universe during parsing attempt of velocities file.
# ? Solution: somehow identify velocities-only file and declare it an auxillary file,
# ? use universe generated with positions?
# TODO Fix dealing with multiple output files with archive_to_universe function, then add back in this test
# def test_multiple_dump(parser):
#     archive = EntryArchive()
#     parser.parse('tests/data/lammps/2_xyz_files/log.lammps', archive, None)

#     sec_systems = archive.run[0].system
#     assert len(sec_systems) == 101
#     assert sec_systems[2].atoms.positions[468][0].magnitude == approx(3.00831)
#     assert sec_systems[-1].atoms.velocities[72][1].magnitude == approx(-4.61496)  # JFR - universe cannot be built without positions


def test_systems(parser) -> None:
    archive = EntryArchive()
    parser.parse(
        'tests/data/lammps/methane_dcd/log.methane_nvt_traj_dcd_thermo_style_custom',
        archive,
        LOGGER,
    )
    sec_systems = archive.data.model_system
    assert len(sec_systems) == 4
    assert np.shape(sec_systems[0].positions) == (1134, 3)
    # assert np.shape(sec_systems[0].velocities) == (1134, 3)
    assert sec_systems[0].n_particles == 1134
    assert sec_systems[0].particle_states[100].chemical_symbol == 'H'
    assert sec_systems[0].particle_states[100].label == 'H'

    assert sec_systems[2].positions[567][1].to('angstrom').magnitude == pytest.approx(
        -58847500000.0
    )
    # TODO: Atomic test data does not have velocities, update testing!
    # assert sec_systems[idx].velocities[idx][idx].to(
    #     'angstrom/ps'
    # ).magnitude == pytest.approx(target_float)
    assert sec_systems[3].cell[0].lattice_vectors[2][2].to(
        'angstrom'
    ).magnitude == pytest.approx(214680000000.0)
    assert sec_systems[3].cell[0].periodic_boundary_conditions == [True, True, True]
    assert sec_systems[0].bond_list[200][0] == np.array([189, 192])
    # TODO: fiugre out why dimensionality isn't set
    # assert sec_systems[0].dimensionality == 3
    assert sec_systems[0].is_molecule() is False


def assert_system_hierarchy(archive: EntryArchive) -> None:
    sec_atoms_group = archive.data.model_system[0].sub_systems
    assert len(sec_atoms_group) == 4
    assert sec_atoms_group[0].particle_states == []
    # TODO comment back in once nested fix is in release
    # assert sec_atoms_group[0].cell == []
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


# TODO: update structure to new schema
# def test_md_atomsgroup(parser):
#     archive = EntryArchive()
#     parser.parse(
#         'tests/data/lammps/polymer_melt/Emin/log.step4.0_minimization', archive, None
#     )

#     # sec_run = archive.run[0]
#     sec_systems = archive.data.model_system

#     # assert len(sec_systems[0].atoms_group) == 1
#     # assert len(sec_systems[0].atoms_group[0].atoms_group) == 100

#     # assert sec_systems[0].atoms_group[0].label == 'group_0'
#     # assert sec_systems[0].atoms_group[0].type == 'molecule_group'
#     # assert sec_systems[0].atoms_group[0].index == 0
#     # assert sec_systems[0].atoms_group[0].composition_formula == '0(100)'
#     # assert sec_systems[0].atoms_group[0].n_atoms == 7200
#     # assert sec_systems[0].atoms_group[0].atom_indices[5] == 5
#     # assert sec_systems[0].atoms_group[0].is_molecule is False

#     # assert sec_systems[0].atoms_group[0].atoms_group[52].label == '0'
#     # assert sec_systems[0].atoms_group[0].atoms_group[52].type == 'molecule'
#     # assert sec_systems[0].atoms_group[0].atoms_group[52].index == 52
#     # assert (
#     #     sec_systems[0].atoms_group[0].atoms_group[52].composition_formula
#     #     == '1(1)2(1)3(1)4(1)5(1)6(1)7(1)8(1)9(1)10(1)'
#     # )
#     # assert sec_systems[0].atoms_group[0].atoms_group[52].n_atoms == 72
#     # assert sec_systems[0].atoms_group[0].atoms_group[52].atom_indices[8] == 3752
#     # assert sec_systems[0].atoms_group[0].atoms_group[52].is_molecule is True

#     # assert (
#     #     sec_systems[0].atoms_group[0].atoms_group[76].atoms_group[7].label == 'group_8'
#     # )
#     # assert (
#     #     sec_systems[0].atoms_group[0].atoms_group[76].atoms_group[7].type
#     #     == 'monomer_group'
#     # )
#     # assert sec_systems[0].atoms_group[0].atoms_group[76].atoms_group[7].index == 7
#     # assert (
#     #     sec_systems[0].atoms_group[0].atoms_group[76].atoms_group[7].composition_formula
#     #     == '8(1)'
#     # )
#     # assert sec_systems[0].atoms_group[0].atoms_group[76].atoms_group[7].n_atoms == 7
#     # assert (
#     #     sec_systems[0].atoms_group[0].atoms_group[76].atoms_group[7].atom_indices[5]
#     #     == 5527
#     # )
#     # assert (
#     #     sec_systems[0].atoms_group[0].atoms_group[76].atoms_group[7].is_molecule
#     #     is False
#     # )

#     # assert (
#     #     sec_systems[0]
#     #     .atoms_group[0]
#     #     .atoms_group[76]
#     #     .atoms_group[7]
#     #     .atoms_group[0]
#     #     .label
#     #     == '8'
#     # )
#     # assert (
#     #     sec_systems[0].atoms_group[0].atoms_group[76].atoms_group[7].atoms_group[0].type
#     #     == 'monomer'
#     # )
#     # assert (
#     #     sec_systems[0]
#     #     .atoms_group[0]
#     #     .atoms_group[76]
#     #     .atoms_group[7]
#     #     .atoms_group[0]
#     #     .index
#     #     == 0
#     # )
#     # assert (
#     #     sec_systems[0]
#     #     .atoms_group[0]
#     #     .atoms_group[76]
#     #     .atoms_group[7]
#     #     .atoms_group[0]
#     #     .composition_formula
#     #     == '1(4)4(2)6(1)'
#     # )
#     # assert (
#     #     sec_systems[0]
#     #     .atoms_group[0]
#     #     .atoms_group[76]
#     #     .atoms_group[7]
#     #     .atoms_group[0]
#     #     .n_atoms
#     #     == 7
#     # )
#     # assert (
#     #     sec_systems[0]
#     #     .atoms_group[0]
#     #     .atoms_group[76]
#     #     .atoms_group[7]
#     #     .atoms_group[0]
#     #     .atom_indices[5]
#     #     == 5527
#     # )
#     # assert (
#     #     sec_systems[0]
#     #     .atoms_group[0]
#     #     .atoms_group[76]
#     #     .atoms_group[7]
#     #     .atoms_group[0]
#     #     .is_molecule
#     #     is False
#     # )


# TODO re-activate when migrating workflow parsing
# def test_geometry_optimization(parser):
#     archive = EntryArchive()
#     parser.parse(
#         'tests/data/lammps/polymer_melt/Emin/log.step4.0_minimization', archive, None
#     )

#     sec_workflow = archive.workflow2

#     assert sec_workflow.method.type == 'atomic'
#     assert sec_workflow.method.method == 'polak_ribiere_conjugant_gradient'

#     assert (
#         sec_workflow.method.convergence_tolerance_energy_difference.magnitude
#         == approx(0.0)
#     )
#     assert sec_workflow.method.convergence_tolerance_energy_difference.units == 'joule'
#     assert sec_workflow.results.final_energy_difference.magnitude == approx(0.0)
#     assert sec_workflow.results.final_energy_difference.units == 'joule'

#     assert sec_workflow.method.convergence_tolerance_force_maximum.magnitude == approx(
#         100
#     )
#     assert sec_workflow.method.convergence_tolerance_force_maximum.units == 'newton'

#     assert sec_workflow.results.final_force_maximum.magnitude == approx(5091750000.0)
#     assert sec_workflow.results.final_force_maximum.units == 'newton'

#     assert sec_workflow.method.optimization_steps_maximum == 10000
#     assert sec_workflow.results.optimization_steps == 160
#     assert len(sec_workflow.results.energies) == 159
#     assert sec_workflow.results.energies[14].magnitude == approx(6.931486093999211e-17)
#     assert sec_workflow.results.energies[14].units == 'joule'
#     assert len(sec_workflow.results.steps) == 159
#     assert sec_workflow.results.steps[22] == 1100
