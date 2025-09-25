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


import os
import tempfile

import numpy as np
import pytest
from nomad.utils import get_logger

from nomad_simulation_parsers.parsers.lammps.parser import LammpsParser
from nomad_simulation_parsers.parsers.lammps.trajectory_parsers import (
    TrajParser,
    TrajParsers,
    XYZTrajParser,
)
from nomad_simulation_parsers.parsers.utils.mdanalysisparser import MDAnalysisParser

LOGGER = get_logger(__name__)


@pytest.fixture(scope='module')
def parser():
    return LammpsParser()


# TODO: add tests for file_parsers functions
# Tests for get_unit() function with different unit types
# Tests for DataParser: regex patterns and section parsing
# Tests for LogParser: command extraction and thermodynamic data parsing
# Tests for file discovery methods (get_traj_files, get_data_files)


# Tests for TrajParser, XYZTrajParser, TrajParsers classes
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


def test_traj_xyz():
    # Synthetic XYZ trajectory content
    xyz_content = """5
Atoms. Timestep: 0
1 4.39861 0.0809956 -1.6196
2 3.65138 0.778109 -1.97822
2 4.72189 -0.655793 -2.40238
2 5.23117 0.689443 -1.27747
2 3.94587 -0.457468 -0.7756
5
Atoms. Timestep: 400
1 4.17634 0.0441698 -1.4592
2 3.33775 0.267888 -2.08495
2 4.74748 -0.845205 -1.78471
2 4.8507 0.87915 -1.4652
2 3.77509 -0.143827 -0.474483
"""

    # XYZTrajParser has no _file_handler attribute, needs a temporary file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.xyz', delete=False) as f:
        f.write(xyz_content)
        temp_file = f.name

    try:
        xyz_parser = XYZTrajParser()
        xyz_parser.mainfile = temp_file
        xyz_parser.logger = LOGGER
        xyz_parser.init_quantities()

        parsers = TrajParsers([xyz_parser])
        n_frames = parsers.eval('n_frames')
        assert n_frames == 2
        positions = xyz_parser.get_positions(1)
        assert positions[2][1] == pytest.approx(-0.845205)

    finally:
        os.unlink(temp_file)


def test_unwrapped_pos():
    # 1_xyz dataset (CG), file type 'custom' -> TrajParser
    traj_parser = TrajParser()
    traj_parser.mainfile = 'tests/data/lammps/1_xyz_files/pos_vel.xyz'
    traj_parser.init_quantities()
    # TODO: add assertion for calculation
    positions = traj_parser.get_positions(1)
    assert positions[452][2] == pytest.approx(5.99898)
    velocities = traj_parser.get_velocities(2)
    assert velocities[457][-2] == pytest.approx(-0.928553)


# TODO Fix dealing with multiple output files (positions and velocities in separate files)
# TODO with archive_to_universe function, then add back in this test


# Tests for: LammpsArchiveWriter, LammpsParser (integration tests)
def test_traj_dcd():
    dcd_parser = MDAnalysisParser(topology_format='DATA', format='DCD')
    dcd_parser.mainfile = 'tests/data/lammps/methane_dcd/data.64xmethane_from_restart'
    dcd_parser.auxilliary_files = ['tests/data/lammps/methane_dcd/64xmethane-nvt.dcd']
    dcd_parser.logger = LOGGER
    dcd_parser.parse()
    # TODO: add assertion for calculation
    positions = dcd_parser.get_positions(56)
    assert np.shape(positions) == (320, 3)
    labels = dcd_parser.get_atom_labels(107)
    assert len(labels) == 320


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

# TODO Add tests that use the full parser fixture and test end-to-end parsing
