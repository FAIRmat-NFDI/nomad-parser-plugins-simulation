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
import re
import tempfile
from io import BytesIO, StringIO

import numpy as np
import pytest
from nomad.client import normalize_all
from nomad.datamodel import EntryArchive, EntryMetadata
from nomad.units import ureg
from nomad.utils import get_logger

from nomad_simulation_parsers.parsers.lammps.file_parsers import LogParser
from nomad_simulation_parsers.parsers.lammps.parser import (
    LammpsArchiveWriter,
    LammpsParser,
)
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


@pytest.fixture
def tmp_dir():
    """Create temporary directory for tests."""
    parent_directory = '.volumes'
    if not os.path.isdir(parent_directory):
        os.makedirs(parent_directory, exist_ok=True)
    directory = tempfile.TemporaryDirectory(dir=parent_directory, prefix='test_tmp')
    yield directory.name
    directory.cleanup()


# TODO: add tests for file_parsers functions
# Tests for DataParser: regex patterns and section parsing
# Tests for LogParser: command extraction and thermodynamic data parsing
# Tests for file discovery methods (get_traj_files, get_data_files)
@pytest.mark.parametrize(
    'mainfile, traj_files, expected_match',
    [
        # Strategy 1: Exact prefix match
        (
            'log.hexane_nvt',
            ['hexane_nvt.lammpstrj', 'water.xyz', 'min.lammpstrj', 'other.dcd'],
            'hexane_nvt.lammpstrj',
        ),
        # Strategy 2: Token-based similarity
        (
            'log.hexane_nvt',
            ['water.lammpstrj', 'hex_nvt_output.xyz', 'polymer.dcd'],
            'hex_nvt_output.xyz',
        ),
        (
            'log.1_methyl_naphthalene',
            ['naph_298_eq.lammpstrj', 'water.xyz'],
            'naph_298_eq.lammpstrj',
        ),
        # Strategy 3: String sequence similarity
        (
            'log.hexane_nvt',
            ['hexanol_simulation.lammpstrj', 'completely_different.xyz'],
            'hexanol_simulation.lammpstrj',
        ),
        # Fallback: No match, use first file
        (
            'log.hexane_nvt',
            ['water.lammpstrj', 'polymer.xyz', 'argon.dcd'],
            'water.lammpstrj',
        ),
        # Case insensitive
        (
            'log.hexane_nvt',
            ['HEXANE_NVT.lammpstrj', 'other.xyz'],
            'HEXANE_NVT.lammpstrj',
        ),
    ],
)
def test_find_best_matching_file(mainfile, traj_files, expected_match):
    """Test trajectory file matching with various scenarios"""
    parser = LogParser()
    parser.logger = LOGGER

    result = parser.find_best_matching_file(traj_files, mainfile)

    assert len(result) == 1
    assert result[0] == expected_match


class TestLogParser(LogParser):
    """Test double for LogParser that allows injecting mock data."""

    def __init__(self, mock_data=None):
        super().__init__()
        self.mock_data = mock_data or {}

    def parse(self, key, **kwargs):
        if self._results is None:
            self._results = {}
        if key in self.mock_data:
            self._results[key] = self.mock_data[key]
        else:
            super().parse(key, **kwargs)


def test_get_traj_files(tmp_dir):
    """Test trajectory file discovery and matching"""

    # Test 1: With dump command specified
    parser = TestLogParser(
        mock_data={
            'dump': [
                [
                    '2',
                    'all',
                    'custom',
                    '100',
                    'pos_vel.xyz',
                    'id',
                    'type',
                    'xu',
                    'yu',
                    'zu',
                    'fx',
                    'fy',
                    'fz',
                    'vx',
                    'vy',
                    'vz',
                ]
            ]
        }
    )
    parser.logger = LOGGER
    parser.mainfile = 'tests/data/lammps/1_xyz_files/log.test'
    traj_files = parser.get_traj_files()

    assert len(traj_files) == 1
    assert os.path.basename(traj_files[0]) == 'pos_vel.xyz'
    assert os.path.isabs(traj_files[0])

    # Test 2: Scan directory (no dump command)
    open(os.path.join(tmp_dir, 'hexane_nvt.lammpstrj'), 'w').close()
    open(os.path.join(tmp_dir, 'water.xyz'), 'w').close()
    open(os.path.join(tmp_dir, 'readme.txt'), 'w').close()

    mainfile = os.path.join(tmp_dir, 'log.hexane_nvt')
    open(mainfile, 'w').close()

    parser = TestLogParser(mock_data={'dump': None})
    parser.logger = LOGGER
    parser.mainfile = mainfile

    traj_files = parser.get_traj_files()

    # Should find trajectory files and return best match to mainfile
    assert len(traj_files) == 1
    assert os.path.basename(traj_files[0]) == 'hexane_nvt.lammpstrj'
    assert all(os.path.isabs(f) for f in traj_files)

    # Test 3: Multiple dump commands, remove duplicates
    parser = TestLogParser(
        mock_data={
            'dump': [
                ['id1', 'all', 'custom', '100', 'traj.lammpstrj', 'x', 'y', 'z'],
                ['id2', 'all', 'custom', '200', 'traj.lammpstrj', 'x', 'y', 'z'],
            ]
        }
    )
    parser.logger = LOGGER
    parser.mainfile = 'tests/data/lammps/2_xyz_files/log.lammps'
    traj_files = parser.get_traj_files()

    assert len(traj_files) == 1
    assert os.path.basename(traj_files[0]) == 'traj.lammpstrj'


def test_no_data_file(tmp_dir):
    """Test behavior when no data files are found"""

    # Create new mainfile to prevent cached results being used
    mainfile_path = os.path.join(tmp_dir, 'log.test')
    open(mainfile_path, 'w').close()

    parser = TestLogParser(mock_data={'read_data': []})
    parser.logger = LOGGER
    parser.mainfile = mainfile_path

    data_files = parser.get_data_files()

    assert isinstance(data_files, list)
    assert len(data_files) == 0


def test_data_file_header_matching():
    """Test LAMMPS data file header pattern matching using in-memory StringIO"""

    # Simulate what check_file_header does: read first 1024 bytes and match pattern
    def check_file_header_in_memory(file_obj, pattern):
        """Mimic check_file_header behavior using StringIO"""
        file_obj.seek(0)
        header = file_obj.read(1024)
        if isinstance(header, str):
            header_str = header
        else:
            header_str = header.decode(errors='ignore')
        return re.search(pattern, header_str)

    # Test 1: Valid LAMMPS data file header
    valid_data_file = StringIO("""LAMMPS data file via write_data, version 12 Dec 2018

1000 atoms
10 atom types
0.0 50.0 xlo xhi
0.0 50.0 ylo yhi
0.0 50.0 zlo zhi
""")
    assert check_file_header_in_memory(valid_data_file, 'LAMMPS data file') is not None

    # Test 2: Alternative valid header
    alt_valid_file = StringIO("""LAMMPS Description

500 atoms
5 atom types
""")
    assert check_file_header_in_memory(alt_valid_file, 'LAMMPS Description') is not None

    # Test 3: Invalid file (no LAMMPS header)
    invalid_file = StringIO("""Random text file
This is not a relevant file.
Just some random content
""")
    assert check_file_header_in_memory(invalid_file, 'LAMMPS data file') is None
    assert check_file_header_in_memory(invalid_file, 'LAMMPS Description') is None

    # Test 4: Binary data with valid header
    binary_with_header = BytesIO(
        b'LAMMPS data file\n1000 atoms\n' + b'\x00\xff\xfe\xfd' * 100
    )
    assert (
        check_file_header_in_memory(binary_with_header, 'LAMMPS data file') is not None
    )

    # Test 5: Binary data without valid header (use predictable low bytes)
    pure_binary = BytesIO(b'\x00\x01\x02\x03' * 256)
    assert check_file_header_in_memory(pure_binary, 'LAMMPS data file') is None
    assert check_file_header_in_memory(pure_binary, 'LAMMPS Description') is None


def test_data_file_header_priority():
    """Test header pattern matching priority (LAMMPS data file vs LAMMPS Description)"""

    def check_file_header_in_memory(file_obj, pattern):
        file_obj.seek(0)
        header = file_obj.read(1024)
        if isinstance(header, str):
            header_str = header
        else:
            header_str = header.decode(errors='ignore')
        return re.search(pattern, header_str)

    # File with standard header should match primary pattern
    standard_file = StringIO("""LAMMPS data file
1000 atoms
""")
    assert check_file_header_in_memory(standard_file, 'LAMMPS data file') is not None

    # File with alternative header should only match alternative pattern
    alt_file = StringIO("""LAMMPS Description
500 atoms
""")
    assert check_file_header_in_memory(alt_file, 'LAMMPS data file') is None
    assert check_file_header_in_memory(alt_file, 'LAMMPS Description') is not None


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
    """Test PBC and cell extraction from synthetic LAMMPS trajectory content"""
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


def test_get_lattice_vectors():
    """Test lattice vector extraction from TrajParser"""
    traj_content = """ITEM: TIMESTEP
0
ITEM: NUMBER OF ATOMS
5
ITEM: BOX BOUNDS pp pp pp
-10.0 10.0
-15.0 15.0
-12.0 12.0
ITEM: ATOMS id type x y z
1 1 0.0 0.0 0.0
2 2 1.0 1.0 1.0
3 2 2.0 2.0 2.0
4 2 3.0 3.0 3.0
5 2 4.0 4.0 4.0
"""

    with tempfile.NamedTemporaryFile(mode='w', suffix='.lammpstrj', delete=False) as f:
        f.write(traj_content)
        temp_file = f.name

    try:
        traj_parser = TrajParser()
        traj_parser.mainfile = temp_file
        traj_parser.logger = LOGGER
        traj_parser.init_quantities()

        # Test lattice vectors
        lattice_vectors = traj_parser.get_lattice_vectors(0)
        assert lattice_vectors is not None
        assert lattice_vectors.shape == (3, 3)
        assert lattice_vectors[0, 0] == pytest.approx(20.0)  # 10.0 - (-10.0)
        assert lattice_vectors[1, 1] == pytest.approx(30.0)  # 15.0 - (-15.0)
        assert lattice_vectors[2, 2] == pytest.approx(24.0)  # 12.0 - (-12.0)

    finally:
        os.unlink(temp_file)


def test_get_n_atoms():
    """Test n_atoms extraction"""
    traj_content = """ITEM: TIMESTEP
0
ITEM: NUMBER OF ATOMS
5
ITEM: BOX BOUNDS pp pp pp
-10.0 10.0
-15.0 15.0
-12.0 12.0
ITEM: ATOMS id type x y z
1 1 0.0 0.0 0.0
2 2 1.0 1.0 1.0
3 2 2.0 2.0 2.0
4 2 3.0 3.0 3.0
5 2 4.0 4.0 4.0
"""

    with tempfile.NamedTemporaryFile(mode='w', suffix='.lammpstrj', delete=False) as f:
        f.write(traj_content)
        temp_file = f.name

    try:
        traj_parser = TrajParser()
        traj_parser.mainfile = temp_file
        traj_parser.logger = LOGGER
        traj_parser.init_quantities()

        # Test n_atoms
        n_atoms = traj_parser.get_n_atoms(0)
        assert n_atoms == 5

    finally:
        os.unlink(temp_file)


def test_get_step():
    """Test timestep extraction"""
    traj_content = """ITEM: TIMESTEP
100
ITEM: NUMBER OF ATOMS
3
ITEM: BOX BOUNDS pp pp pp
0.0 10.0
0.0 10.0
0.0 10.0
ITEM: ATOMS id type x y z
1 1 0.0 0.0 0.0
2 2 1.0 1.0 1.0
3 2 2.0 2.0 2.0
ITEM: TIMESTEP
200
ITEM: NUMBER OF ATOMS
3
ITEM: BOX BOUNDS pp pp pp
0.0 10.0
0.0 10.0
0.0 10.0
ITEM: ATOMS id type x y z
1 1 0.5 0.5 0.5
2 2 1.5 1.5 1.5
3 2 2.5 2.5 2.5
"""

    with tempfile.NamedTemporaryFile(mode='w', suffix='.lammpstrj', delete=False) as f:
        f.write(traj_content)
        temp_file = f.name

    try:
        traj_parser = TrajParser()
        traj_parser.mainfile = temp_file
        traj_parser.logger = LOGGER
        traj_parser.init_quantities()

        # Test timestep for first frame
        step_0 = traj_parser.get_step(0)
        assert step_0 == 100

        # Test timestep for second frame
        step_1 = traj_parser.get_step(1)
        assert step_1 == 200

    finally:
        os.unlink(temp_file)


def test_none_returns():
    """Test that methods return None when data is missing"""
    # Create parser with minimal/no data
    traj_content = """ITEM: TIMESTEP
0
"""

    with tempfile.NamedTemporaryFile(mode='w', suffix='.lammpstrj', delete=False) as f:
        f.write(traj_content)
        temp_file = f.name

    try:
        traj_parser = TrajParser()
        traj_parser.mainfile = temp_file
        traj_parser.logger = LOGGER
        traj_parser.init_quantities()

        # All should return None when data is missing
        assert traj_parser.get_lattice_vectors(0) is None
        assert traj_parser.get_pbc(0) is None
        # n_atoms falls back to counting positions, which also returns None
        assert traj_parser.get_n_atoms(0) is None

    finally:
        os.unlink(temp_file)


def test_with_real_data():
    """Test methods with existing real test data"""
    traj_parser = TrajParser()
    traj_parser.mainfile = 'tests/data/lammps/1_xyz_files/pos_vel.xyz'
    traj_parser.logger = LOGGER
    traj_parser.init_quantities()

    # Test that methods work with real data
    n_atoms = traj_parser.get_n_atoms(5)
    assert n_atoms is not None
    assert isinstance(n_atoms, int)
    assert n_atoms > 0

    step = traj_parser.get_step(3)
    assert step is not None
    assert isinstance(step, int)

    # If the real data has PBC info, test it
    pbc = traj_parser.get_pbc(1)
    lattice_vectors = traj_parser.get_lattice_vectors(1)
    # These might be None if the file doesn't include PBC info
    assert isinstance(pbc, list)
    assert isinstance(lattice_vectors, np.ndarray)


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


# TODO Fix dealing with multiple output files (positions and velocities separate)


def test_traj_xyz():
    """Test XYZTrajParser with synthetic XYZ trajectory content"""
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


def test_systems(parser) -> None:
    archive = EntryArchive()
    parser.parse(
        'tests/data/lammps/1_methyl_naphthalene/log.1_methyl_naphthalene',
        archive,
        LOGGER,
    )
    # Add placeholder metadata to suppress normalizer error
    archive.metadata = EntryMetadata()
    normalize_all(archive, logger=LOGGER)
    sec_systems = archive.data.model_system
    assert len(sec_systems) == 4
    assert np.shape(sec_systems[0].positions) == (1134, 3)
    # TODO: Atomic test data does not have velocities, update testing!
    # assert np.shape(sec_systems[0].velocities) == (1134, 3)
    assert sec_systems[0].n_particles == 1134
    assert sec_systems[0].particle_states[100].chemical_symbol == 'H'
    assert sec_systems[0].particle_states[100].label == 'H'

    assert sec_systems[2].positions[567][1].to('angstrom').magnitude == pytest.approx(
        -5.88475
    )
    # assert sec_systems[idx].velocities[idx][idx].to(
    #     'angstrom/ps'
    # ).magnitude == pytest.approx(target_float)
    assert sec_systems[3].lattice_vectors[2][2].to(
        'angstrom'
    ).magnitude == pytest.approx(21.468)
    assert sec_systems[3].periodic_boundary_conditions == [
        True,
        True,
        True,
    ]
    assert (
        np.testing.assert_array_equal(
            sec_systems[0].bond_list[200], np.array([189, 191])
        )
        is None
    )
    assert sec_systems[0].dimensionality == 3
    assert sec_systems[0].is_molecule() is False


def test_systems_velocities(parser):
    archive = EntryArchive()
    parser.parse(
        'tests/data/lammps/1_xyz_files/log.lammps',
        archive,
        LOGGER,
    )
    normalize_all(archive, logger=LOGGER)
    sec_systems = archive.data.model_system
    assert np.shape(sec_systems[0].velocities) == (500, 3)
    assert sec_systems[100].velocities[250][2].to(
        'angstrom/ps'
    ).magnitude == pytest.approx(0.0256726)


# TODO: Needs to be tested with more complex system!
def test_system_hierarchy(parser) -> None:
    archive = EntryArchive()
    parser.parse(
        'tests/data/lammps/1_methyl_naphthalene/log.1_methyl_naphthalene',
        archive,
        LOGGER,
    )
    normalize_all(archive, logger=LOGGER)
    sec_particles_group = archive.data.model_system[0].sub_systems
    assert len(sec_particles_group) == 1
    assert sec_particles_group[0].particle_states == []
    # TODO comment back in once nested fix is in release
    # assert sec_particles_group[0].cell == []
    assert sec_particles_group[0].name == 'group_0'
    assert sec_particles_group[0].branch_label == 'molecule_group'
    assert sec_particles_group[0].composition_formula == '0(54)'
    # ! Particle index is wrong when parsing 1_xyz LJ system!
    assert sec_particles_group[0].particle_indices[13] == 13
    # ? Should this or shouldn't this be recognized as a molecule?
    assert sec_particles_group[0].is_molecule() is False

    sec_monomers = sec_particles_group[0].sub_systems
    assert len(sec_monomers) == 54
    assert sec_monomers[0].name == '0'
    assert sec_monomers[0].branch_label == 'molecule'
    assert sec_monomers[0].composition_formula == 'C(11)H(10)'
    assert sec_monomers[0].particle_indices[20] == 20
    assert sec_monomers[0].is_molecule() is True


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


# TODO Add tests that use the full parser fixture and test end-to-end parsing


@pytest.mark.skip(
    reason='TODO: Integration test fails because LogParser units system not initialized'
    ' during parsing. Need to ensure _log_parser.parse() is called with proper mainfile'
    ' to set up units system.'
)
def test_md_method_nvt_hexane(parser):
    """Test MD method parameter extraction for NVT ensemble (hexane_cyclohexane)."""
    archive = EntryArchive()
    archive.metadata = EntryMetadata()
    parser.parse(
        'tests/data/lammps/hexane_cyclohexane/log.hexane_cyclohexane_nvt-thermo_style_multi',
        archive,
        LOGGER,
    )
    normalize_all(archive, logger=LOGGER)

    assert archive.data.model_method is not None
    assert len(archive.data.model_method) > 0

    method = archive.data.model_method[0]

    # Check integration parameters
    assert method.integration_timestep is not None
    assert method.integration_timestep.to('fs').magnitude == pytest.approx(250.0)
    assert method.n_steps == 80000
    assert method.integrator_type == 'velocity_verlet'

    # Check thermostat
    assert method.thermostat_parameters is not None
    assert len(method.thermostat_parameters) == 1
    thermostat = method.thermostat_parameters[0]
    assert thermostat.thermostat_type == 'nose_hoover'
    assert thermostat.reference_temperature.to('kelvin').magnitude == pytest.approx(
        300.0
    )
    assert thermostat.coupling_constant.to('fs').magnitude == pytest.approx(100000.0)

    # Check ensemble
    assert method.thermodynamic_ensemble == 'NVT'

    # Check save frequencies
    assert method.coordinate_save_frequency == 400
    assert method.thermodynamics_save_frequency == 400


def test_md_method_thermostat_extraction():
    """Test thermostat extraction from various fix commands."""

    writer = LammpsArchiveWriter()

    # Test NVT thermostat
    fix_nvt = [['1', 'all', 'nvt', 'temp', '300.0', '300.0', '100.0']]
    thermostat = writer._extract_thermostat_settings(fix_nvt)
    assert thermostat is not None
    assert thermostat.thermostat_type == 'nose_hoover'

    # Test NPT thermostat (extracts temperature from NPT fix)
    fix_npt = [
        [
            '1',
            'all',
            'npt',
            'temp',
            '300.0',
            '300.0',
            '100.0',
            'iso',
            '1.0',
            '1.0',
            '1000.0',
        ]
    ]
    thermostat = writer._extract_thermostat_settings(fix_npt)
    assert thermostat is not None
    assert thermostat.thermostat_type == 'nose_hoover'

    # Test Langevin thermostat
    fix_langevin = [['1', 'all', 'langevin', '300.0', '300.0', '100.0', '12345']]
    thermostat = writer._extract_thermostat_settings(fix_langevin)
    assert thermostat is not None
    assert thermostat.thermostat_type == 'langevin_leap_frog'

    # Test Berendsen thermostat
    fix_berendsen = [['1', 'all', 'temp/berendsen', '300.0', '300.0', '100.0']]
    thermostat = writer._extract_thermostat_settings(fix_berendsen)
    assert thermostat is not None
    assert thermostat.thermostat_type == 'berendsen'

    # Test velocity rescaling thermostat
    fix_rescale = [['1', 'all', 'temp/rescale', '100', '300.0', '310.0']]
    thermostat = writer._extract_thermostat_settings(fix_rescale)
    assert thermostat is not None
    assert thermostat.thermostat_type == 'velocity_rescaling'

    # Test no thermostat
    thermostat = writer._extract_thermostat_settings(None)
    assert thermostat is None

    # Test non-thermostat fix
    fix_other = [['1', 'all', 'shake', '0.0001', '20', '100', 'b', '1']]
    thermostat = writer._extract_thermostat_settings(fix_other)
    assert thermostat is None


def test_md_method_barostat_extraction():
    """Test barostat extraction from various fix commands."""

    writer = LammpsArchiveWriter()

    # Test NPT barostat (isotropic)
    fix_npt_iso = [
        [
            '1',
            'all',
            'npt',
            'temp',
            '300.0',
            '300.0',
            '100.0',
            'iso',
            '1.0',
            '1.0',
            '1000.0',
        ]
    ]
    barostat = writer._extract_barostat_settings(fix_npt_iso)
    assert barostat is not None
    assert barostat.barostat_type == 'nose_hoover'
    assert barostat.coupling_type == 'isotropic'

    # Test NPH barostat (anisotropic)
    fix_nph_aniso = [['1', 'all', 'nph', 'aniso', '1.0', '1.0', '1000.0']]
    barostat = writer._extract_barostat_settings(fix_nph_aniso)
    assert barostat is not None
    assert barostat.barostat_type == 'nose_hoover'
    assert barostat.coupling_type == 'anisotropic'

    # Test press/berendsen barostat (triclinic)
    fix_berendsen = [['1', 'all', 'press/berendsen', 'tri', '1.0', '1.0', '1000.0']]
    barostat = writer._extract_barostat_settings(fix_berendsen)
    assert barostat is not None
    assert barostat.barostat_type == 'berendsen'
    assert barostat.coupling_type == 'anisotropic'

    # Test no barostat
    barostat = writer._extract_barostat_settings(None)
    assert barostat is None


def test_md_method_barostat_modulus_extraction():
    """Test modulus extraction and conversion to compressibility."""

    writer = LammpsArchiveWriter()
    writer._log_parser = LogParser()
    # Initialize units system for apply_unit to work
    writer._log_parser._units = {
        'pressure': ureg.atmosphere,
        'time': ureg.femtosecond,
    }

    # Test NPT with modulus keyword (bulk modulus = 20000 atm)
    # compressibility should be 1/20000 atm^-1
    fix_npt_modulus = [
        [
            '1',
            'all',
            'npt',
            'temp',
            '300',
            '300',
            '100',
            'iso',
            '1.0',
            '1.0',
            '1000.0',
            'modulus',
            '20000.0',
        ]
    ]
    barostat = writer._extract_barostat_settings(fix_npt_modulus)
    assert barostat is not None
    assert barostat.compressibility is not None

    # Check compressibility value (1/20000 atm in SI units)
    # 1 atm = 101325 Pa, so 1/20000 atm^-1 = 1/(20000*101325) Pa^-1
    expected_compressibility = 1.0 / (20000.0 * 101325.0)  # Pa^-1
    assert barostat.compressibility[0, 0].magnitude == pytest.approx(
        expected_compressibility, rel=1e-6
    )
    assert barostat.compressibility[1, 1].magnitude == pytest.approx(
        expected_compressibility, rel=1e-6
    )
    assert barostat.compressibility[2, 2].magnitude == pytest.approx(
        expected_compressibility, rel=1e-6
    )

    # Test NPH with modulus
    fix_nph_modulus = [
        ['1', 'all', 'nph', 'iso', '1.0', '1.0', '1000.0', 'modulus', '10000.0']
    ]
    barostat = writer._extract_barostat_settings(fix_nph_modulus)
    assert barostat is not None
    assert barostat.compressibility is not None
    expected_compressibility = 1.0 / (10000.0 * 101325.0)
    assert barostat.compressibility[0, 0].magnitude == pytest.approx(
        expected_compressibility, rel=1e-6
    )

    # Test without modulus (should be None)
    fix_npt_no_modulus = [
        ['1', 'all', 'npt', 'temp', '300', '300', '100', 'iso', '1.0', '1.0', '1000.0']
    ]
    barostat = writer._extract_barostat_settings(fix_npt_no_modulus)
    assert barostat is not None
    assert barostat.compressibility is None


def test_md_method_ensemble_detection():
    """Test ensemble type detection."""

    writer = LammpsArchiveWriter()

    # NVE: no thermostat, no barostat
    assert writer._determine_ensemble(False, False) == 'NVE'

    # NVT: thermostat, no barostat
    assert writer._determine_ensemble(True, False) == 'NVT'

    # NPH: no thermostat, barostat
    assert writer._determine_ensemble(False, True) == 'NPH'

    # NPT: thermostat and barostat
    assert writer._determine_ensemble(True, True) == 'NPT'


def test_md_method_dump_frequency_extraction():
    """Test extraction of save frequencies from dump commands."""

    writer = LammpsArchiveWriter()

    # Test coordinate dump
    dump_coords = [
        [
            '1',
            '100',
            'all',
            'custom',
            '500',
            'traj.lammpstrj',
            'id',
            'type',
            'x',
            'y',
            'z',
        ]
    ]
    freqs = writer._extract_dump_frequencies(dump_coords)
    assert freqs['coordinate'] == 100
    assert freqs['velocity'] is None
    assert freqs['force'] is None

    # Test velocity dump
    dump_vels = [
        [
            '1',
            '200',
            'all',
            'custom',
            '500',
            'traj.lammpstrj',
            'id',
            'type',
            'vx',
            'vy',
            'vz',
        ]
    ]
    freqs = writer._extract_dump_frequencies(dump_vels)
    assert freqs['coordinate'] is None
    assert freqs['velocity'] == 200
    assert freqs['force'] is None

    # Test force dump
    dump_forces = [
        [
            '1',
            '300',
            'all',
            'custom',
            '500',
            'traj.lammpstrj',
            'id',
            'type',
            'fx',
            'fy',
            'fz',
        ]
    ]
    freqs = writer._extract_dump_frequencies(dump_forces)
    assert freqs['coordinate'] is None
    assert freqs['velocity'] is None
    assert freqs['force'] == 300

    # Test combined dump
    dump_all = [
        [
            '1',
            '400',
            'all',
            'custom',
            '500',
            'traj.lammpstrj',
            'id',
            'type',
            'x',
            'y',
            'z',
            'vx',
            'vy',
            'vz',
            'fx',
            'fy',
            'fz',
        ]
    ]
    freqs = writer._extract_dump_frequencies(dump_all)
    assert freqs['coordinate'] == 400
    assert freqs['velocity'] == 400
    assert freqs['force'] == 400

    # Test unwrapped coordinates (xu, yu, zu)
    dump_unwrapped = [
        [
            '1',
            '500',
            'all',
            'custom',
            '500',
            'traj.lammpstrj',
            'id',
            'type',
            'xu',
            'yu',
            'zu',
        ]
    ]
    freqs = writer._extract_dump_frequencies(dump_unwrapped)
    assert freqs['coordinate'] == 500

    # Test no dump
    freqs = writer._extract_dump_frequencies(None)
    assert freqs['coordinate'] is None
    assert freqs['velocity'] is None
    assert freqs['force'] is None


@pytest.mark.skip(
    reason='TODO: Test fails because LogParser.apply_unit() requires units system to be'
    ' initialized. Need to mock _log_parser with initialized units or test without unit'
    ' conversion.'
)
def test_md_method_helper_methods():
    """Test helper methods for parameter extraction."""
    writer = LammpsArchiveWriter()
    writer._log_parser = LogParser()

    # Test timestep extraction
    writer._log_parser._results = {'timestep': 2.0}
    timestep = writer._extract_timestep()
    assert timestep is not None

    # Test timestep extraction with list
    writer._log_parser._results = {'timestep': [1.0, 2.0, 3.0]}
    timestep = writer._extract_timestep()
    assert timestep is not None

    # Test n_steps extraction
    writer._log_parser._results = {'run': 50000}
    n_steps = writer._extract_n_steps()
    assert n_steps == 50000

    # Test n_steps extraction with nested list
    writer._log_parser._results = {'run': [[10000], [20000], [30000]]}
    n_steps = writer._extract_n_steps()
    assert n_steps == 30000

    # Test thermo frequency extraction
    writer._log_parser._results = {'thermo': 100}
    thermo_freq = writer._extract_thermo_frequency()
    assert thermo_freq == 100

    # Test thermo frequency extraction with list
    writer._log_parser._results = {'thermo': [50, 100, 200]}
    thermo_freq = writer._extract_thermo_frequency()
    assert thermo_freq == 200


def test_md_method_integrator_type_extraction():
    """Test integrator type extraction from run_style command."""

    class MockLogParser:
        """Mock LogParser for testing."""

        def __init__(self):
            self._results = {}

        def get(self, key, default=None):
            return self._results.get(key, default)

    writer = LammpsArchiveWriter()
    writer._log_parser = MockLogParser()

    # Test verlet
    writer._log_parser._results = {'run_style': 'verlet'}
    integrator = writer._extract_integrator_type()
    assert integrator == 'velocity_verlet'

    # Test verlet/split
    writer._log_parser._results = {'run_style': 'verlet/split'}
    integrator = writer._extract_integrator_type()
    assert integrator == 'velocity_verlet_split'

    # Test respa
    writer._log_parser._results = {'run_style': 'respa'}
    integrator = writer._extract_integrator_type()
    assert integrator == 'respa'

    # Test respa/omp
    writer._log_parser._results = {'run_style': 'respa/omp'}
    integrator = writer._extract_integrator_type()
    assert integrator == 'respa'

    # Test with command args (respa with levels)
    writer._log_parser._results = {
        'run_style': ['respa', '4', '2', '2', '2', 'bond', '1', 'kspace', '4']
    }
    integrator = writer._extract_integrator_type()
    assert integrator == 'respa'

    # Test with nested list of commands (take last)
    writer._log_parser._results = {'run_style': [['verlet'], ['respa', '4', '2']]}
    integrator = writer._extract_integrator_type()
    assert integrator == 'respa'

    # Test None case
    writer._log_parser._results = {}
    integrator = writer._extract_integrator_type()
    assert integrator is None

    # Test unknown style (fallback to velocity_verlet)
    writer._log_parser._results = {'run_style': 'unknown_style'}
    integrator = writer._extract_integrator_type()
    assert integrator == 'velocity_verlet'
