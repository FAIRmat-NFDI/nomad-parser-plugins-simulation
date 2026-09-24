import os
import re
import tempfile
from io import BytesIO, StringIO

import numpy as np
import pytest
from nomad.utils import get_logger

from nomad_simulation_parsers.parsers.lammps.file_parsers import DataParser, LogParser
from nomad_simulation_parsers.parsers.lammps.trajectory_parsers import (
    TrajParser,
    TrajParsers,
    XYZTrajParser,
)
from tests.parsers.common import approx, assert_approx

LOGGER = get_logger(__name__)


@pytest.mark.unit
class TestLogParser:
    @pytest.fixture
    def parser(self):
        parser = LogParser()
        parser.mainfile = __file__
        return parser

    def test_maps_boundary_to_periodic_flags(self, parser):
        parser._results = {'boundary': ['p', 'f', 'p']}

        assert parser.get_pbc() == [True, False, True]

    def test_maps_sampling_and_thermostat_settings(self, parser):
        parser._results = {
                'units': ['real'],
                'fix': [['1', 'all', 'nvt', 'temp', '300', '300', '100']],
            }
        

        assert parser.get_sampling_method() == ('molecular_dynamics', 'nvt')
        settings = parser.get_thermostat_settings()
        assert settings['target_T'].to('K').magnitude == approx(300)
        assert settings['thermostat_tau'].to('fs').magnitude == approx(100)


    def test_reads_workflow(self, parser):
        parser.mainfile = (
            'tests/data/lammps/hexane_cyclohexane/'
            'log.hexane_cyclohexane_nvt'
        )
        parser.parse('units')
        parser.parse('fix')

        assert parser.get_sampling_method() == (
            'molecular_dynamics',
            'nvt',
        )
        thermostat = parser.get_thermostat_settings()
        assert thermostat['target_T'].to('K').magnitude == approx(300.0)
        assert thermostat['thermostat_tau'].to('fs').magnitude == approx(100.0)

    def test_converts_thermodynamic_data_from_mainfile(self, parser):
        source = 'tests/data/lammps/hexane_cyclohexane/log.hexane_cyclohexane_nvt'
        parser.mainfile = source
        parser.parse('units')
        units = parser.units

        parser.mainfile = f'{source}-thermo_style_multi'
        parser.parse('thermo_data')
        # The auxiliary thermo log does not repeat the units command. Reuse
        # the units declared by the corresponding main log.
        parser._units = units

        thermo_data = parser.get_thermodynamic_data()

        assert len(thermo_data['Step']) == 201
        assert thermo_data['Temp'][103].to('K').magnitude == approx(291.4591)
        assert thermo_data['Press'][56].to('Pa').magnitude == approx(
            -77642135.4975
        )
        assert thermo_data['TotEng'][21].to('J').magnitude == approx(
            8.866891968814006e-18
        )
        assert thermo_data['KinEng'][21].to('J').magnitude == approx(
            3.900969749768176e-18
        )
        assert thermo_data['PotEng'][21].to('J').magnitude == approx(
            4.96592221904583e-18
        )
        assert thermo_data['E_bond'][21].to('J').magnitude == approx(
            1.917461089369534e-18
        )
        assert thermo_data['E_angle'][21].to('J').magnitude == approx(
            1.7659645080599946e-18
        )

    def test_converts_provided_thermodynamic_data(self, parser):
        parser._results = {
            'units': ['real'],
            'thermo_data': {
                'Step': [0, 10],
                'Temp': [300.0, 310.0],
                'TotEng': [-5.0, -4.0],
                'Press': [1.0, 2.0],
            },
        }

        thermo_data = parser.get_thermodynamic_data()

        assert thermo_data['Step'] == [0, 10]
        assert_approx(
            thermo_data['Temp'].to('K').magnitude, [300.0, 310.0]
        )
        assert_approx(
            thermo_data['TotEng'].to('J').magnitude,
            [-3.47384767e-20, -2.77907814e-20],
            rtol=1e-7,
        )
        assert_approx(
            thermo_data['Press'].to('atm').magnitude, [1.0, 2.0]
        )

    @pytest.mark.parametrize(
        ('mainfile', 'traj_files', 'expected_match'),
        [
            (
                'log.hexane_nvt',
                ['hexane_nvt.lammpstrj', 'water.xyz'],
                'hexane_nvt.lammpstrj',
            ),
            (
                'log.hexane_nvt',
                ['water.lammpstrj', 'hex_nvt_output.xyz'],
                'hex_nvt_output.xyz',
            ),
            (
                'log.1_methyl_naphthalene',
                ['naph_298_eq.lammpstrj', 'water.xyz'],
                'naph_298_eq.lammpstrj',
            ),
            (
                'log.hexane_nvt',
                ['hexanol_simulation.lammpstrj', 'completely_different.xyz'],
                'hexanol_simulation.lammpstrj',
            ),
            (
                'log.hexane_nvt',
                ['water.lammpstrj', 'polymer.xyz', 'argon.dcd'],
                'water.lammpstrj',
            ),
            (
                'log.hexane_nvt',
                ['HEXANE_NVT.lammpstrj', 'other.xyz'],
                'HEXANE_NVT.lammpstrj',
            ),
        ],
    )
    def test_find_best_matching_file(
        self, parser, mainfile, traj_files, expected_match
    ):
        """Test trajectory file matching with various scenarios"""
        parser.logger = LOGGER

        result = parser.find_best_matching_file(traj_files, mainfile)

        assert len(result) == 1
        assert result[0] == expected_match

    @pytest.mark.parametrize(
        ('mainfile', 'aux_files', 'dump', 'expected_file'),
        [
            (
                'tests/data/lammps/1_xyz_files/log.test',
                None,
                [
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
                ],
                'pos_vel.xyz',
            ),
            (
                'log.hexane_nvt',
                ['hexane_nvt.lammpstrj', 'water.xyz', 'readme.txt'],
                None,
                'hexane_nvt.lammpstrj',
            ),
            (
                'tests/data/lammps/2_xyz_files/log.lammps',
                None,
                [
                    ['id1', 'all', 'custom', '100', 'traj.lammpstrj', 'x', 'y', 'z'],
                    ['id2', 'all', 'custom', '200', 'traj.lammpstrj', 'x', 'y', 'z'],
                ],
                'traj.lammpstrj',
            ),
            (
                'log.lammps',
                ['trajectory.lammpstrj'],
                [
                    [
                        '1',
                        'all',
                        'custom',
                        '100',
                        'trajectory.lammpstrj',
                        'id',
                        'x',
                        'y',
                        'z',
                    ]
                ],
                'trajectory.lammpstrj',
            ),
        ],
    )
    def test_get_traj_files(  # noqa: PLR0913, PLR0917
        self, parser, tmp_dir, mainfile, aux_files, dump, expected_file
    ):
        if aux_files:
            for filename in aux_files:
                open(os.path.join(tmp_dir, filename), 'w').close()
            mainfile = os.path.join(tmp_dir, mainfile)
            open(mainfile, 'w').close()

        parser.mainfile = mainfile
        parser.logger = LOGGER
        parser._results = {'dump': dump}
        traj_files = parser.get_traj_files()

        assert len(traj_files) == 1
        assert os.path.basename(traj_files[0]) == expected_file
        assert all(os.path.isabs(filename) for filename in traj_files)

    @pytest.mark.parametrize(
        ('read_data', 'expected_files'),
        [([], []), (['structure.data'], ['structure.data'])],
    )
    def test_get_data_files(
        self, parser, tmp_dir, read_data, expected_files
    ):
        mainfile = os.path.join(tmp_dir, 'log.test')
        open(mainfile, 'w').close()
        for filename in read_data:
            open(os.path.join(tmp_dir, filename), 'w').close()

        parser.mainfile = mainfile
        parser._results = {'read_data': read_data}
        parser.logger = LOGGER

        data_files = parser.get_data_files()

        assert [os.path.basename(filename) for filename in data_files] == expected_files

    @pytest.mark.parametrize(
        ('content', 'pattern', 'expected', 'binary'),
        [
            (
                'LAMMPS data file via write_data, version 12 Dec 2018\n'
                '\n1000 atoms\n',
                'LAMMPS data file',
                True,
                False,
            ),
            ('LAMMPS Description\n\n500 atoms\n', 'LAMMPS Description', True, False),
            ('Random text file\nNot a LAMMPS file.', 'LAMMPS data file', False, False),
            (
                'Random text file\nNot a LAMMPS file.',
                'LAMMPS Description',
                False,
                False,
            ),
            (
                b'LAMMPS data file\n1000 atoms\n' + b'\x00\xff\xfe\xfd' * 100,
                'LAMMPS data file',
                True,
                True,
            ),
            (b'\x00\x01\x02\x03' * 256, 'LAMMPS data file', False, True),
            (b'\x00\x01\x02\x03' * 256, 'LAMMPS Description', False, True),
        ],
    )
    def test_data_file_header_matching(self, content, pattern, expected, binary):
        file_obj = BytesIO(content) if binary else StringIO(content)
        file_obj.seek(0)
        header = file_obj.read(1024)
        header = header if isinstance(header, str) else header.decode(errors='ignore')

        assert (re.search(pattern, header) is not None) is expected

    def test_data_file_header_priority(
        self,
    ):
        """Test header pattern matching priority.

        Covers the standard and alternative LAMMPS data-file headers.
        """

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
        assert (
            check_file_header_in_memory(standard_file, 'LAMMPS data file') is not None
        )

        # File with alternative header should only match alternative pattern
        alt_file = StringIO("""LAMMPS Description
    500 atoms
    """)
        assert check_file_header_in_memory(alt_file, 'LAMMPS data file') is None
        assert check_file_header_in_memory(alt_file, 'LAMMPS Description') is not None

    def test_parses_commands_and_units(self, parser, tmp_path):
        log = tmp_path / 'log.lammps'
        log.write_text(
            'LAMMPS (12 Dec 2018)\n'
            'units real\n'
            'boundary p p f\n'
            'dump 1 all custom 100 trajectory.lammpstrj id type x y z\n'
            'read_data structure.data\n'
        )

        parser.mainfile = str(log)
        parser.parse()

        assert parser.get('units') == ['real']
        assert parser.get('boundary') == [['p', 'p', 'f']]
        assert parser.get('dump')[0][4] == 'trajectory.lammpstrj'
        assert parser.get('read_data') == ['structure.data']
        assert str(parser.units['distance']) == 'angstrom'

    def test_extracts_interaction_styles_and_coefficients(self, parser):
        parser._results = {
                'pair_style': [['lj/cut', '2.5']],
                'pair_coeff': [[1, 1, 1.0, 1.0, 2.5]],
            }
        

        assert parser.get_interactions() == [('lj/cut 2.5', [[1, 1, 1.0, 1.0, 2.5]])]

    def test_reset_clears_cached_results(self, parser, tmp_path):
        log = tmp_path / 'log.lammps'
        log.write_text('LAMMPS (12 Dec 2018)\nunits real\n')
        parser.mainfile = str(log)
        parser.parse()
        assert parser.get('units') == ['real']

        parser.reset()

        assert parser._units is None


@pytest.mark.unit
class TestDataParser:
    def test_reads_headers_sections_and_interactions(self, tmp_path):
        data_file = tmp_path / 'data.test'
        data_file.write_text(
            'LAMMPS data file via test\n\n'
            '2 atoms\n1 bonds\n2 atom types\n1 bond types\n\n'
            'Masses\n\n'
            '1 12.011\n2 1.008\n\n'
            'Atoms # full\n\n'
            '1 1 1 -0.1 0.0 0.0 0.0\n'
            '2 1 2 0.1 1.0 0.0 0.0\n\n'
            'Bonds\n\n'
            '1 1 1 2\n\n'
            'Bond Coeffs\n\n'
            '1 300.0 1.0\n'
        )

        parser = DataParser()
        parser.mainfile = str(data_file)
        parser.init_quantities()
        parser.parse()

        assert parser.get('atoms') == [2]
        assert parser.get('bonds') == [1]
        assert_approx(
            parser.get('Masses')[0][1], [[1, 12.011], [2, 1.008]]
        )
        assert parser.get('Atoms')[0][0] == 'full'
        assert_approx(parser.get('Bonds')[0][1], [[1, 1, 1, 2]])
        assert_approx(parser.get_interactions()[0][1], [[1, 300.0, 1.0]])


@pytest.mark.unit
class TestTrajParser:
    @pytest.fixture
    def parser(self):
        return TrajParser()

    def test_extracts_frames_and_vectors(self, parser, tmp_path):
        trajectory = tmp_path / 'trajectory.lammpstrj'
        trajectory.write_text(
            'ITEM: TIMESTEP\n0\n'
            'ITEM: NUMBER OF ATOMS\n2\n'
            'ITEM: BOX BOUNDS pp pp pp\n0 10\n0 20\n0 30\n'
            'ITEM: ATOMS id type x y z vx vy vz fx fy fz\n'
            '2 1 2 3 4 0.2 0.3 0.4 2 3 4\n'
            '1 1 1 2 3 0.1 0.2 0.3 1 2 3\n'
        )

        parser.mainfile = str(trajectory)
        parser.init_quantities()

        assert parser.get_step(0) == 0
        assert parser.get_n_atoms(0) == 2
        assert_approx(parser.get_positions(0), [[1, 2, 3], [2, 3, 4]])
        assert_approx(
            parser.get_velocities(0), [[0.1, 0.2, 0.3], [0.2, 0.3, 0.4]]
        )
        assert_approx(parser.get_forces(0), [[1, 2, 3], [2, 3, 4]])
        assert_approx(parser.get_lattice_vectors(0), np.diag([10, 20, 30]))
        assert parser.get_pbc(0) == [True, True, True]

    def test_extracts_atom_labels_and_trajectory_properties(self, parser, tmp_path):
        trajectory = tmp_path / 'trajectory.lammpstrj'
        trajectory.write_text(
            'ITEM: TIMESTEP\n0\n'
            'ITEM: NUMBER OF ATOMS\n2\n'
            'ITEM: BOX BOUNDS pp ff pp\n0 10\n0 20\n0 30\n'
            'ITEM: ATOMS id type x y z\n'
            '1 1 0 1 2\n'
            '2 2 3 4 5\n'
        )

        parser.mainfile = str(trajectory)
        parser.init_quantities()
        parser.parse()
        parser.masses = np.array([[1, 12.011], [2, 1.008]])

        assert parser.with_trajectory is True
        assert parser.n_frames == 1
        assert parser.get_atom_labels(0) == ['C', 'H']
        assert parser.get_pbc(0) == [True, False, True]
        pbc, cell = parser.get_pbc_cell('pp ff pp 0 10 0 20 0 30')
        assert pbc == [True, False, True]
        assert_approx(cell, np.diag([10, 20, 30]))

    def test_reset_clears_trajectory_state(self, parser):
        parser.masses = np.array([[1, 12.011]])
        parser._chemical_symbols = {1: 'C'}

        parser.reset()

        assert parser.masses is None
        assert parser.chemical_symbols is None

    def test_get_lattice_vectors(self, parser):
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

        with tempfile.NamedTemporaryFile(
            mode='w', suffix='.lammpstrj', delete=False
        ) as f:
            f.write(traj_content)
            temp_file = f.name

        try:
            parser.mainfile = temp_file
            parser.logger = LOGGER
            parser.init_quantities()

            # Test lattice vectors
            lattice_vectors = parser.get_lattice_vectors(0)
            assert lattice_vectors is not None
            assert lattice_vectors.shape == (3, 3)
            assert lattice_vectors[0, 0] == approx(20.0)  # 10.0 - (-10.0)
            assert lattice_vectors[1, 1] == approx(30.0)  # 15.0 - (-15.0)
            assert lattice_vectors[2, 2] == approx(24.0)  # 12.0 - (-12.0)

        finally:
            os.unlink(temp_file)

    def test_get_n_atoms(self, parser):
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

        with tempfile.NamedTemporaryFile(
            mode='w', suffix='.lammpstrj', delete=False
        ) as f:
            f.write(traj_content)
            temp_file = f.name

        try:
            parser.mainfile = temp_file
            parser.logger = LOGGER
            parser.init_quantities()

            # Test n_atoms
            n_atoms = parser.get_n_atoms(0)
            assert n_atoms == 5

        finally:
            os.unlink(temp_file)

    def test_get_step(self, parser):
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

        with tempfile.NamedTemporaryFile(
            mode='w', suffix='.lammpstrj', delete=False
        ) as f:
            f.write(traj_content)
            temp_file = f.name

        try:
            parser.mainfile = temp_file
            parser.logger = LOGGER
            parser.init_quantities()

            # Test timestep for first frame
            step_0 = parser.get_step(0)
            assert step_0 == 100

            # Test timestep for second frame
            step_1 = parser.get_step(1)
            assert step_1 == 200

        finally:
            os.unlink(temp_file)

    def test_none_returns(self, parser):
        """Test that methods return None when data is missing"""
        # Create parser with minimal/no data
        traj_content = """ITEM: TIMESTEP
    0
    """

        with tempfile.NamedTemporaryFile(
            mode='w', suffix='.lammpstrj', delete=False
        ) as f:
            f.write(traj_content)
            temp_file = f.name

        try:
            parser.mainfile = temp_file
            parser.logger = LOGGER
            parser.init_quantities()

            # All should return None when data is missing
            assert parser.get_lattice_vectors(0) is None
            assert parser.get_pbc(0) is None
            # n_atoms falls back to counting positions, which also returns None
            assert parser.get_n_atoms(0) is None

        finally:
            os.unlink(temp_file)

    def test_with_real_data(self, parser):
        """Test methods with existing real test data"""
        parser.mainfile = 'tests/data/lammps/1_xyz_files/pos_vel.xyz'
        parser.logger = LOGGER
        parser.init_quantities()

        # Test that methods work with real data
        n_atoms = parser.get_n_atoms(5)
        assert n_atoms is not None
        assert isinstance(n_atoms, int)
        assert n_atoms > 0

        step = parser.get_step(3)
        assert step is not None
        assert isinstance(step, int)

        # If the real data has PBC info, test it
        pbc = parser.get_pbc(1)
        lattice_vectors = parser.get_lattice_vectors(1)
        # These might be None if the file doesn't include PBC info
        assert isinstance(pbc, list)
        assert isinstance(lattice_vectors, np.ndarray)

    def test_unwrapped_pos(self, parser):
        # 1_xyz dataset (CG), file type 'custom' -> TrajParser
        parser.mainfile = 'tests/data/lammps/1_xyz_files/pos_vel.xyz'
        parser.init_quantities()
        # TODO: add assertion for calculation
        positions = parser.get_positions(1)
        assert positions[452][2] == approx(5.99898)
        velocities = parser.get_velocities(2)
        assert velocities[457][-2] == approx(-0.928553)



@pytest.mark.unit
class TestXYZTrajParser:
    def test_traj_xyz(self):
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

            assert xyz_parser.n_frames == 2
            positions = xyz_parser.get_positions(1)
            assert positions[2][1] == approx(-0.845205)

        finally:
            os.unlink(temp_file)


class StubTrajectoryParser:
    def __init__(self, n_frames=None, n_atoms=None):
        self.logger = None
        self.n_frames = n_frames
        self._n_atoms = n_atoms

    def parse(self):
        return None

    def get_n_atoms(self, index):
        return self._n_atoms


@pytest.mark.unit
class TestTrajParsers:
    @pytest.mark.parametrize(
        ('description', 'content', 'expected_pbc', 'expected_cell'),
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
                np.diag([26.7069, 28.8056, 24.8952]),
            )
        ],
    )
    def test_pbc_cell_extraction(
        self, description, content, expected_pbc, expected_cell
    ):
        """Test PBC and cell extraction from synthetic LAMMPS trajectory content"""
        parser = TrajParser()
        parser._mainfile_contents = 'dummy'
        parser._file_handler = content.encode('utf-8')
        parser.init_quantities()

        parsers = TrajParsers([parser])
        pbc_cell = parsers.eval('pbc_cell')
        assert len(pbc_cell) == 1, f'{description} - pbc_cell not extracted'

        pbc, cell = pbc_cell[0]
        assert pbc == expected_pbc, f'{description} - wrong PBC'
        assert cell == approx(expected_cell), f'{description} - wrong cell'

    def test_eval_uses_parser_with_available_quantity(self):
        parsers = TrajParsers(
            [
                StubTrajectoryParser(n_frames=None, n_atoms=None),
                StubTrajectoryParser(n_frames=3, n_atoms=12),
            ]
        )

        assert parsers.eval('n_frames') == 3
        assert parsers.eval('get_n_atoms', 0) == 12

    def test_indexes_configured_parsers(self):
        first = StubTrajectoryParser(n_frames=1)
        second = StubTrajectoryParser(n_frames=2)
        parsers = TrajParsers([first, second])

        assert parsers[0] is first
        assert parsers[1] is second

    def test_dispatches_xyz_trajectory(self, tmp_path):
        trajectory = tmp_path / 'trajectory.xyz'
        trajectory.write_text(
            '5\nAtoms. Timestep: 0\n'
            '1 4.39861 0.0809956 -1.6196\n'
            '2 3.65138 0.778109 -1.97822\n'
            '2 4.72189 -0.655793 -2.40238\n'
            '2 5.23117 0.689443 -1.27747\n'
            '2 3.94587 -0.457468 -0.7756\n'
            '5\nAtoms. Timestep: 400\n'
            '1 4.17634 0.0441698 -1.4592\n'
            '2 3.33775 0.267888 -2.08495\n'
            '2 4.74748 -0.845205 -1.78471\n'
            '2 4.8507 0.87915 -1.4652\n'
            '2 3.77509 -0.143827 -0.474483\n'
        )
        parser = XYZTrajParser()
        parser.mainfile = str(trajectory)
        parser.init_quantities()
        parsers = TrajParsers([parser])

        assert parsers.eval('n_frames') == 2
        assert parser.get_positions(1)[2][1] == approx(-0.845205)


