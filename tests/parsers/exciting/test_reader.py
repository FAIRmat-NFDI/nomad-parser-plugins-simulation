import numpy as np
import pytest
from nomad.units import ureg

from nomad_simulation_parsers.parsers.exciting.eigval_parser import (
    EigvalFileParser,
    str_to_eigenvalues,
)
from nomad_simulation_parsers.parsers.exciting.info_parser import (
    InfoFileParser,
    str_to_array,
    str_to_energy_dict,
)


@pytest.mark.unit
class TestInfoReader:
    def test_extracts_initialization_data_from_minimal_info_out(self, tmp_path):
        mainfile = tmp_path / 'INFO.OUT'
        mainfile.write_text(
            'EXCITING v1.2.3 started\n'
            'Starting initialization\n'
            'Lattice vectors (cartesian) :\n'
            ' 1.0 0.0 0.0\n'
            ' 0.0 1.0 0.0\n'
            ' 0.0 0.0 1.0\n'
            'Spin treatment : spin-polarised\n'
            'Ending initialization\n'
        )

        parser = InfoFileParser()
        parser.mainfile = str(mainfile)
        source = parser.to_dict()

        assert source['program_version'] == 'v1.2.3'
        assert source['initialization']['x_exciting_spin_treatment'] == 'spin-polarised'
        lattice_vectors = source['initialization']['lattice_vectors']
        np.testing.assert_allclose(
            lattice_vectors.to('bohr').magnitude,
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        )

    def test_keeps_header_data_when_initialization_is_truncated(self, tmp_path):
        mainfile = tmp_path / 'INFO.OUT'
        mainfile.write_text(
            'EXCITING v1.2.3 started\n'
            'Starting initialization\n'
            'Lattice vectors (cartesian) :\n'
            ' 1.0 0.0 0.0\n'
        )

        parser = InfoFileParser()
        parser.mainfile = str(mainfile)
        source = parser.to_dict()

        assert source == {'program_version': 'v1.2.3'}

    def test_ignores_malformed_numeric_initialization_value(self, tmp_path):
        mainfile = tmp_path / 'INFO.OUT'
        mainfile.write_text(
            'EXCITING v1.2.3 started\n'
            'Starting initialization\n'
            'Unit cell volume : not-a-number\n'
            'Ending initialization\n'
        )

        parser = InfoFileParser()
        parser.mainfile = str(mainfile)
        source = parser.to_dict()

        assert source['program_version'] == 'v1.2.3'
        assert 'x_exciting_unit_cell_volume' not in source['initialization']

    @pytest.mark.parametrize('newline', ['\n', '\r\n'])
    def test_preserves_initialization_values_across_line_endings(
        self, tmp_path, newline
    ):
        mainfile = tmp_path / 'INFO.OUT'
        mainfile.write_text(
            newline.join(
                [
                    'EXCITING v1.2.3 started',
                    'Starting initialization',
                    'Spin treatment : spin-polarised',
                    'Ending initialization',
                    '',
                ]
            ),
            newline='',
        )

        parser = InfoFileParser()
        parser.mainfile = str(mainfile)
        source = parser.to_dict()

        assert source['initialization']['x_exciting_spin_treatment'] == 'spin-polarised'

    @pytest.mark.parametrize(
        ('line', 'key', 'expected', 'unit'),
        [
            ('Smearing scheme : Gaussian', 'smearing_kind', 'Gaussian', None),
            (
                'Smearing scheme : Extended linear tetrahedron method',
                'smearing_kind',
                ['Extended', 'linear', 'tetrahedron', 'method'],
                None,
            ),
            ('Smearing width : 0.00100000', 'smearing_width', 0.001, 'hartree'),
            ('k-point grid : 6 6 4', 'kpoint_grid', [6, 6, 4], None),
            ('Total number of k-points : 24', 'x_exciting_number_kpoints', 24, None),
            (
                'R^MT_min * |G+k|_max (rgkmax) : 7.00000000',
                'x_exciting_rgkmax',
                7.0,
                None,
            ),
            (
                'Total number of atoms per unit cell : 2',
                'x_exciting_number_of_atoms',
                2,
                None,
            ),
            (
                'Unit cell volume : 76.1170172117',
                'x_exciting_unit_cell_volume',
                76.1170172117,
                'bohr ** 3',
            ),
        ],
        ids=[
            'smearing-kind-single-word',
            'smearing-kind-multi-word',
            'smearing-width',
            'k-point-grid',
            'k-points-reduced',
            'rgkmax',
            'number-of-atoms',
            'unit-cell-volume',
        ],
    )
    def test_extracts_method_and_system_settings(
        self, tmp_path, line, key, expected, unit
    ):
        mainfile = tmp_path / 'INFO.OUT'
        mainfile.write_text(
            'EXCITING SODIUM started\n'
            'Starting initialization\n'
            f' {line}\n'
            'Ending initialization\n'
        )

        parser = InfoFileParser()
        parser.mainfile = str(mainfile)
        value = parser.to_dict()['initialization'][key]

        if unit is not None:
            value = value.to(unit).magnitude
        if isinstance(value, np.ndarray):
            value = value.tolist()
        if isinstance(expected, float):
            assert value == pytest.approx(expected)
        else:
            assert value == expected

    def test_extracts_version_hash_and_species_mass(self, tmp_path):
        mainfile = tmp_path / 'INFO.OUT'
        mainfile.write_text(
            'EXCITING SODIUM started\n'
            'version hash id: c71435148380004635d32fdbcc20adbad59a105\n'
            'Starting initialization\n'
            ' Species :    1 (C)\n'
            '     atomic mass                            :   21894.16673000\n'
            '     atomic positions (lattice) :\n'
            '       1 :   0.00000000  0.00000000  0.00000000\n'
            '\n'
            'Ending initialization\n'
        )

        parser = InfoFileParser()
        parser.mainfile = str(mainfile)
        source = parser.to_dict()

        assert source['hash_id'] == 'c71435148380004635d32fdbcc20adbad59a105'
        species = source['initialization']['species'][0]
        assert species['atomic_mass'].to('electron_mass').magnitude == pytest.approx(
            21894.16673
        )

    def test_reads_charge_distance_in_elementary_charges(self, tmp_path):
        mainfile = tmp_path / 'INFO.OUT'
        mainfile.write_text(
            'EXCITING SODIUM started\n'
            'Self-consistent loop started\n'
            '+ SCF iteration number :    1\n'
            ' Charge distance                   (target) :  0.752460E-02  '
            '( 0.100000E-04)\n'
            '\n'
            '++++++++++++++++++++\n'
            'Groundstate module stopped\n'
        )

        parser = InfoFileParser()
        parser.mainfile = str(mainfile)
        scf_step = parser.to_dict()['groundstate']['scf_iteration'][0]

        np.testing.assert_allclose(
            scf_step['x_exciting_charge_convergence'].to('elementary_charge').magnitude,
            [0.752460e-02, 0.100000e-04],
        )

    def test_reads_array_and_energy_helpers(self):
        np.testing.assert_allclose(str_to_array('1 2 3\n4 5 6'), [[1, 2, 3], [4, 5, 6]])
        energies = str_to_energy_dict('total: -1.5\nkinetic: 0.5')
        assert energies['total'].to('hartree').magnitude == pytest.approx(-1.5)
        assert energies['kinetic'].units == ureg.hartree


@pytest.mark.unit
class TestEigvalReader:
    def test_extracts_k_point_eigenvalues_and_occupancies(self, tmp_path):
        mainfile = tmp_path / 'EIGVAL.OUT'
        mainfile.write_text(
            '     1 : nkpt\n'
            '     2 : nstsv\n\n'
            '     1   0.000000000 0.000000000 0.000000000 : k-point, vkl\n'
            ' (state, eigenvalue and occupancy below)\n'
            '     1 -1.000000000 2.000000000\n'
            '     2  0.500000000 0.000000000\n\n\n'
        )

        parser = EigvalFileParser()
        parser.mainfile = str(mainfile)
        source = parser.to_dict()

        assert source['n_k_points'] == 1
        assert source['n_states'] == 2
        np.testing.assert_allclose(source['k_points'], [[0.0, 0.0, 0.0]])
        np.testing.assert_allclose(
            source['eigenvalues_occupancies'][0]['eigenvalues'], [[-1.0, 0.5]]
        )
        np.testing.assert_allclose(
            source['eigenvalues_occupancies'][0]['occupancies'], [[2.0, 0.0]]
        )

    def test_reads_spin_resolved_eigenvalues_and_occupancies(self):
        source = '1 -1.0 2.0\n2 0.5 0.0\n\n'

        parsed = str_to_eigenvalues(source)

        np.testing.assert_allclose(parsed['eigenvalues'], [[-1.0, 0.5]])
        np.testing.assert_allclose(parsed['occupancies'], [[2.0, 0.0]])

    def test_declares_expected_quantities(self):
        parser = EigvalFileParser()
        parser.init_quantities()

        names = {quantity.name for quantity in parser._quantities}
        assert names == {
            'k_points',
            'eigenvalues_occupancies',
            'n_k_points',
            'n_states',
        }
