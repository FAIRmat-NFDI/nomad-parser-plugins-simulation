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
