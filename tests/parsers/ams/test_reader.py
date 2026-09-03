from pathlib import Path

import pytest

from nomad_simulation_parsers.parsers.ams.file_parser import OutParser, RKFParser

DATA_DIR = Path(__file__).resolve().parents[2] / 'data' / 'ams'

MINIMAL_OUTPUT = """
* r2024 2024-01-01
SINGLE POINT CALCULATION *
Index Symbol       x (bohr)       y (bohr)       z (bohr)
      1      H     0.0     0.0     0.0
      2      H     0.0     0.0     1.4
Lattice vectors (bohr)
      1    5.0     0.0     0.0
      2    0.0     5.0     0.0
      3    0.0     0.0     5.0
Total System Charge             0.00000
Energy (hartree)            -1.123
Timing
"""


def read_source(path):
    parser = OutParser()
    parser.mainfile = str(path)
    return parser.to_dict()


def read_rkf_source(path):
    return RKFParser(mainfile=str(path)).parse().results


@pytest.mark.unit
class TestOutReader:
    def test_extracts_header_inputs_and_results(self, tmp_path):
        mainfile = tmp_path / 'minimal.out'
        mainfile.write_text(MINIMAL_OUTPUT)

        source = read_source(mainfile)

        assert source['program_version'] == '2024 2024-01-01'
        single_point = source['single_point']
        assert single_point['labels_positions'][0] == ['H', 'H']
        assert single_point['labels_positions'][1].units == 'atomic_unit_of_length'
        assert single_point['labels_positions'][1].magnitude[1, 2] == 1.4
        assert single_point['lattice_vectors'].magnitude[2, 2] == 5.0
        assert single_point['total_charge'].magnitude == 0.0
        assert single_point['energy_total'].units == 'hartree'
        assert single_point['energy_total'].magnitude == -1.123


@pytest.mark.unit
class TestRKFReader:
    def test_extracts_adf_single_point_results(self):
        source = read_rkf_source(
            DATA_DIR / 'adf_SP' / 'adf_SP.results' / 'ams.rkf'
        )['single_point']

        assert source['program_version'] == '98925 2021-11-25'
        assert source['program_x_ams_name'] == 'ADF'
        assert source['model_parameters'] == {
            'spin': False,
            'dft_potential': {'LDA': 'VWN'},
        }
        assert source['labels_positions'][0] == ['O', 'O']
        assert source['labels_positions'][1].units == 'atomic_unit_of_length'
        assert source['labels_positions'][1].magnitude[1, 0] == pytest.approx(
            -3.1867335262123166
        )
        assert source['energy_total'].to('hartree').magnitude == pytest.approx(
            -0.25256548879865254
        )
        assert source['energies']['xc'].to('hartree').magnitude == pytest.approx(
            -0.25252261177051805
        )
        assert source['fermi_energy'].to('hartree').magnitude == pytest.approx(
            -0.28159925821297915
        )
        assert source['atomic_charges'][1] == pytest.approx(-1.2850914998807639e-08)
