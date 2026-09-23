import pytest
from netCDF4 import Dataset

from nomad_simulation_parsers.parsers.yambo.file_parsers import (
    InputParser,
    MainfileParser,
    NetCDFParser,
)


@pytest.mark.unit
class TestMainfileReader:
    def test_reads_report_values(self, tmp_path):
        mainfile = tmp_path / 'r-example'
        mainfile.write_text(
            '\nCORE Variables Setup\n'
            'Energies & Occupations\n'
            '  [X] Fermi Level : 5.110763 [eV]\n'
            '  [X] Valence Band Max : 0.000000 [eV]\n'
            '  [X] Conduction Band Min : 3.878048 [eV]\n'
            '  Filled Bands : 8\n'
            '  Empty Bands : 9 100\n'
            '  Electronic Temp. : 0.000000 300.000000 [K]\n'
            '  Bosonic Temp. : 0.000000 150.000000 [K]\n'
            '  Finite Temperature mode: yes\n'
            '  El. density : 0.125 [electrons/bohr^3]\n'
            '  Indirect Gaps : 1.100000 2.200000 [eV]\n'
            '  Direct Gaps : 1.300000 2.400000 [eV]\n'
            '  Indirect Gap : 1.100000 [eV]\n'
            '  Direct Gap : 1.300000 [eV]\n'
            '  Direct Gap localized at k-point : 2\n'
            '  Indirect Gap between k-points : 1 2\n'
            '[0]\n'
        )
        parser = MainfileParser()
        parser.mainfile = str(mainfile)
        source = parser.to_dict()
        source = source['core_variables_setup']['energies_occupations']

        assert source['fermi'].magnitude == pytest.approx(5.110763)
        assert source['valence'].magnitude == pytest.approx(0.0)
        assert source['conduction'].magnitude == pytest.approx(3.878048)
        assert source['x_yambo_filled_bands'].tolist() == [1, 8]
        assert source['x_yambo_empty_bands'].tolist() == [9, 100]
        assert source['x_yambo_electronic_temperature'].magnitude == pytest.approx(300)
        assert source['x_yambo_bosonic_temperature'].magnitude == pytest.approx(150)
        assert source['x_yambo_finite_temperature_mode'] is True
        assert source['x_yambo_electronic_density'] == pytest.approx(0.125)
        assert source['x_yambo_indirect_gaps'].magnitude.tolist() == [1.1, 2.2]
        assert source['x_yambo_direct_gaps'].magnitude.tolist() == [1.3, 2.4]
        assert source['x_yambo_indirect_gap'].magnitude == pytest.approx(1.1)
        assert source['x_yambo_direct_gap'].magnitude == pytest.approx(1.3)
        assert source['x_yambo_direct_gap_kpoint'] == 2
        assert source['x_yambo_indirect_gap_kpoints'].tolist() == [1, 2]


@pytest.mark.unit
class TestInputReader:
    def test_reads_input_key_values(self, tmp_path):
        mainfile = tmp_path / 'input'
        mainfile.write_text(
            '\nFFTGvecs = 10 Ry\n'
        )
        parser = InputParser()
        parser.mainfile = str(mainfile)
        source = parser.to_dict()

        assert source['key_value'][0] == ['FFTGvecs', 10, 'Ry']


@pytest.mark.unit
class TestNetCDFReader:
    def test_reads_variables(self, tmp_path):
        mainfile = tmp_path / 'sample.nc'
        with Dataset(mainfile, 'w') as dataset:
            dataset.createDimension('points', 2)
            values = dataset.createVariable('K-POINTS', 'f8', ('points',))
            values[:] = [0.0, 0.5]

        parser = NetCDFParser()
        parser.mainfile = str(mainfile)
        parser.parse()

        assert parser.results['K-POINTS'].tolist() == [0.0, 0.5]
