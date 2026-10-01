from pathlib import Path

import pytest
from nomad.units import ureg

from nomad_simulation_parsers.parsers.octopus.file_parser import (
    ControlParser,
    EigenvalueParser,
    InfoParser,
    InpParser,
    LogParser,
    OutParser,
)
from tests.parsers.common import approx, assert_approx


@pytest.mark.unit
class TestInpReader:
    def test_reads_input_options_and_coordinates(self, tmp_path):
        mainfile = tmp_path / 'inp'
        mainfile.write_text(
            'Units = eV_Angstrom\n'
            'SpinComponents = spin_polarized\n'
            '%Coordinates\n'
            'Si | 0.0 | 0.0 | 0.0\n'
            'Si | 1.0 | 1.0 | 1.0\n'
            '%\n'
        )

        parser = InpParser()
        parser.mainfile = str(mainfile)

        assert parser.info['Units'] == 'eV_Angstrom'
        symbols, coordinates = parser.get_coordinates()
        assert symbols == ['Si', 'Si']
        assert_approx(coordinates, [[0, 0, 0], [1, 1, 1]])

    def test_evaluates_octopus_constants_and_units(self):
        parser = InpParser()
        parser._info = {}
        parser._keys_mapping = {}

        assert parser.evaluate_value('yes') is True
        assert parser.evaluate_value('2*angstrom') == approx(
            2 / (1 * ureg.bohr).to('angstrom').magnitude
        )


@pytest.mark.unit
class TestInfoReader:
    def test_reads_scf_kpoints(self):
        parser = InfoParser()
        parser.mainfile = str(
            Path(__file__).resolve().parents[2]
            / 'data'
            / 'octopus'
            / 'Si_scf'
            / 'static'
            / 'info'
        )
        sampling = parser.to_dict()['brillouin_zone_sampling']

        assert sampling['kgrid'].tolist() == [4, 4, 4]
        assert sampling['n_kpoints'] == 64
        assert sampling['kpoints'][11][2] == approx(0.25)


@pytest.mark.unit
class TestEigenvalueReader:
    def test_reads_eigenvalue_blocks(self):
        parser = EigenvalueParser()
        parser.mainfile = str(
            Path(__file__).resolve().parents[2]
            / 'data'
            / 'octopus'
            / 'Si_scf'
            / 'static'
            / 'info'
        )
        eigenvalues = parser.to_dict()['eigenvalues']
        entries = [entry for entry in eigenvalues['eigenvalues'] if entry is not None]

        assert eigenvalues['unit'] == 'hartree'
        assert len(entries) == 18
        assert entries[0][0].tolist() == [0.0, 0.0, 0.0]
        assert entries[0][1].shape == (8, 1)
        assert entries[0][2][0][0] == approx(2.0)


@pytest.mark.unit
class TestControlReader:
    def test_reads_control_lines_and_comments(self, tmp_path):
        mainfile = tmp_path / 'control'
        mainfile.write_text('Spacing = 0.283459 # default\nSpinComponents = 2\n')

        parser = ControlParser()
        parser.mainfile = str(mainfile)

        assert parser.info['Spacing'] == approx(0.283459)
        assert parser.info['SpinComponents'] == 2.0


@pytest.mark.unit
class TestLogReader:
    def test_reads_coordinates_from_parser_log(self):
        parser = LogParser()
        parser.mainfile = str(
            Path(__file__).resolve().parents[2]
            / 'data'
            / 'octopus'
            / 'Si_scf'
            / 'exec'
            / 'parser.log'
        )

        symbols, coordinates = parser.get_coordinates()

        assert symbols == ['Si', 'Si', 'Si', 'Si']
        assert coordinates.shape == (4, 3)


@pytest.mark.unit
class TestOutReader:
    def test_reads_output_header_grid_and_theory_level(self):
        parser = OutParser()
        parser.mainfile = str(
            Path(__file__).resolve().parents[2]
            / 'data'
            / 'octopus'
            / 'Si_scf'
            / 'stdout.txt'
        )
        source = parser.to_dict()

        assert source['header']['options'][0] == ['Version', 'wolfi']
        assert source['grid']['npbc'] == 3
        assert source['theory_level']['theory_level'] == 'dft'
        assert source['theory_level']['correlation'] == 'Perdew & Zunger (Modified)'
