import numpy as np
import pytest

from nomad_simulation_parsers.parsers.orca.text_parser import OutReader

# Minimal scaffolding the nested block regexes require: a Single Point section
# wrapping an `ORCA SCF` block. Individual feature blocks are appended per test.
SINGLE_POINT_SCF = '* Single Point Calculation *\n****\nORCA SCF\n' + '-' * 40 + '\n'


def _unwrap(value):
    return value._results if hasattr(value, '_results') else value or {}


def _read(tmp_path, fragment: str) -> dict:
    mainfile = tmp_path / 'mol.out'
    mainfile.write_text(fragment)
    reader = OutReader()
    reader.mainfile = str(mainfile)
    return reader.to_dict()


def _self_consistent(source: dict) -> dict:
    return _unwrap(_unwrap(source.get('single_point')).get('self_consistent'))


def _scf_settings(source: dict) -> dict:
    return _unwrap(_self_consistent(source).get('scf_settings'))


@pytest.mark.unit
class TestOrcaReader:
    def test_extracts_scf_settings(self, tmp_path):
        source = _read(
            tmp_path,
            SINGLE_POINT_SCF + 'SCF SETTINGS\n------------\n'
            'Hartree-Fock type      HFTyp           .... RHF\n'
            'Total Charge           Charge          .... 0\n'
            'Multiplicity           Mult            .... 1\n'
            '------------\n',
        )
        scf = _scf_settings(source)
        assert scf['hf_type'] == 'RHF'
        assert scf['total_charge'] == 0.0
        assert scf['multiplicity'] == 1.0

    def test_extracts_orbital_energies_columns(self, tmp_path):
        source = _read(
            tmp_path,
            SINGLE_POINT_SCF + 'ORBITAL ENERGIES\n----------------\n'
            '  NO   OCC          E(Eh)            E(eV) \n'
            '   0   2.0000     -20.550000      -559.0000 \n'
            '   1   0.0000       0.500000        13.6000 \n'
            '\n',
        )
        # str_operation yields rows of [NO, OCC, E(Eh), E(eV)]; repeats -> list.
        table = _self_consistent(source)['orbital_energies'][0]
        np.testing.assert_allclose(table[0], [0.0, 2.0, -20.55, -559.0])
        assert table[1][1] == 0.0  # the unoccupied orbital's occupation

    def test_extracts_cartesian_coordinates(self, tmp_path):
        source = _read(
            tmp_path,
            SINGLE_POINT_SCF
            + 'CARTESIAN COORDINATES (ANGSTROEM)\n---------------------------------\n'
            '  O      0.000000    0.000000    0.117790\n'
            '  H      0.000000    0.755450   -0.471160\n'
            '\n',
        )
        coordinates = _unwrap(source.get('single_point'))['cartesian_coordinates']
        assert coordinates[0] == 'O'
        assert coordinates[3] == pytest.approx(0.11779)

    def test_keeps_partial_data_when_scf_block_is_truncated(self, tmp_path):
        # No closing `------------`, so the scf_settings block never closes and is
        # not captured -- a truncated file must not raise.
        source = _read(
            tmp_path,
            SINGLE_POINT_SCF + 'SCF SETTINGS\n------------\n'
            'Hartree-Fock type      HFTyp           .... RHF\n',
        )
        assert 'hf_type' not in _scf_settings(source)

    @pytest.mark.parametrize('newline', ['\n', '\r\n'])
    def test_line_ending_invariance(self, tmp_path, newline):
        fragment = (
            SINGLE_POINT_SCF + 'SCF SETTINGS\n------------\n'
            'Total Charge           Charge          .... 0\n'
            'Multiplicity           Mult            .... 1\n'
            '------------\n'
        ).replace('\n', newline)
        scf = _scf_settings(_read(tmp_path, fragment))
        assert scf['total_charge'] == 0.0
        assert scf['multiplicity'] == 1.0

    def test_ignores_malformed_numeric_charge(self, tmp_path):
        source = _read(
            tmp_path,
            SINGLE_POINT_SCF + 'SCF SETTINGS\n------------\n'
            'Total Charge           Charge          .... not-a-number\n'
            'Multiplicity           Mult            .... 1\n'
            '------------\n',
        )
        scf = _scf_settings(source)
        assert 'total_charge' not in scf  # malformed value is not captured
        assert scf['multiplicity'] == 1.0
