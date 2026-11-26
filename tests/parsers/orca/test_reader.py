import numpy as np
import pytest

from nomad_simulation_parsers.parsers.orca.text_parser import OutReader

# Minimal scaffolding the nested block regexes require: a Single Point section
# wrapping an `ORCA SCF` block. Individual feature blocks are appended per test.
SINGLE_POINT_SCF = '* Single Point Calculation *\n****\nORCA SCF\n' + '-' * 40 + '\n'
GEOMETRY_OPTIMIZATION = '* Geometry Optimization Run *\n****\n'


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

    @pytest.mark.parametrize(
        ('fragment', 'headers', 'rows'),
        [
            pytest.param(
                SINGLE_POINT_SCF + 'SCF ITERATIONS\n--------------\n'
                'ITER Energy Delta-E Max-DP RMS-DP [F,P] Damp\n'
                '  0 -76.28 0.00 0.04 0.005 0.29 0.85\n'
                'ITER Energy Delta-E Grad Rot Max-DP RMS-DP\n'
                '  1 -76.31 -0.03 0.001 0.002 0.003 0.0004\n',
                [
                    'Energy Delta-E Max-DP RMS-DP [F,P] Damp',
                    'Energy Delta-E Grad Rot Max-DP RMS-DP',
                ],
                [
                    [-76.28, 0.0, 0.04, 0.005, 0.29, 0.85],
                    [-76.31, -0.03, 0.001, 0.002, 0.003, 0.0004],
                ],
                id='orca4',
            ),
            pytest.param(
                '* Single Point Calculation *\n****\n'
                '----------------------\nSHARK INTEGRAL PACKAGE\n'
                '----------------------\n'
                '----------------------------------------'
                '----------------------------------------\n'
                'Iteration Energy (Eh) Delta-E RMSDP MaxDP DIISErr Damp Time(sec)\n'
                '  1 -76.28 0.00e+00 5.0e-03 4.0e-02 2.9e-01 0.85 0.2\n'
                'Iteration Energy (Eh) Delta-E RMSDP MaxDP MaxGrad Time(sec)\n'
                '  2 -76.31 -3.0e-02 4.0e-04 3.0e-03 1.0e-03 0.1\n',
                [
                    'Energy (Eh) Delta-E RMSDP MaxDP DIISErr Damp Time(sec)',
                    'Energy (Eh) Delta-E RMSDP MaxDP MaxGrad Time(sec)',
                ],
                [
                    [-76.28, 0.0, 0.005, 0.04, 0.29, 0.85, 0.2],
                    [-76.31, -0.03, 0.0004, 0.003, 0.001, 0.1],
                ],
                id='orca6',
            ),
        ],
    )
    def test_extracts_scf_iteration_headers_and_rows(
        self, tmp_path, fragment, headers, rows
    ):
        source = _read(
            tmp_path,
            fragment
            + '********************\nTotal SCF time: 0 days 0 hours 0 min 1 sec\n',
        )
        blocks = _self_consistent(source)['scf_iterations']
        for block, header, row in zip(blocks, headers, rows, strict=True):
            parsed = _unwrap(block)
            assert list(parsed['header']) == header.split()
            np.testing.assert_allclose(parsed['values'], [row])

    def test_extracts_geometry_convergence_values_for_each_cycle(self, tmp_path):
        source = _read(
            tmp_path,
            GEOMETRY_OPTIMIZATION + 'Max. Gradient TolMAXG .... 1.0000e-04 Eh/bohr\n'
            '* GEOMETRY OPTIMIZATION CYCLE   1 *\n****\n'
            'Energy change  -0.0006358496  0.0000010000  NO\n'
            'RMS gradient    0.0024928050  0.0000300000  NO\n'
            'MAX gradient    0.0042659240  0.0001000000  NO\n'
            'RMS step        0.0088383489  0.0006000000  NO\n'
            'MAX step        0.0148658788  0.0010000000  NO\n'
            '* GEOMETRY OPTIMIZATION CYCLE   2 *\n****\n'
            'Energy change  -0.0000000369  0.0000010000  YES\n'
            'RMS gradient    0.0000100000  0.0000300000  YES\n'
            'MAX gradient    0.0000177027  0.0001000000  YES\n'
            'RMS step        0.0000304880  0.0006000000  YES\n'
            'MAX step        0.0000369167  0.0010000000  YES\n'
            '*** OPTIMIZATION RUN DONE ***\n',
        )
        geometry = _unwrap(source['geometry_optimization'])
        first, second = map(_unwrap, geometry['cycle'])
        label, tolerance = geometry['max_gradient_tol']
        assert label.magnitude == 'TolMAXG'
        assert float(tolerance.magnitude) == pytest.approx(1e-4)
        expected_values = {
            'energy_change': [-0.0006358496, -0.0000000369],
            'rms_gradient': [0.0024928050, 0.0000100000],
            'max_gradient': [0.0042659240, 0.0000177027],
            'rms_displacement': [0.0088383489, 0.0000304880],
            'max_displacement': [0.0148658788, 0.0000369167],
        }
        for key, expected in expected_values.items():
            assert [
                cycle[f'geom_opt_{key}'].magnitude for cycle in (first, second)
            ] == pytest.approx(expected, abs=1e-12)

    @pytest.mark.parametrize(
        ('banner', 'converged', 'unconverged'),
        [
            (
                '*** THE OPTIMIZATION HAS CONVERGED ***',
                ['THE', 'OPTIMIZATION', 'HAS', 'CONVERGED'],
                None,
            ),
            (
                'ERROR !!!\nThe optimization did not converge but reached the maximum\n'
                'number of optimization cycles.',
                None,
                ['The', 'optimization', 'did', 'not', 'converge'],
            ),
            (
                'The optimization has not yet converged - '
                'more geometry cycles are needed',
                None,
                None,
            ),
        ],
        ids=['converged', 'unconverged', 'in-progress'],
    )
    def test_extracts_geometry_optimization_banners(
        self, tmp_path, banner, converged, unconverged
    ):
        source = _read(tmp_path, GEOMETRY_OPTIMIZATION + banner + '\n')
        geometry = _unwrap(source['geometry_optimization'])
        assert geometry.get('is_converged') == converged
        assert geometry.get('is_not_converged') == unconverged
