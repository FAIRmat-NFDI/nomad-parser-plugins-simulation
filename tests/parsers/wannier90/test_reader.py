import numpy as np
import pytest

from nomad_simulation_parsers.parsers.wannier90.file_parsers import (
    HrParser,
    WInParser,
    WOutParser,
)


def read(parser, path, text):
    path.write_text(text)
    parser.mainfile = str(path)
    return parser.to_dict()


@pytest.mark.unit
class TestWOutParser:
    def test_reads_structure_kmesh_and_run_settings(self, tmp_path):
        source = read(
            WOutParser(),
            tmp_path / 'mock.wout',
            """
| Release: 3.1.0
a_1 -1.0 0.0 0.0
a_2 0.0 1.0 0.0
a_3 0.0 0.0 2.0
b_1 0.0 1.0 0.0
b_2 1.0 0.0 0.0
b_3 0.0 0.0 0.5
Fractional Coordinate
| Si 1 0.0 0.0 0.0 | 1.0 2.0 3.0 |
| O 1 0.25 0.25 0.25 | 4.0 5.0 6.0 |
PROJECTIONS
K-POINT GRID
Total points = 2
Grid size = 2 x 1 x 1
| 1 0.0 0.0 0.0 |
| 2 0.5 0.0 0.0 |
- MAIN
| Number of Wannier Functions : 2
| Number of input Bloch states : 4
| Total number of iterations : 10
| Convergence tolerence : 1.0E-12
""",
        )

        assert source['version'] == '3.1.0'
        assert len(source['lattice_vectors']) == 3
        assert source['lattice_vectors'][2].tolist() == [0, 0, 2]
        assert len(source['reciprocal_lattice_vectors']) == 3
        assert source['structure']['labels'][:2] == ['Si', 'O']
        assert source['structure']['positions'][1].tolist() == [4, 5, 6]
        assert source['k_mesh']['n_points'] == 2
        assert source['k_mesh']['grid'].tolist() == [2, 1, 1]
        assert source['k_mesh']['k_points'][1].tolist() == [0.5, 0.0, 0.0]
        assert source['Nwannier'] == 2
        assert source['Nband'] == 4
        assert source['Niter'] == 10
        assert source['conv_tol'] == pytest.approx(1.0e-12)

    def test_truncated_wout_keeps_available_header(self, tmp_path):
        assert (
            read(WOutParser(), tmp_path / 'truncated.wout', '|  Release: 3.1.0\n')[
                'version'
            ]
            == '3.1.0'
        )


@pytest.mark.unit
class TestWInParser:
    def test_reads_projections_and_fermi_energy(self, tmp_path):
        source = read(
            WInParser(),
            tmp_path / 'mock.win',
            """
fermi_energy = 5.25
begin projections
Cu:dx2-y2:x=1,0,0
end projections
""",
        )

        assert source['energy_fermi'] == pytest.approx(5.25)
        assert source['projections'] == [['Cu', 'dx2-y2', 'x=1,0,0']]


@pytest.mark.unit
class TestHrParser:
    def test_reads_hr_numeric_rows_with_exponents(self, tmp_path):
        source = read(
            HrParser(),
            tmp_path / 'mock_hr.dat',
            """
written on 1Jan2025 at 00:00:00
1
1
1 1 1 1 2.5E-01 -3.0E-01
""",
        )

        assert source['degeneracy_factors'][:3].tolist() == [1, 1, 1]
        assert len(source['hoppings']) == 1
        assert np.asarray(source['hoppings'][0]).tolist() == [
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            0.25,
            -0.3,
        ]
