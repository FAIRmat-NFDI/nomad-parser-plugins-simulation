import numpy as np
import pytest

from nomad_simulation_parsers.parsers.lobster.file_parser import (
    CHARGEParser,
    COXPCARParser,
    ICOXPLISTParser,
    OutParser,
)
from nomad_simulation_parsers.parsers.vasp.doscar_parser import DOSCARFileParser
from tests.parsers.common import assert_approx


@pytest.mark.unit
class TestOutParser:
    def test_extracts_metadata_basis_and_spillings(self, tmp_path):
        mainfile = tmp_path / 'lobsterout'
        mainfile.write_text(
            'LOBSTER v5.1.1\n'
            'starting on host node on 2025-01-02 at 03:04:05 UTC\n'
            'detecting used PAW program... VASP\n'
            'setting up local basis functions...\n'
            'Fe (3d) 4s 4p\n\n'
            'spillings:\n'
            'abs. total spilling: 1.23%\n'
            'abs. charge spilling: 0.45%\n'
            'finished in 2 seconds\n'
        )

        parser = OutParser()
        parser.mainfile = str(mainfile)
        source = parser.to_dict()

        assert source['program_version'] == '5.1.1'
        assert source['datetime'] == '2025-01-02 at 03:04:05'
        assert source['x_lobster_code'] == 'VASP'
        assert source['x_lobster_basis']['x_lobster_basis_species'] == [
            ['Fe', '3d', '4s', '4p']
        ]
        assert source['spilling'][0]['abs_total_spilling'] == pytest.approx(1.23)
        assert source['finished'] == 2


@pytest.mark.unit
class TestICOXPLISTParser:
    def test_extracts_icoxp_list_data(self, tmp_path):
        mainfile = tmp_path / 'ICOHPLIST.lobster'
        mainfile.write_text('1 Fe1 Fe2 2.830 -0.10 0 0 1\n')

        parser = ICOXPLISTParser()
        parser.mainfile = str(mainfile)
        source = parser.to_dict()

        assert len(source['data']) == 1
        np.testing.assert_array_equal(
            source['data'][0][:4, 0], ['1', 'Fe1', 'Fe2', '2.830']
        )
        assert source['data'][0][4, 0] == '-0.10'


@pytest.mark.unit
class TestCOXPCARParser:
    def test_extracts_non_spin_polarized_coxpcar_data(self, tmp_path):
        mainfile = tmp_path / 'COHPCAR.lobster'
        mainfile.write_text(
            'No.1:Fe1->Fe2(2.830)\n-1.0 1.0 2.0 3.0 4.0\n 0.0 5.0 6.0 7.0 8.0\n'
        )

        parser = COXPCARParser()
        parser.mainfile = str(mainfile)
        source = parser.to_dict()

        assert source['bond_pairs'] == [('1', 'Fe1', 'Fe2', '2.830')]
        assert_approx(source['energy'], [-1.0, 0.0])
        assert_approx(source['total_coxp'], [[1.0, 5.0]])
        assert_approx(source['total_icoxp'], [[2.0, 6.0]])
        assert_approx(source['pair_coxp'], [[[3.0, 7.0]]])
        assert_approx(source['pair_icoxp'], [[[4.0, 8.0]]])


@pytest.mark.unit
class TestDOSCARFileParser:
    def test_extracts_doscar_total_density_of_states(self, tmp_path):
        mainfile = tmp_path / 'DOSCAR.lobster'
        mainfile.write_text(
            'header 1\n'
            'header 2\n'
            'header 3\n'
            'header 4\n'
            'header 5\n'
            '-10.0 10.0 3 0.0 0.0\n'
            '-1.0 1.0 2.0 3.0 4.0\n'
            ' 0.0 5.0 6.0 7.0 8.0\n'
            ' 1.0 9.0 10.0 11.0 12.0\n'
        )

        parser = DOSCARFileParser()
        parser.mainfile = str(mainfile)
        source = parser.to_dict()

        assert source['e_fermi'] == 0.0
        assert_approx(source['energies'], [-1.0, 0.0, 1.0])
        assert_approx(
            source['total_dos'],
            [[1.0, 5.0, 9.0], [2.0, 6.0, 10.0], [3.0, 7.0, 11.0], [4.0, 8.0, 12.0]],
        )


@pytest.mark.unit
class TestCHARGEParser:
    def test_extracts_charge_table(self, tmp_path):
        mainfile = tmp_path / 'CHARGE.lobster'
        mainfile.write_text('1 Fe 0.1 0.2\n2 Fe -0.1 -0.2\ntotal 0.0 0.0\n')

        parser = CHARGEParser()
        parser.mainfile = str(mainfile)
        source = parser.to_dict()

        np.testing.assert_array_equal(source['charges']['indices'], [1, 2])
        assert_approx(source['charges']['mulliken'], [0.1, -0.1])
        assert_approx(source['charges']['loewdin'], [0.2, -0.2])
        assert_approx(source['total'], [0.0, 0.0])
