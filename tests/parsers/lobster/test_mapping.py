import numpy as np
import pytest

from nomad_simulation_parsers.parsers.lobster.parser import (
    LobsterCHARGEParser,
    LobsterCOXPCARParser,
    LobsterDOSCARParser,
    LobsterICOXPLISTParser,
    LobsterMainfileParser,
)
from tests.parsers.common import assert_approx


@pytest.fixture(scope='module')
def mainfile_parser():
    return LobsterMainfileParser()


@pytest.fixture(scope='module')
def coxpcar_parser():
    return LobsterCOXPCARParser()


@pytest.fixture(scope='module')
def icoxplist_parser():
    return LobsterICOXPLISTParser()


@pytest.fixture(scope='module')
def charge_parser():
    return LobsterCHARGEParser()


@pytest.fixture(scope='module')
def doscar_parser():
    return LobsterDOSCARParser()


@pytest.mark.unit
class TestLobsterMainfileParserMapping:
    @pytest.mark.parametrize(
        ('basis', 'expected'),
        [
            ('pbevaspfit2015', 'pbeVaspFit2015'),
            ('PBEVASPfit2015', 'PBEVASPFit2015'),
            ('bunge', 'Bunge'),
            ('KOGA', 'KOGA'),
        ],
    )
    def test_maps_basis_names_and_spillings(self, mainfile_parser, basis, expected):
        assert mainfile_parser.to_basis_set(basis) == expected
        assert_approx(
            mainfile_parser.get_spilling(
                [{'abs_total_spilling': 1.2}, {'abs_total_spilling': 3.4}],
                type='total',
            ),
            [1.2, 3.4],
        )

    def test_maps_lobster_datetime_to_unix_time(self, mainfile_parser):
        assert mainfile_parser.to_unix_time('2025-01-02 at 03:04:05') == pytest.approx(
            1735787045.0
        )


@pytest.mark.unit
class TestLobsterCOXPCARParserMapping:
    def test_maps_coxp_pair_helpers(self, coxpcar_parser):
        source = {
            'coxp_pairs': [
                ('1', 'Fe1', 'Fe2', '2.83'),
                ('1', 'Fe1', 'Fe2', '2.83'),
                ('3', 'Fe2', 'Fe1', '2.83'),
            ]
        }

        assert coxpcar_parser.get_bond_label(source['coxp_pairs'][0]) == 'Fe1 -> Fe2'
        assert coxpcar_parser.get_orbital_pairs(
            ('1', 'Fe1', 'Fe2', '2.83'), source
        ) == [1]
        assert coxpcar_parser.get_atom_label(0, source) == 'Fe1'
        assert coxpcar_parser.get_atom_label(0, source, atom=1) == 'Fe2'

    def test_maps_atom_values_and_orbital_values(self, coxpcar_parser):
        values = np.array([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]])
        pairs = [
            ('1', 'Fe1', 'Fe2', '2.83'),
            ('1', 'Fe1[1]', 'Fe2', '2.83'),
            ('2', 'Fe2', 'Fe1', '2.83'),
        ]
        source = {'pair_icoxp': values}

        assert_approx(
            coxpcar_parser.get_atom_value(values, pairs), [[1.0, 2.0], [5.0, 6.0]]
        )
        assert_approx(coxpcar_parser.get_orbital_value(1, source), [3.0, 4.0])


@pytest.mark.unit
class TestLobsterICOXPLISTParserMapping:
    @staticmethod
    def source():
        return {
            'data': [
                np.array(
                    [
                        ['1', '1', '2'],
                        ['Fe1', 'Fe1', 'Fe2'],
                        ['Fe2', 'Fe2', 'Fe1'],
                        ['2.83', '2.83', '2.45'],
                        ['0', '0', '1'],
                        ['0', '0', '0'],
                        ['0', '0', '0'],
                        ['-0.10', '-0.20', '-0.30'],
                    ]
                )
            ]
        }

    def test_maps_unique_atom_pairs_and_bond_properties(self, icoxplist_parser):
        source = self.source()

        assert icoxplist_parser.get_atom_pair_indices(source) == [0, 2]
        assert icoxplist_parser.get_atom_labels(source) == ['Fe1', 'Fe2']
        assert icoxplist_parser.get_atom_labels(source, atom=1) == ['Fe2', 'Fe1']
        assert_approx(icoxplist_parser.get_atom_distances(source), [2.83, 2.45])
        assert_approx(icoxplist_parser.get_translations(source), [[0, 0, 0], [1, 0, 0]])
        assert_approx(
            icoxplist_parser.get_integrated_coxp_at_fermi_level(source),
            [[-0.10, -0.30]],
        )

    def test_maps_optional_bond_indices_and_missing_translations(
        self, icoxplist_parser
    ):
        source = self.source()
        source['data'][0][4] = ['translation', 'translation', 'translation']
        source['data'][0][-1] = ['1', '2', '3']

        assert icoxplist_parser.get_translations(source) is None
        assert_approx(icoxplist_parser.get_bonds(source), [1, 2, 3])


@pytest.mark.unit
class TestLobsterCHARGEParserMapping:
    def test_maps_charge_contributions(self, charge_parser):
        source = {
            'charges': {
                'symbols': ['Fe', 'Fe'],
                'mulliken': np.array([0.1, -0.1]),
                'loewdin': np.array([0.2, -0.2]),
            },
            'total': np.array([0.0, 0.0]),
        }

        assert charge_parser.get_contributions(source['charges'], 'mulliken') == [
            {'symbol': 'Fe', 'value': 0.1},
            {'symbol': 'Fe', 'value': -0.1},
        ]
        assert charge_parser.get_charges(source)[0]['kind'] == 'mulliken'
        assert charge_parser.get_charges(source)[1]['contributions'][1][
            'value'
        ] == pytest.approx(-0.2)


@pytest.mark.unit
class TestLobsterDOSCARParserMapping:
    def test_maps_spin_dos_and_clips_negative_values(self, doscar_parser):
        dos = doscar_parser.get_dos(
            np.array([[-1.0, 2.0], [3.0, -4.0]]),
            [],
        )

        assert len(dos) == 1
        assert dos[0]['spin'] is None
        assert_approx(dos[0]['dos'], [0.0, 2.0])
        assert_approx(dos[0]['integrated'], [3.0, -4.0])
