import numpy as np
import pytest

from nomad_simulation_parsers.parsers.wannier90.parser import (
    WBandTextParser,
    WDosTextParser,
    WHrTextParser,
    WInTextParser,
    WOutTextParser,
)
from tests.parsers.common import approx, assert_approx


@pytest.mark.unit
class TestWOutTextParser:
    @pytest.fixture
    def parser(self):
        return WOutTextParser()

    def test_maps_lattice_vectors_and_periodic_boundary_conditions(self, parser):
        vectors = [np.array([1, 0, 0]), np.array([0, 2, 0]), np.array([0, 0, 3])]
        assert_approx(parser.get_lattice_vectors([np.zeros(3), *vectors]), vectors)
        assert parser.get_pbc([np.ones(3)]) == [True, True, True]
        assert parser.get_pbc(None) == [False, False, False]

    def test_maps_localization_and_kpoints(self, parser):
        assert parser.is_maximally_localized(1000)
        assert not parser.is_maximally_localized(0)
        assert parser.is_maximally_localized(0, default=2)
        assert_approx(
            parser.get_kpoints(np.array([[0, 0, 0], [1, 1, 1], [2, 2, 2]])),
            [[0, 0, 0], [2, 2, 2]],
        )

    def test_maps_k_line_path(self, parser):
        result = parser.get_k_line_path(
            {
                'high_symm_name': [['G', 'X'], ['X', 'G']],
                'high_symm_value': [[0, 0, 0, 0.5, 0, 0], [0.5, 0, 0, 0, 0, 0]],
            }
        )
        assert result['names'] == ['G', 'X', 'G']


@pytest.mark.unit
class TestWBandTextParser:
    @pytest.fixture
    def parser(self):
        return WBandTextParser()

    def test_maps_band_data(self, parser):
        assert_approx(
            parser.get_data(np.array([[0.0, 1.0], [0.5, 2.0]])),
            [[1.0], [2.0]],
        )

    def test_maps_band_structure_with_metadata(self, parser):
        parser._highest_occupied = 1.5
        parser._k_path = {'names': ['G'], 'values': [[0, 0, 0]]}
        mapped_band = parser.get_band_structure(np.array([[0.0, 1.0], [0.5, 2.0]]))
        assert_approx(mapped_band['value'], [[1.0], [2.0]])
        assert mapped_band['highest_occupied'] == approx(1.5)

    def test_maps_band_structure_without_optional_metadata(self, parser):
        mapped_band = parser.get_band_structure(np.array([[0.0, 1.0], [0.5, 2.0]]))
        assert list(mapped_band) == ['value']


@pytest.mark.unit
class TestWDosTextParser:
    @pytest.fixture
    def parser(self):
        return WDosTextParser()

    def test_maps_dos_array_with_energy_origin(self, parser):
        parser._energies_origin = 11.375
        mapped_dos = parser.get_dos(np.array([[1.0, 2.0], [3.0, 4.0]]))
        assert_approx(mapped_dos['energies'], [1.0, 3.0])
        assert_approx(mapped_dos['value'], [2.0, 4.0])
        assert mapped_dos['energies_origin'] == approx(11.375)

    def test_maps_dos_array_without_energy_origin(self, parser):
        mapped_dos = parser.get_dos(np.array([[1.0, 2.0]]))
        assert_approx(mapped_dos['energies'], [1.0])
        assert_approx(mapped_dos['value'], [2.0])
        assert 'energies_origin' not in mapped_dos


@pytest.mark.unit
class TestWHrTextParser:
    @pytest.fixture
    def parser(self):
        return WHrTextParser()

    def test_maps_hoppings_and_crystal_fields(self, parser):
        mapped = parser.get_hoppings(
            {
                'degeneracy_factors': np.array([1, 1, 1]),
                'hoppings': np.array([1, 0, 0, 1, 1, 2.5, -0.25]),
                'n_orbitals': 1,
            },
            ws=True,
        )
        assert_approx(mapped['crystal_fields'], [2.5])
        assert mapped['hoppings'].shape == (1, 1, 1)
        assert mapped['hoppings'][0, 0, 0] == approx(2.5 - 0.25j)

    def test_maps_hoppings_without_wigner_seitz_points(self, parser):
        mapped = parser.get_hoppings(
            {
                'degeneracy_factors': np.array([1, 1, 1]),
                'hoppings': np.array([1, 0, 0, 1, 1, 2.5, -0.25]),
                'n_orbitals': 1,
            }
        )
        assert 'ws_points' not in mapped
        assert 'n_ws_points' not in mapped


@pytest.mark.unit
class TestWInTextParser:
    @pytest.fixture
    def parser(self):
        return WInTextParser()

    def test_maps_projection_and_orbital_quantum_numbers(self, parser):
        assert parser.get_projections([['Cu', 'dx2-y2']]) == [
            {'projection': ['Cu', 'dx2-y2']}
        ]

    def test_maps_branch_indices_by_integer_and_label(self, parser):
        positions = [np.array([0.0, 0.0, 0.0]), np.array([1.0, 0.0, 0.0])]
        labels = ['Cu', 'O']
        lattice = np.eye(3)
        assert parser.get_branch_label_indices(1, positions, labels, lattice) == {
            'label': '',
            'indices': [1],
        }
        assert parser.get_branch_label_indices('O', positions, labels, lattice) == {
            'label': 'O',
            'indices': [1],
        }

    def test_maps_branch_indices_by_cartesian_and_fractional_coordinates(self, parser):
        positions = [np.array([0.0, 0.0, 0.0]), np.array([1.0, 2.0, 3.0])]
        labels = ['Cu', 'O']
        lattice = np.diag([1.0, 2.0, 3.0])
        assert parser.get_branch_label_indices(
            'c=1,2,3', positions, labels, lattice
        ) == {'label': 'O', 'indices': [1]}
        assert parser.get_branch_label_indices(
            'f=1,1,1', positions, labels, lattice
        ) == {'label': 'O', 'indices': [1]}
        assert parser.get_branch_label_indices(None, positions, labels, lattice) is None

    @pytest.mark.parametrize(
        ('symbol', 'quantum_numbers'),
        [
            ('s', (0, 0)),
            ('px', (1, -1)),
            ('py', (1, 0)),
            ('pz', (1, 1)),
            ('dz2', (2, 0)),
            ('dxz', (2, 1)),
            ('dyz', (2, -1)),
            ('dx2-y2', (2, 2)),
            ('dxy', (2, -2)),
            ('fz3', (3, 0)),
            ('fxz2', (3, 1)),
            ('fyz2', (3, -1)),
            ('fz(x2-y2)', (3, 2)),
            ('fxyz', (3, -2)),
            ('fx(x2-3y2)', (3, 3)),
            ('fy(3x2-y2)', (3, -3)),
        ],
    )
    def test_maps_all_wannier_orbital_symbols(self, parser, symbol, quantum_numbers):
        assert parser._get_quantum_numbers_from_symbol(symbol) == quantum_numbers

    def test_unknown_symbol(self, parser):
        assert parser._get_quantum_numbers_from_symbol('unknown') is None

    def test_calculates_ml_from_orbital_position(self, parser):
        assert parser._calculate_ml_from_position(2, 0) == -2
        assert parser._calculate_ml_from_position(2, 3) == 1

    def test_maps_l_based_and_symbol_based_orbital_states(self, parser):
        assert parser.get_orbitals_state('l=2,mr=1') == [
            {'spin_orbit_state': {'l_quantum_number': 2, 'ml_quantum_number': -1}}
        ]
        assert parser.get_orbitals_state('px;dxy') == [
            {'spin_orbit_state': {'l_quantum_number': 1, 'ml_quantum_number': -1}},
            {'spin_orbit_state': {'l_quantum_number': 2, 'ml_quantum_number': -2}},
        ]
        assert parser.get_orbitals_state(None) is None
