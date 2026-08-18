"""Tests for common parser utilities in parsers/utils/general.py"""

import numpy as np
import pytest
from nomad.units import ureg

from nomad_simulation_parsers.parsers.utils.general import (
    OCCUPATION_THRESHOLD,
    calculate_band_gap_from_occupations,
    create_mapping_table,
)


class TestOccupationThreshold:
    """Test the common OCCUPATION_THRESHOLD constant."""

    def test_threshold_value(self):
        """Verify occupation threshold is 0.5 as per DFT convention."""
        assert OCCUPATION_THRESHOLD == 0.5


class _Quantity:
    def __init__(self, name, sub_parser=None):
        self.name = name
        self.sub_parser = sub_parser


class _FileParser:
    def __init__(self, quantities):
        self.quantities = quantities


class _Path:
    def __init__(self, path):
        self.absolute_path = path


class _Mapper:
    def __init__(self, paths, function_name=None):
        self.source = None
        self.function_args = [_Path(path) for path in paths]
        self.mappers = []
        self.function_name = function_name


class _MappingFunctions:
    def transform(self, source):
        return source.get('unused')


class _MappingFunctionsWithMultipleAccessStyles:
    def transform(self, data):
        return data.get('get_value'), data['dict_value'], data.object_value


class _MappingFunctionsWithHelper:
    def transform(self, data):
        return self.extract(data)

    def extract(self, data):
        values = data.get('values', {})
        return values['indirect_value']


def test_create_mapping_table_compares_file_quantities_with_mapper_paths():
    file_parser = _FileParser(
        [
            _Quantity(
                'eigenvalues_occupancies',
                _FileParser(
                    [_Quantity('eigenvalues'), _Quantity('occupancies')]
                ),
            ),
            _Quantity('n_states'),
            _Quantity('unused'),
        ]
    )
    archive_parser = _FileParser([])
    archive_parser.mapper = _Mapper(['eigenvalues_occupancies', 'n_states'])

    assert create_mapping_table(file_parser, archive_parser) == [
        {
            'quantity': 'eigenvalues_occupancies',
            'mapped': True,
            'mapping': ['eigenvalues_occupancies'],
        },
        {
            'quantity': 'eigenvalues_occupancies.eigenvalues',
            'mapped': True,
            'mapping': ['eigenvalues_occupancies'],
        },
        {
            'quantity': 'eigenvalues_occupancies.occupancies',
            'mapped': True,
            'mapping': ['eigenvalues_occupancies'],
        },
        {'quantity': 'n_states', 'mapped': True, 'mapping': ['n_states']},
        {'quantity': 'unused', 'mapped': False, 'mapping': []},
    ]


def test_create_mapping_table_inspects_mapping_functions():
    file_parser = _FileParser([_Quantity('unused')])
    archive_parser = _FileParser([])
    archive_parser.mapper = _Mapper([], function_name='transform')

    assert create_mapping_table(
        file_parser, archive_parser, function_objects=_MappingFunctions()
    ) == [
        {
            'quantity': 'unused',
            'mapped': True,
            'mapping': ['transform(unused)'],
        }
    ]


def test_create_mapping_table_inspects_mapping_function_access_styles():
    file_parser = _FileParser(
        [_Quantity('get_value'), _Quantity('dict_value'), _Quantity('object_value')]
    )
    archive_parser = _FileParser([])
    archive_parser.mapper = _Mapper([], function_name='transform')

    assert create_mapping_table(
        file_parser,
        archive_parser,
        function_objects=_MappingFunctionsWithMultipleAccessStyles(),
    ) == [
        {
            'quantity': 'get_value',
            'mapped': True,
            'mapping': ['transform(get_value)'],
        },
        {
            'quantity': 'dict_value',
            'mapped': True,
            'mapping': ['transform(dict_value)'],
        },
        {
            'quantity': 'object_value',
            'mapped': True,
            'mapping': ['transform(object_value)'],
        },
    ]


def test_create_mapping_table_inspects_helper_function_accesses():
    file_parser = _FileParser([_Quantity('indirect_value')])
    archive_parser = _FileParser([])
    archive_parser.mapper = _Mapper([], function_name='transform')

    assert create_mapping_table(
        file_parser,
        archive_parser,
        function_objects=_MappingFunctionsWithHelper(),
    ) == [
        {
            'quantity': 'indirect_value',
            'mapped': True,
            'mapping': ['transform(indirect_value)'],
        }
    ]


class TestCalculateBandGapFromOccupations:
    """Test band gap calculation utility."""

    def test_simple_gap(self):
        """Test basic band gap calculation for semiconductor."""
        eigenvalues = np.array([-5.0, -4.0, -3.0, 3.0, 4.0])
        occupations = np.array([1.0, 1.0, 1.0, 0.0, 0.0])

        result = calculate_band_gap_from_occupations(eigenvalues, occupations)

        assert result is not None
        assert 'value' in result
        assert result['value'] == pytest.approx(6.0)  # Gap from -3.0 to 3.0
        assert 'spin_channel' not in result

    def test_partial_occupation_threshold_gap(self):
        """Partial occupations are split at the threshold, yielding a finite gap.

        The threshold-based algorithm classifies states as occupied or
        unoccupied at 0.5, so a metal with fractional occupations reports the
        energy distance between the states straddling the threshold rather
        than a zero gap.
        """
        eigenvalues = np.array([-5.0, -4.0, -3.0, -2.0, -1.0])
        occupations = np.array([1.0, 1.0, 0.8, 0.2, 0.0])

        result = calculate_band_gap_from_occupations(eigenvalues, occupations)

        assert result is not None
        # Occupied: [-5, -4, -3] (occ >= 0.5)
        # Unoccupied: [-2, -1] (occ < 0.5)
        # Gap = -2 - (-3) = 1.0
        assert result['value'] == pytest.approx(1.0)

    def test_negative_gap_forced_to_zero(self):
        """Test that negative gaps (metals with partial occupations) are zeroed."""
        eigenvalues = np.array([-5.0, -4.0, -3.0, -2.0])
        occupations = np.array([1.0, 1.0, 0.6, 0.4])

        result = calculate_band_gap_from_occupations(eigenvalues, occupations)

        assert result is not None
        # Occupied: [-5, -4, -3] (occ >= 0.5)
        # Unoccupied: [-2] (occ < 0.5)
        # Gap would be -2 - (-3) = 1.0 (still positive)
        assert result['value'] >= 0.0

    def test_with_spin_channel(self):
        """Test spin channel labeling."""
        eigenvalues = np.array([-5.0, -4.0, 3.0, 4.0])
        occupations = np.array([1.0, 1.0, 0.0, 0.0])

        result = calculate_band_gap_from_occupations(
            eigenvalues, occupations, spin_channel=1
        )

        assert result is not None
        assert result['value'] == pytest.approx(7.0)  # Gap from -4.0 to 3.0
        assert result['spin_channel'] == 1

    def test_with_pint_units(self):
        """Test handling of pint quantities (eigenvalues with units)."""
        eigenvalues = np.array([-5.0, -4.0, 3.0, 4.0]) * ureg.eV
        occupations = np.array([1.0, 1.0, 0.0, 0.0])

        result = calculate_band_gap_from_occupations(eigenvalues, occupations)

        assert result is not None
        assert hasattr(result['value'], 'magnitude')
        assert result['value'].magnitude == pytest.approx(7.0)
        assert result['value'].units == ureg.eV

    def test_with_energy_units_parameter(self):
        """Test applying units via energy_units parameter."""
        eigenvalues = np.array([-5.0, -4.0, 3.0, 4.0])  # Unitless
        occupations = np.array([1.0, 1.0, 0.0, 0.0])

        result = calculate_band_gap_from_occupations(
            eigenvalues, occupations, energy_units=ureg.hartree
        )

        assert result is not None
        assert hasattr(result['value'], 'magnitude')
        assert result['value'].magnitude == pytest.approx(7.0)
        assert result['value'].units == ureg.hartree

    def test_custom_occupation_threshold(self):
        """Test using non-default occupation threshold."""
        eigenvalues = np.array([-5.0, -4.0, -3.0, 3.0, 4.0])
        occupations = np.array([1.0, 0.7, 0.6, 0.4, 0.0])

        # Default threshold (0.5): occupied=[-5,-4,-3], unoccupied=[3,4], gap=6.0
        result_default = calculate_band_gap_from_occupations(eigenvalues, occupations)
        assert result_default['value'] == pytest.approx(6.0)

        # Custom threshold (0.65): occupied=[-5,-4], unoccupied=[-3,3,4], gap=1.0
        result_custom = calculate_band_gap_from_occupations(
            eigenvalues, occupations, occupation_threshold=0.65
        )
        assert result_custom['value'] == pytest.approx(1.0)

    def test_none_inputs(self):
        """Test that None inputs return None."""
        assert calculate_band_gap_from_occupations(None, None) is None
        assert calculate_band_gap_from_occupations(None, np.array([1.0, 0.0])) is None
        assert calculate_band_gap_from_occupations(np.array([1.0, 0.0]), None) is None

    def test_empty_arrays(self):
        """Test that empty arrays return None."""
        assert calculate_band_gap_from_occupations(np.array([]), np.array([])) is None

    def test_shape_mismatch(self):
        """Test that mismatched shapes return None."""
        eigenvalues = np.array([1.0, 2.0, 3.0])
        occupations = np.array([1.0, 0.0])

        result = calculate_band_gap_from_occupations(eigenvalues, occupations)
        assert result is None

    def test_all_occupied(self):
        """Test case where all states are occupied (no gap)."""
        eigenvalues = np.array([-5.0, -4.0, -3.0])
        occupations = np.array([1.0, 1.0, 1.0])

        result = calculate_band_gap_from_occupations(eigenvalues, occupations)
        assert result is None  # No unoccupied states

    def test_all_unoccupied(self):
        """Test case where all states are unoccupied (no gap)."""
        eigenvalues = np.array([-5.0, -4.0, -3.0])
        occupations = np.array([0.0, 0.0, 0.0])

        result = calculate_band_gap_from_occupations(eigenvalues, occupations)
        assert result is None  # No occupied states

    def test_multidimensional_arrays(self):
        """Test handling of 2D arrays (k-points x bands)."""
        # Multiple k-points, flattened energy/occupation arrays
        eigenvalues = np.array(
            [
                [-5.0, -4.0, 3.0, 4.0],  # k-point 1
                [-5.1, -3.9, 3.1, 4.1],  # k-point 2
            ]
        ).flatten()
        occupations = np.array(
            [
                [1.0, 1.0, 0.0, 0.0],
                [1.0, 1.0, 0.0, 0.0],
            ]
        ).flatten()

        result = calculate_band_gap_from_occupations(eigenvalues, occupations)

        assert result is not None
        # VBM should be max(-4.0, -3.9) = -3.9
        # CBM should be min(3.0, 3.1) = 3.0
        # Gap = 3.0 - (-3.9) = 6.9
        assert result['value'] == pytest.approx(6.9)

    def test_spin_polarized_data(self):
        """Test typical spin-polarized eigenvalues (occupation 0-2 range)."""
        # Spin-up channel
        eigenvalues_up = np.array([-5.0, -4.0, 3.0, 4.0])
        occupations_up = np.array([2.0, 2.0, 0.0, 0.0])  # Fully occupied/empty

        result_up = calculate_band_gap_from_occupations(
            eigenvalues_up, occupations_up, spin_channel=0
        )

        assert result_up is not None
        assert result_up['value'] == pytest.approx(7.0)
        assert result_up['spin_channel'] == 0

        # Spin-down channel with different gap
        eigenvalues_down = np.array([-4.5, -3.5, 3.5, 4.5])
        occupations_down = np.array([2.0, 2.0, 0.0, 0.0])

        result_down = calculate_band_gap_from_occupations(
            eigenvalues_down, occupations_down, spin_channel=1
        )

        assert result_down is not None
        assert result_down['value'] == pytest.approx(7.0)
        assert result_down['spin_channel'] == 1

    def test_fractional_occupations(self):
        """Test with fractional occupations (metallic or finite temperature)."""
        eigenvalues = np.array([-5.0, -4.0, -3.0, -2.0, 3.0])
        occupations = np.array([1.0, 1.0, 0.9, 0.1, 0.0])

        result = calculate_band_gap_from_occupations(eigenvalues, occupations)

        assert result is not None
        # Occupied (>= 0.5): [-5, -4, -3]
        # Unoccupied (< 0.5): [-2, 3]
        # Gap = -2 - (-3) = 1.0
        assert result['value'] == pytest.approx(1.0)
