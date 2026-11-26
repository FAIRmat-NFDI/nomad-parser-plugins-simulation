import numpy as np
import pytest
from nomad.units import ureg

from nomad_simulation_parsers.parsers.crystal.parser import CrystalOutputParser


@pytest.mark.unit
def test_maps_atomic_numbers_and_labels():
    parser = CrystalOutputParser()
    source = {
        'labels_positions': np.array(
            [['1', 'T', '238', 'Sr', '0.0', '0.0', '0.0']], dtype=str
        )
    }

    atoms = parser.get_atoms(source)

    assert atoms == [{'label': 'Sr', 'number': 38}]


@pytest.mark.unit
def test_maps_outputs_with_explicit_units():
    parser = CrystalOutputParser()
    parser._data = {'dimensionality': 3}
    source = {
        'lattice_parameters': np.eye(3) * 5.0,
        'labels_positions': np.array(
            [['1', 'T', '14', 'Si', '0.0', '0.0', '0.0']], dtype=str
        ),
        'energy_total': -10.0 * ureg.hartree,
        'forces': np.array([[1.0, 2.0, 3.0, 4.0, 5.0]]),
    }

    outputs = parser.get_outputs(source)

    assert outputs[0]['energy'].to('hartree').magnitude == pytest.approx(-10.0)
    np.testing.assert_allclose(
        outputs[0]['forces'].to('hartree / bohr').magnitude, [[3.0, 4.0, 5.0]]
    )


@pytest.mark.unit
@pytest.mark.parametrize(
    'source, expected',
    [
        ({'exchange': 'BECKE', 'correlation': 'PZ'}, ['GGA_X_B88', 'LDA_C_PZ']),
        ({'correlation': 'PZ', 'exchange': 'BECKE'}, ['GGA_X_B88', 'LDA_C_PZ']),
        ({'exchange_correlation': 'PBE0'}, ['HYB_GGA_XC_PBEH']),
        ({'exchange_correlation': 'PBEXC'}, ['GGA_C_PBE', 'GGA_X_PBE']),
        (
            {'exchange': 'PBE', 'exchange_correlation': 'PBE'},
            ['GGA_X_PBE', 'GGA_C_PBE'],
        ),
        ({'exchange': 'XYZ'}, []),
        ({}, []),
    ],
    ids=[
        'exchange_then_correlation',
        'correlation_key_first_in_source',
        'combined_functional',
        'multi_component_entry',
        'duplicate_across_keys',
        'unknown_functional',
        'empty_source',
    ],
)
def test_xc_functionals_order_is_deterministic(source, expected):
    """Components follow the fixed key order, deduplicated, never set order."""
    parser = CrystalOutputParser()

    assert [xc['name'] for xc in parser.get_xc_functionals(source)] == expected
