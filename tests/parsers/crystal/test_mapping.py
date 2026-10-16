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
