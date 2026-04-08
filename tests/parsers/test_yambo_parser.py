import numpy as np
from pytest import approx

from nomad_simulation_parsers.parsers.yambo.parser import YamboMainfileParser


def test_get_outputs_sets_highest_occupied_from_valence_conduction():
    parser = YamboMainfileParser()

    energies_occupations = {
        'eigenenergies': {
            'kpoints': np.array([[0.0, 0.0, 0.0]]),
            'energies': np.array([1.0, 2.0]),
        },
        'valence_conduction': [1.2, 2.3],
    }

    outputs = parser.get_outputs(energies_occupations, [], {})

    assert outputs
    assert outputs[0].get('highest_occupied') is not None
    assert outputs[0]['highest_occupied'].to('eV').magnitude == approx(1.2)


def test_get_outputs_emits_eigenvalues_payload():
    parser = YamboMainfileParser()

    energies_occupations = {
        'eigenenergies': {
            'kpoints': np.array([[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]]),
            'energies': np.array([1.0, 2.0, 1.5, 2.5]),
        }
    }

    outputs = parser.get_outputs(energies_occupations, [], {})

    assert outputs
    assert outputs[0].get('eigenvalues') is not None
    assert len(outputs[0]['eigenvalues']) == 1
    assert outputs[0]['eigenvalues'][0].get('energies') is not None
