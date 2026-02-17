import numpy as np
import pytest
from nomad.datamodel import EntryArchive
from nomad.utils import get_logger

from nomad_simulation_parsers.parsers.exciting.parser import ExcitingParser

LOGGER = get_logger(__name__)

CASES = {
    'C_minimal': {
        'workflow_name': 'SinglePoint',
        'target_groups': {
            'convergence_targets': {
                'PotentialConvergenceTarget': {
                    'threshold_type': 'rms',
                    'threshold': (4.3597447222071695e-24, 'joule'),
                },
                'EnergyConvergenceTarget': {
                    'threshold_type': 'absolute',
                    'threshold': (4.3597447222071695e-24, 'joule'),
                },
                'ChargeConvergenceTarget': {
                    'threshold_type': 'absolute',
                    'threshold': (1.602176634e-24, 'coulomb'),
                },
            }
        },
        'scf_expectations': {
            'delta_energies_total': {'len': 10, 'last': (4.3213092127361915e-25, 'joule')},
            'delta_potential_rms': {'len': 10, 'last': (1.41636334713761e-26, 'joule')},
            'delta_density_rms': {'len': 10, 'last': (3.77314199483634e-28, 'coulomb')},
            'delta_force_abs': None,
        },
    },
    'GaO_sodium': {
        'workflow_name': 'GeometryOptimization',
        'target_groups': {
            'convergence_targets': {
                'ForceConvergenceTarget': {
                    'threshold_type': 'maximum',
                    'threshold': (4.119361749119496e-09, 'newton'),
                },
            },
            'single_point_convergence_targets': {
                'PotentialConvergenceTarget': {
                    'threshold_type': 'rms',
                    'threshold': (4.3597447222071695e-24, 'joule'),
                },
                'EnergyConvergenceTarget': {
                    'threshold_type': 'absolute',
                    'threshold': (4.35974472220717e-22, 'joule'),
                },
                'ChargeConvergenceTarget': {
                    'threshold_type': 'absolute',
                    'threshold': (1.602176634e-24, 'coulomb'),
                },
                'ForceConvergenceTarget': {
                    'threshold_type': 'absolute',
                    'threshold': (4.1193617491194954e-12, 'newton'),
                },
            },
        },
        'scf_expectations': {
            'delta_energies_total': {'len': 24, 'last': (1.4699010123263138e-23, 'joule')},
            'delta_potential_rms': {'len': 24, 'last': (2.13237730209986e-25, 'joule')},
            'delta_density_rms': {'len': 24, 'last': (6.666224386382819e-27, 'coulomb')},
            'delta_force_abs': {'len': 24, 'last': (1.047809093229533e-13, 'newton')},
        },
    },
}


def _assert_quantity_close(quantity, expected_value: float, unit: str) -> None:
    assert np.isclose(quantity.to(unit).magnitude, expected_value, rtol=1e-12, atol=0.0)


@pytest.fixture(scope='module')
def parser():
    return ExcitingParser()


@pytest.fixture(params=[pytest.param(case, id=case) for case in CASES])
def parsed_archive(request, parser):
    case = request.param
    archive = EntryArchive()
    parser.parse(f'tests/data/exciting/{case}/INFO.OUT', archive, LOGGER)
    return case, archive


def test_parse_file(parsed_archive):
    """Test basic parsing without crashing for the explicit SCF fixtures."""
    _, archive = parsed_archive
    assert archive.workflow2 is not None


def test_convergence_target_thresholds_and_types(parsed_archive):
    case, archive = parsed_archive
    expected = CASES[case]
    method = archive.workflow2.method
    assert method is not None
    assert archive.workflow2.m_def.name == expected['workflow_name']

    for group_name, expected_targets in expected['target_groups'].items():
        group = getattr(method, group_name)
        assert len(group) == len(expected_targets)
        actual_targets = {t.m_def.name: t for t in group}
        assert set(actual_targets.keys()) == set(expected_targets.keys())
        for target_name, target_expectation in expected_targets.items():
            target = actual_targets[target_name]
            assert target.threshold_type == target_expectation['threshold_type']
            threshold_value, threshold_unit = target_expectation['threshold']
            _assert_quantity_close(target.threshold, threshold_value, threshold_unit)


def test_scf_quantities_relevant_for_convergence_targets(parsed_archive):
    case, archive = parsed_archive
    expected = CASES[case]['scf_expectations']
    scf_steps = archive.data.outputs[0].scf_steps

    for quantity_name, quantity_expectation in expected.items():
        value = getattr(scf_steps, quantity_name)
        if quantity_expectation is None:
            assert value is None
            continue
        assert value is not None
        assert len(value) == quantity_expectation['len']
        expected_last_value, expected_last_unit = quantity_expectation['last']
        _assert_quantity_close(value[-1], expected_last_value, expected_last_unit)
