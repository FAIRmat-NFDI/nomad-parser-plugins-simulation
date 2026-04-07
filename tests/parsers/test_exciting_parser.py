import tempfile
import zipfile
from pathlib import Path

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
                    'threshold': (1e-05, 'coulomb'),
                },
            }
        },
        'scf_expectations': {
            'delta_energies_total': {
                'len': 10,
                'last': (4.3213092127361915e-25, 'joule'),
            },
            'delta_potential_rms': {'len': 10, 'last': (1.41636334713761e-26, 'joule')},
            'delta_density_rms': {'len': 10, 'last': (2.35501e-09, 'coulomb')},
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
                    'threshold': (1e-05, 'coulomb'),
                },
                'ForceConvergenceTarget': {
                    'threshold_type': 'absolute',
                    'threshold': (4.1193617491194954e-12, 'newton'),
                },
            },
        },
        'scf_expectations': {
            'delta_energies_total': {
                'len': 24,
                'last': (1.4699010123263138e-23, 'joule'),
            },
            'delta_potential_rms': {'len': 24, 'last': (2.13237730209986e-25, 'joule')},
            'delta_density_rms': {
                'len': 24,
                'last': (4.16073e-08, 'coulomb'),
            },
            'delta_force_abs': {'len': 24, 'last': (1.047809093229533e-13, 'newton')},
        },
    },
}


def _assert_quantity_close(quantity, expected_value: float, unit: str) -> None:
    # Handle both plain floats and Pint Quantities (flexible_unit behavior)
    magnitude = quantity.to(unit).magnitude if hasattr(quantity, 'to') else quantity
    assert np.isclose(magnitude, expected_value, rtol=1e-12, atol=0.0)


@pytest.fixture(scope='module')
def parser():
    return ExcitingParser()


@pytest.fixture(params=[pytest.param(case, id=case) for case in CASES])
def parsed_archive(request, parser):
    case = request.param
    archive = EntryArchive()
    mainfile = (
        Path(__file__).resolve().parents[1] / 'data' / 'exciting' / case / 'INFO.OUT'
    )
    parser.parse(str(mainfile), archive, LOGGER)
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


def test_electronic_outputs_mapping(parsed_archive):
    case, archive = parsed_archive
    if case != 'C_minimal':
        pytest.skip('electronic output assertions use C_minimal fixture')

    outputs = archive.data.outputs
    assert outputs is not None
    assert len(outputs) > 0

    output = outputs[0]

    if output.electronic_dos:
        sec_dos = output.electronic_dos[0]
        assert sec_dos.value is not None
        assert sec_dos.energies is not None
        assert sec_dos.energies.points is not None

    if output.electronic_band_structures:
        sec_band_structure = output.electronic_band_structures[0]
        assert sec_band_structure.value is not None

    if output.electronic_band_gaps:
        sec_gap = output.electronic_band_gaps[0]
        assert sec_gap.value is not None


def test_system_fundamental_quantities_mapping(parsed_archive):
    """System gate: parser should populate core model_system quantities used by normalizer."""
    _, archive = parsed_archive

    simulation = archive.data
    assert simulation is not None
    assert simulation.model_system is not None
    assert len(simulation.model_system) > 0

    representative = next(
        (s for s in simulation.model_system if getattr(s, 'is_representative', False)),
        simulation.model_system[0],
    )
    assert representative.positions is not None
    assert representative.lattice_vectors is not None
    assert representative.periodic_boundary_conditions is not None

    if representative.particle_states:
        assert all(
            getattr(state, 'chemical_symbol', None) is not None
            for state in representative.particle_states
        )


def test_outputs_contract_for_normalizer(parsed_archive):
    """Outputs gate: mapped outputs should include normalizer-required payloads when present."""
    _, archive = parsed_archive

    outputs = archive.data.outputs
    assert outputs is not None
    assert len(outputs) > 0
    output = outputs[0]

    # Non-electronic core output used by normalizer pipeline.
    assert output.total_energies or output.scf_steps is not None

    # Electronic outputs: when present, require payload completeness expected by
    # results normalizer compatibility mapping.
    if output.electronic_dos:
        dos = output.electronic_dos[0]
        assert dos.value is not None
        assert dos.energies is not None
        assert dos.energies.points is not None

    if output.electronic_band_structures:
        bs = output.electronic_band_structures[0]
        assert bs.value is not None
        assert bs.k_path is not None
        assert getattr(bs.k_path, 'points', None) is not None


def test_root_test_data_exciting_zip_populates_reference_energy_fields(parser):
    """Parser scope: exciting zip should expose reference-energy fields on outputs."""
    root_dir = Path(__file__).resolve().parents[3]
    zip_path = root_dir / 'test_data' / 'Si_gw-exciting.zip'
    if not zip_path.is_file():
        pytest.skip('Si_gw-exciting.zip fixture not available in repository root test_data.')

    archive = EntryArchive()
    with tempfile.TemporaryDirectory() as tmpdir:
        with zipfile.ZipFile(zip_path) as zf:
            zf.extractall(tmpdir)
        extracted = Path(tmpdir)
        mainfiles = sorted(extracted.rglob('GW_INFO.OUT')) + sorted(
            extracted.rglob('INFO.OUT')
        )
        assert mainfiles
        parser.parse(str(mainfiles[0]), archive, LOGGER)

    outputs = archive.data.outputs
    assert outputs
    output = outputs[0]

    # Legacy-equivalent migration behavior: parser should provide a reference
    # energy for BS and DOS normalization paths.
    assert output.electronic_band_structures
    bs = output.electronic_band_structures[0]
    assert bs.highest_occupied is not None

    assert output.electronic_dos
    dos = output.electronic_dos[0]
    assert dos.energies is not None
    assert dos.energies.points is not None
    assert dos.energies_origin is not None
