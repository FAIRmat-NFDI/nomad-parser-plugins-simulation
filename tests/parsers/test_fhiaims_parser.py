from nomad.datamodel import EntryArchive
from nomad.utils import get_logger
from pytest import approx, mark

from nomad_simulation_parsers.parsers.fhiaims.parser import FHIAimsParser

LOGGER = get_logger(__name__)


def test_parse_file():
    parser = FHIAimsParser()
    archive = EntryArchive()
    parser.parse('tests/data/fhiaims/Si_geomopt/out.out', archive, LOGGER)


def test_workflow_convergence_targets():
    parser = FHIAimsParser()
    archive = EntryArchive()
    parser.parse('tests/data/fhiaims/Si_geomopt/out.out', archive, LOGGER)

    workflow = archive.workflow2
    assert workflow is not None
    assert workflow.m_def.name == 'GeometryOptimization'
    assert workflow.method is not None
    assert workflow.method.optimization_method is not None
    assert workflow.method.optimization_method.strip() == 'Modified BFGS'

    targets = workflow.method.convergence_targets
    assert targets is not None
    assert len(targets) == 1
    force_target = targets[0]
    assert force_target.m_def.name == 'ForceConvergenceTarget'
    assert force_target.threshold_type == 'maximum'
    assert force_target.threshold.to('eV/angstrom').magnitude == approx(0.01)

    sp_targets = workflow.method.single_point_convergence_targets
    assert sp_targets is not None
    assert len(sp_targets) == 1
    energy_target = sp_targets[0]
    assert energy_target.m_def.name == 'EnergyConvergenceTarget'
    assert energy_target.threshold_type == 'absolute'
    assert energy_target.threshold.to('eV').magnitude == approx(1.0e-6)


def test_scf_steps_quantities():
    parser = FHIAimsParser()
    archive = EntryArchive()
    parser.parse('tests/data/fhiaims/Si_geomopt/out.out', archive, LOGGER)

    outputs = archive.data.outputs
    assert outputs is not None
    assert len(outputs) == 5

    expected_n_scf = [11, 7, 8, 6, 7]
    for output, n_scf in zip(outputs, expected_n_scf):
        scf_steps = output.scf_steps
        assert scf_steps is not None
        assert len(scf_steps.delta_energies_total) == n_scf
        assert len(scf_steps.delta_density_rms) == n_scf
        assert len(scf_steps.durations) == n_scf

    # Explicitly check last SCF-step deltas in first geometry-optimization step
    first_steps = outputs[1].scf_steps
    assert first_steps.delta_energies_total[-1].to('eV').magnitude == approx(7.477e-09)
    assert first_steps.delta_density_rms[-1].to('coulomb').magnitude == approx(
        6.375e-08 * 1.602176634e-19
    )


@mark.parametrize(
    'k_offset_line,expected_offset',
    [
        (None, [0.0, 0.0, 0.0]),
        ('  k_offset                           0.5 0.25 0.0\n', [0.5, 0.25, 0.0]),
    ],
    ids=['default_offset', 'explicit_offset'],
)
def test_k_mesh(tmp_path, k_offset_line, expected_offset):
    """Test k-mesh parsing with and without explicit k_offset."""
    source_path = 'tests/data/fhiaims/Si_geomopt/out.out'

    if k_offset_line is None:
        # Use original file (no k_offset)
        test_file = source_path
    else:
        # Inject k_offset line after k_grid
        with open(source_path, encoding='utf-8') as f:
            content = f.read()
        modified_content = content.replace(
            '  Found k-point grid:         8         8         8\n',
            '  Found k-point grid:         8         8         8\n' + k_offset_line,
            1,
        )
        assert modified_content != content, 'Failed to inject k_offset line'
        test_file = tmp_path / 'out_with_k_offset.out'
        test_file.write_text(modified_content, encoding='utf-8')

    parser = FHIAimsParser()
    archive = EntryArchive()
    parser.parse(str(test_file), archive, LOGGER)

    # Check DFT section exists
    assert archive.data.model_method is not None
    assert len(archive.data.model_method) == 1
    dft = archive.data.model_method[0]
    assert dft.m_def.name == 'DFT'

    # Check NumericalSettings/KSpace exists
    assert dft.numerical_settings is not None
    assert len(dft.numerical_settings) == 1
    k_space = dft.numerical_settings[0]
    assert k_space.m_def.name == 'KSpace'

    # Check KSpace.k_mesh exists
    assert k_space.k_mesh is not None
    assert len(k_space.k_mesh) == 1
    k_mesh = k_space.k_mesh[0]
    assert k_mesh.m_def.name == 'KMesh'

<<<<<<< HEAD
    # Check k-grid values
    assert k_mesh.grid is not None
    assert list(k_mesh.grid) == [8, 8, 8]

    # Check k_offset values (default or explicit)
    assert k_mesh.offset is not None
    assert list(k_mesh.offset) == approx(expected_offset)
=======
    if output.electronic_band_gaps:
        assert output.electronic_band_gaps[0].value is not None

    if output.electronic_eigenvalues:
        assert output.electronic_band_structures is not None
        assert output.electronic_band_structures[0].value is not None
>>>>>>> 2d1c3b0 (patch outputs and octopus reader)
