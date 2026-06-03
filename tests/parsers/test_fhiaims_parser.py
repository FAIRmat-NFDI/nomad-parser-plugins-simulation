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
    # Filter for KSpace (may also contain SelfConsistency criteria)
    k_spaces = [ns for ns in dft.numerical_settings if ns.m_def.name == 'KSpace']
    assert len(k_spaces) == 1
    k_space = k_spaces[0]

    # Check KSpace.k_mesh exists
    assert k_space.k_mesh is not None
    assert len(k_space.k_mesh) == 1
    k_mesh = k_space.k_mesh[0]
    assert k_mesh.m_def.name == 'KMesh'

    # Check k-grid values
    assert k_mesh.grid is not None
    assert list(k_mesh.grid) == [8, 8, 8]

    # Check k_offset values (default or explicit)
    assert k_mesh.offset is not None
    assert list(k_mesh.offset) == approx(expected_offset)


def test_scf_convergence_criteria():
    """Test extraction of SCF convergence criteria."""
    parser = FHIAimsParser()
    archive = EntryArchive()
    parser.parse('tests/data/fhiaims/Si_geomopt/out.out', archive, LOGGER)

    # Check DFT section exists
    assert archive.data.model_method is not None
    assert len(archive.data.model_method) == 1
    dft = archive.data.model_method[0]
    assert dft.m_def.name == 'DFT'

    # Check NumericalSettings contains SelfConsistency sections
    assert dft.numerical_settings is not None
    # Should have at least: KSpace + 3 SelfConsistency (energy, density, eigenvalues)
    assert len(dft.numerical_settings) >= 4

    # Filter for SelfConsistency sections
    scf_criteria = [
        ns for ns in dft.numerical_settings if ns.m_def.name == 'SelfConsistency'
    ]
    assert len(scf_criteria) == 3

    # Check energy convergence criterion (distinguished by name)
    energy_criterion = next(
        (sc for sc in scf_criteria if sc.name == 'total_energy_change'),
        None,
    )
    assert energy_criterion is not None
    assert energy_criterion.threshold_change.magnitude == approx(1.0e-6)
    assert str(energy_criterion.threshold_change.units) == 'electron_volt'

    # Check density convergence criterion (distinguished by name)
    density_criterion = next(
        (sc for sc in scf_criteria if sc.name == 'charge_density_change'),
        None,
    )
    assert density_criterion is not None
    assert density_criterion.threshold_change == approx(1.0e-5)

    # Check eigenvalues convergence criterion (distinguished by name)
    eigenvalues_criterion = next(
        (sc for sc in scf_criteria if sc.name == 'sum_eigenvalues_change'),
        None,
    )
    assert eigenvalues_criterion is not None
    assert eigenvalues_criterion.threshold_change.magnitude == approx(1.0e-3)
    assert str(eigenvalues_criterion.threshold_change.units) == 'electron_volt'


def test_total_energies():
    """Test that total energies are extracted and mapped correctly."""
    parser = FHIAimsParser()
    archive = EntryArchive()
    parser.parse('tests/data/fhiaims/Si_geomopt/out.out', archive, LOGGER)

    # Check outputs exist
    outputs = archive.data.outputs
    assert outputs is not None
    assert len(outputs) == 5  # Geometry optimization has 5 steps

    # Check each output has total_energies
    for i, output in enumerate(outputs):
        assert output.total_energies is not None, f'Output {i} has no total_energies'
        assert len(output.total_energies) == 1, (
            f'Output {i} should have exactly 1 total_energy'
        )

        total_energy = output.total_energies[0]

        # Check that energy value exists and is a Quantity
        assert total_energy.value is not None, f'Output {i} total_energy has no value'
        assert hasattr(total_energy.value, 'magnitude'), (
            f'Output {i} energy value is not a Quantity'
        )

        # Energy should be negative for this Si system
        energy_ev = total_energy.value.to('eV').magnitude
        assert energy_ev < 0, f'Output {i} energy should be negative, got {energy_ev}'

        # Check energy is in reasonable range for Si (~-15696 eV for 2 atoms)
        assert -16000 < energy_ev < -15600, (
            f'Output {i} energy {energy_ev} eV out of expected range'
        )

        # Check that contributions/components exist
        assert total_energy.contributions is not None, (
            f'Output {i} has no energy contributions'
        )
        assert len(total_energy.contributions) > 0, (
            f'Output {i} should have energy contributions'
        )

    # Explicitly check first step energy value (from test data analysis)
    # First output should have "Total energy uncorrected" from compact form
    # which is -15696.1246183848 eV
    first_energy = outputs[0].total_energies[0].value.to('eV').magnitude
    assert first_energy == approx(-15696.1246183848, abs=0.01)

    # Check that energy contributions have names and values
    first_contributions = outputs[0].total_energies[0].contributions
    has_eigenvalues = any(
        'eigenvalues' in contrib.name.lower()
        for contrib in first_contributions
        if contrib.name
    )
    assert has_eigenvalues, 'Should have eigenvalues contribution'

    # All contributions should have a value
    for contrib in first_contributions:
        assert contrib.value is not None, f'Contribution {contrib.name} has no value'
