from nomad.datamodel import EntryArchive
from nomad.utils import get_logger
from pytest import approx

from nomad_simulation_parsers.parsers.abinit.parser import AbinitParser

LOGGER = get_logger(__name__)


def _parse(mainfile: str) -> EntryArchive:
    parser = AbinitParser()
    archive = EntryArchive()
    parser.parse(mainfile, archive, LOGGER)
    return archive


def test_parse_file():
    _parse('tests/data/abinit/Fe/Fe.out')


def test_parse_file_has_core_sections_and_outputs():
    archive = _parse('tests/data/abinit/Fe/Fe.out')

    assert len(archive.data.model_system or []) > 0
    assert len(archive.data.model_method or []) > 0
    assert len(archive.data.outputs or []) > 0

    output = archive.data.outputs[0]
    assert len(output.total_energies or []) > 0
    assert len(output.total_forces or []) > 0
    assert len(output.electronic_band_structures or []) > 0

    band_structure = output.electronic_band_structures[0]
    assert band_structure.value is not None

    if output.electronic_band_gaps:
        band_gap = output.electronic_band_gaps[0]
        assert band_gap.value is not None

    if output.electronic_dos:
        dos = output.electronic_dos[0]
        assert dos.value is not None
        assert dos.energies is not None
        assert dos.energies.points is not None


def test_convergence_targets_parsing():
    archive = _parse('tests/data/abinit/Fe/Fe.out')

    if archive.workflow2 is not None:
        method = archive.workflow2.method
        if hasattr(method, 'convergence_targets') and method.convergence_targets:
            assert len(method.convergence_targets) > 0
            for target in method.convergence_targets:
                assert target.threshold is not None
                assert target.threshold_type in [
                    'absolute',
                    'maximum',
                    'rms',
                    'relative',
                ]

            target_names = [t.m_def.name for t in method.convergence_targets]
            assert any('Energy' in name for name in target_names) or any(
                'Force' in name for name in target_names
            )


def test_scf_steps_parsing():
    archive = _parse('tests/data/abinit/Fe/Fe.out')

    outputs = archive.data.outputs
    assert outputs is not None
    assert len(outputs) == 2

    first_output_steps = outputs[0].scf_steps
    second_output_steps = outputs[1].scf_steps

    assert first_output_steps is not None
    assert second_output_steps is not None

    assert len(first_output_steps.energies_total) == 16
    assert len(second_output_steps.energies_total) == 30

    assert first_output_steps.energies_total[0].to('hartree').magnitude == approx(-23.45196)
    assert first_output_steps.energies_total[-1].to('hartree').magnitude == approx(
        -24.6617073
    )
    first_delta_last = first_output_steps.delta_energies_total[-1].to('hartree').magnitude
    assert first_delta_last == approx(1.243e-13)

    second_delta_last = (
        second_output_steps.delta_energies_total[-1].to('hartree').magnitude
    )
    assert second_delta_last == approx(0.0)


def test_single_point_workflow_convergence_section():
    archive = _parse('tests/data/abinit/Fe/Fe.out')

    assert archive.workflow2 is not None
    assert archive.workflow2.m_def.name == 'SinglePoint'
    assert archive.workflow2.method is None


def test_geometry_optimization_workflow_convergence_section():
    archive = _parse('tests/data/abinit/H2/H2.out')

    workflow = archive.workflow2
    assert workflow is not None
    assert workflow.m_def.name == 'GeometryOptimization'
    assert workflow.method is not None

    method = workflow.method
    assert method.optimization_method == 'bfgs'

    targets = method.convergence_targets
    assert targets is not None
    assert len(targets) == 2

    targets_by_name = {target.m_def.name: target for target in targets}
    assert set(targets_by_name.keys()) == {
        'EnergyConvergenceTarget',
        'ForceConvergenceTarget',
    }

    energy_target = targets_by_name['EnergyConvergenceTarget']
    force_target = targets_by_name['ForceConvergenceTarget']
    assert energy_target.threshold_type == 'relative'
    assert energy_target.threshold == approx(0.0)
    assert force_target.threshold_type == 'maximum'
    force_threshold = (
        force_target.threshold.to('newton').magnitude
        if hasattr(force_target.threshold, 'to')
        else force_target.threshold
    )
    assert force_threshold == approx(5.0e-4)
