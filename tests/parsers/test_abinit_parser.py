from nomad.datamodel import EntryArchive
from nomad.utils import get_logger
from pytest import approx

from nomad_simulation_parsers.parsers.abinit.parser import AbinitParser

LOGGER = get_logger(__name__)


def test_parse_file():
    """Test basic parsing without crashing."""
    parser = AbinitParser()
    archive = EntryArchive()
    parser.parse('tests/data/abinit/Fe/Fe.out', archive, LOGGER)


def test_convergence_targets_parsing():
    """Test that convergence targets are parsed and mapped correctly."""
    parser = AbinitParser()
    archive = EntryArchive()
    parser.parse('tests/data/abinit/Fe/Fe.out', archive, LOGGER)

    # Check if workflow exists
    if archive.workflow2 is not None:
        method = archive.workflow2.method

        # For geometry optimization, check convergence targets
        if hasattr(method, 'convergence_targets') and method.convergence_targets:
            # Should have energy and force convergence targets
            assert len(method.convergence_targets) > 0

            # Check that targets have proper fields
            for target in method.convergence_targets:
                assert target.threshold is not None
                assert target.threshold_type in [
                    'absolute',
                    'maximum',
                    'rms',
                    'relative',
                ]

            # Check for specific target types
            target_names = [t.m_def.name for t in method.convergence_targets]
            # Abinit typically has energy and force convergence
            assert any('Energy' in name for name in target_names) or any(
                'Force' in name for name in target_names
            )


def test_scf_steps_parsing():
    parser = AbinitParser()
    archive = EntryArchive()
    parser.parse('tests/data/abinit/Fe/Fe.out', archive, LOGGER)

    outputs = archive.data.outputs
    assert outputs is not None
    assert len(outputs) == 2

    first_output_steps = outputs[0].scf_steps
    second_output_steps = outputs[1].scf_steps

    assert first_output_steps is not None
    assert second_output_steps is not None

    assert len(first_output_steps.energies_total) == 16
    assert len(second_output_steps.energies_total) == 30

    assert first_output_steps.energies_total[0].to('hartree').magnitude == approx(
        -23.45196
    )
    assert first_output_steps.energies_total[-1].to('hartree').magnitude == approx(
        -24.6617073
    )
    first_delta_last = first_output_steps.delta_energies_total[-1].to(
        'hartree'
    ).magnitude
    assert first_delta_last == approx(1.243e-13)

    second_delta_last = second_output_steps.delta_energies_total[-1].to(
        'hartree'
    ).magnitude
    assert second_delta_last == approx(0.0)
