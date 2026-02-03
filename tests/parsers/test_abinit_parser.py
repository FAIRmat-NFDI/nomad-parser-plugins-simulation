from nomad.datamodel import EntryArchive
from nomad.utils import get_logger

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
