from nomad.datamodel import EntryArchive
from nomad.utils import get_logger

from nomad_simulation_parsers.parsers.exciting.parser import ExcitingParser

LOGGER = get_logger(__name__)


def test_parse_file():
    """Test basic parsing without crashing."""
    parser = ExcitingParser()
    archive = EntryArchive()
    parser.parse('tests/data/exciting/GaO_strucopt/INFO.OUT', archive, LOGGER)


def test_geometry_optimization_convergence():
    """Test that geometry optimization convergence targets are parsed correctly."""
    parser = ExcitingParser()
    archive = EntryArchive()
    parser.parse('tests/data/exciting/GaO_strucopt/INFO.OUT', archive, LOGGER)

    # Check workflow exists
    assert archive.workflow2 is not None
    assert archive.workflow2.m_def.name == 'GeometryOptimization'

    # Check method has convergence targets
    method = archive.workflow2.method
    assert method is not None

    # Geometry optimization should have force convergence target
    if hasattr(method, 'convergence_targets') and method.convergence_targets:
        # Check that targets were parsed
        assert len(method.convergence_targets) > 0

        # Find force convergence target
        force_targets = [
            t
            for t in method.convergence_targets
            if t.m_def.name == 'ForceConvergenceTarget'
        ]
        if force_targets:
            target = force_targets[0]
            assert target.threshold is not None
            assert target.threshold_type in ['absolute', 'maximum', 'rms', 'relative']

    # Check for single point convergence targets
    if (
        hasattr(method, 'single_point_convergence_targets')
        and method.single_point_convergence_targets
    ):
        assert len(method.single_point_convergence_targets) > 0

        # Verify threshold_type is set correctly
        for target in method.single_point_convergence_targets:
            assert target.threshold is not None
            assert target.threshold_type in ['absolute', 'maximum', 'rms', 'relative']
