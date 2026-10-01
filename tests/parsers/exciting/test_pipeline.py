import pytest

from tests.parsers.common import SimulationParserPipelineTestSuite


@pytest.mark.large_fixture
class TestCMinimalPipeline(SimulationParserPipelineTestSuite):
    archive_fixture = 'c_minimal_archive'
    expected_program_name = 'exciting'
