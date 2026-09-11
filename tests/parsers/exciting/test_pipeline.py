from tests.parsers.common import SimulationParserPipelineTestSuite


class TestCMinimalPipeline(SimulationParserPipelineTestSuite):
    archive_fixture = 'c_minimal_archive'
    expected_program_name = 'exciting'
