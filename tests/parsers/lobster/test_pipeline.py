from tests.parsers.common import SimulationParserPipelineTestSuite


class TestFePipeline(SimulationParserPipelineTestSuite):
    archive_fixture = 'fe_archive'
    expected_program_name = 'LOBSTER'
