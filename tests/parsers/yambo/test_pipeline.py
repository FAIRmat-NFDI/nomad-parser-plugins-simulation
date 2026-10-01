from tests.parsers.common import SimulationParserPipelineTestSuite


class TestYamboPipeline(SimulationParserPipelineTestSuite):
    archive_fixture = 'hbn_archive'
    expected_program_name = 'YAMBO'
