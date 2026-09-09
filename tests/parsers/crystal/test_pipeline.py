from tests.parsers.common import SimulationParserPipelineTestSuite


class TestCrystalParserPipelineSuite(SimulationParserPipelineTestSuite):
    archive_fixture = 'single_point_archive'
    expected_program_name = 'Crystal'
