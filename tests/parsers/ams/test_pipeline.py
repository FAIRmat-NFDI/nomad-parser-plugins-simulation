from tests.parsers.common import SimulationParserPipelineTestSuite


class TestAMSParserPipelineSuite(SimulationParserPipelineTestSuite):
    archive_fixture = 'scf_archive'
    expected_program_name = 'AMS'
