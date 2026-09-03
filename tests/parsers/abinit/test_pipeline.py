from tests.parsers.common import SimulationParserPipelineTestSuite


class TestAbinitParserPipelineSuite(SimulationParserPipelineTestSuite):
    archive_fixture = 'fe_archive'
    expected_program_name = 'ABINIT'
