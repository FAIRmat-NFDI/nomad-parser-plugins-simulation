from tests.parsers.common import SimulationParserPipelineTestSuite


class TestH5MDPipeline(SimulationParserPipelineTestSuite):
    archive_fixture = 'h5md_archive'
    expected_program_name = 'OpenMM'
