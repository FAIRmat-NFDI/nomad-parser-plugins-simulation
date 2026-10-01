from tests.parsers.common import SimulationParserPipelineTestSuite


class TestFe2Pipeline(SimulationParserPipelineTestSuite):
    archive_fixture = 'fe2_archive'
    expected_program_name = 'GPAW'
