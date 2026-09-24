from tests.parsers.common import SimulationParserPipelineTestSuite


class TestGromacsWaterPipeline(SimulationParserPipelineTestSuite):
    archive_fixture = 'water_archive'
    expected_program_name = 'GROMACS'
