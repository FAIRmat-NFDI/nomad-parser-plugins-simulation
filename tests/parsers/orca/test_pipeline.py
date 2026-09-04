from tests.parsers.common import SimulationParserPipelineTestSuite


class TestRIMP2WaterPipeline(SimulationParserPipelineTestSuite):
    archive_fixture = 'ri_mp2_water_archive'
    expected_program_name = 'ORCA'
