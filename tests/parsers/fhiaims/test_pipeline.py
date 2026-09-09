from tests.parsers.common import SimulationParserPipelineTestSuite


class TestSiGeometryOptimizationPipeline(SimulationParserPipelineTestSuite):
    archive_fixture = 'si_geomopt_archive'
    expected_program_name = 'FHI-aims'
