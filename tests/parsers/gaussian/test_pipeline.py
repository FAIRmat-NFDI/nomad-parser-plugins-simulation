from tests.parsers.common import SimulationParserPipelineTestSuite


class TestAlSCFPipeline(SimulationParserPipelineTestSuite):
    archive_fixture = 'al_scf_archive'
    expected_program_name = 'Gaussian'
