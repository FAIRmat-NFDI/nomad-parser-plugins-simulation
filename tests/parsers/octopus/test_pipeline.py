from tests.parsers.common import SimulationParserPipelineTestSuite


class TestSiScfPipeline(SimulationParserPipelineTestSuite):
    archive_fixture = 'si_scf_archive'
    expected_program_name = 'Octopus'
