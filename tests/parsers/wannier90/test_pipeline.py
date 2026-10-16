from tests.parsers.common import SimulationParserPipelineTestSuite


class TestWannier90Pipeline(SimulationParserPipelineTestSuite):
    archive_fixture = 'lco_mlwf_archive'
    expected_program_name = 'Wannier90'
