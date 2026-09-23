from tests.parsers.common import SimulationParserPipelineTestSuite


class TestPhonopyPipeline(SimulationParserPipelineTestSuite):
    archive_fixture = 'vasp_phonopy_archive'
    expected_program_name = 'Phonopy'
