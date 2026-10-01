from tests.parsers.common import SimulationParserPipelineTestSuite


class TestPWSCFPipeline(SimulationParserPipelineTestSuite):
    archive_fixture = 'pwscf_archive'
    expected_program_name = 'Quantum Espresso'
