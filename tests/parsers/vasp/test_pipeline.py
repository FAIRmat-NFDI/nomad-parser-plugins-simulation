from tests.parsers.common import SimulationParserPipelineTestSuite


class TestVASPOutcarPipeline(SimulationParserPipelineTestSuite):
    archive_fixture = 'agac_relax_outcar_archive'
    expected_program_name = 'VASP'
