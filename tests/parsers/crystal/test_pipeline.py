from tests.parsers.common import SimulationParserPipelineTestSuite


class TestCrystalParserPipelineSuite(SimulationParserPipelineTestSuite):
    archive_fixture = 'single_point_archive'
    expected_program_name = 'Crystal'
    composition_formula = 'Si(2)'


class TestCrystalParserPipelineSuiteGhosts(SimulationParserPipelineTestSuite):
    archive_fixture = 'ghosts_archive'
    expected_program_name = 'Crystal'
    composition_formula = 'Fe(2)O(47)Sr(16)Ti(14)'
