import pytest

from tests.parsers.common import SimulationParserPipelineTestSuite


@pytest.mark.large_fixture
class TestPhonopyPipeline(SimulationParserPipelineTestSuite):
    archive_fixture = 'vasp_phonopy_archive'
    expected_program_name = 'Phonopy'


class TestCp2kPhonopyPipeline(SimulationParserPipelineTestSuite):
    archive_fixture = 'cp2k_hexagonal_noncanonical_phonopy_archive'
    expected_program_name = 'Phonopy'
