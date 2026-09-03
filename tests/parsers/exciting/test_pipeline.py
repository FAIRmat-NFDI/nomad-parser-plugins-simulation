import pytest

from tests.parsers.common import SimulationParserPipelineTestSuite


# TODO: Remove this skip once EntryMetadata.auxiliary_files is available.
@pytest.mark.skip(reason='requires EntryMetadata.auxiliary_files in nomad.datamodel')
def test_exciting_parser_records_parsed_text_blocks(c_minimal_archive):
    files = c_minimal_archive.metadata.auxiliary_files
    blocks = [block for file in files for block in file.parsed_blocks]

    assert blocks
    assert {file.file_name.rsplit('/', maxsplit=1)[-1] for file in files} == {
        'INFO.OUT',
        'EIGVAL.OUT',
    }
    assert all(block.start < block.end for block in blocks)


class TestCMinimalPipeline(SimulationParserPipelineTestSuite):
    archive_fixture = 'c_minimal_archive'
    expected_program_name = 'exciting'
