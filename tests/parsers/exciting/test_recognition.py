import pytest

from nomad_simulation_parsers.parsers import exciting_parser


@pytest.mark.unit
class TestExcitingRecognition:
    def test_recognizes_exciting_output(self, tmp_path):
        contents = 'EXCITING 2025 started\nheader\nAll units are atomic \n'
        mainfile = tmp_path / 'INFO.OUT'
        mainfile.write_text(contents)
        parser = exciting_parser.load()

        assert (
            parser.is_mainfile(str(mainfile), 'text/plain', contents.encode(), contents)
            is True
        )

    @pytest.mark.parametrize(
        'contents',
        [
            'not an exciting output\n',
            'EXCITING 2025 started\nunits are SI\n',
        ],
    )
    def test_rejects_non_exciting_output(self, tmp_path, contents):
        mainfile = tmp_path / 'INFO.OUT'
        mainfile.write_text(contents)
        parser = exciting_parser.load()

        assert not parser.is_mainfile(
            str(mainfile), 'text/plain', contents.encode(), contents
        )
