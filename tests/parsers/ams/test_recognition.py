import pytest

from nomad_simulation_parsers.parsers import ams_parser


@pytest.mark.unit
class TestAMSRecognition:
    def test_recognizes_ams_header(self, tmp_path):
        contents = (
            '\n *                              |     A M S     | '
            '                             *\n'
        )
        mainfile = tmp_path / 'calculation.out'
        mainfile.write_text(contents)
        parser = ams_parser.load()

        assert (
            parser.is_mainfile(str(mainfile), 'text/plain', contents.encode(), contents)
            is True
        )

    def test_rejects_non_ams_or_unsupported_compression(self, tmp_path):
        contents = 'This is not an AMS output.\n'
        mainfile = tmp_path / 'calculation.out'
        mainfile.write_text(contents)
        parser = ams_parser.load()

        assert not parser.is_mainfile(
            str(mainfile), 'text/plain', contents.encode(), contents
        )
        assert not parser.is_mainfile(
            str(mainfile), 'text/plain', contents.encode(), contents, compression='zip'
        )
