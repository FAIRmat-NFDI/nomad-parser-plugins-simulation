import pytest

from nomad_simulation_parsers.parsers import yambo_parser


@pytest.mark.unit
class TestYamboRecognition:
    def test_recognizes_yambo_report(self, tmp_path):
        contents = (
            'Version 5.0.4 Revision 1\nMPI+HDF5_IO Build\nhttp://www.yambo-code.org\n'
        )
        mainfile = tmp_path / 'r-example'
        mainfile.write_text(contents)

        assert yambo_parser.load().is_mainfile(
            str(mainfile), 'text/plain', contents.encode(), contents
        )

    @pytest.mark.parametrize(
        'contents', ['', 'not a YAMBO report\n', 'Version 5.0.4\n']
    )
    def test_rejects_non_yambo_output(self, tmp_path, contents):
        mainfile = tmp_path / 'output'
        mainfile.write_text(contents)

        assert not yambo_parser.load().is_mainfile(
            str(mainfile), 'text/plain', contents.encode(), contents
        )
