import pytest

from nomad_simulation_parsers.parsers import yambo_parser


@pytest.mark.unit
class TestYamboRecognition:
    @pytest.mark.parametrize(
        'contents', [
            'Version 5.0.4 Revision 1\nMPI+HDF5_IO Build\nhttp://www.yambo-code.org\n',
            ('Version 5.3.0 Revision 23927 Hash (prev commit) 1730222ea\n'
            'Branch is \n' 'Serial+HDF5_IO Build\n' 'http://www.yambo-code.eu\n'),
        ]
    )
    def test_recognizes_yambo_report(self, tmp_path, contents):
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
