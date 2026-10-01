import pytest

from nomad_simulation_parsers.parsers import lobster_parser


@pytest.mark.unit
class TestLobsterRecognition:
    def test_recognizes_lobster_output(self, tmp_path):
        contents = (
            'LOBSTER v5.1.1\nstarting on host node on 2025-01-02 at 03:04:05 UTC\n'
        )
        mainfile = tmp_path / 'lobsterout'
        mainfile.write_text(contents)
        parser = lobster_parser.load()

        assert parser.is_mainfile(
            str(mainfile), 'text/plain', contents.encode(), contents
        )

    @pytest.mark.parametrize(
        'contents',
        [
            'not a LOBSTER output\n',
            'LOBSTER output without a version\n',
            'LOBSTER vX\nstarted on host node\n',
        ],
    )
    def test_rejects_non_lobster_output(self, tmp_path, contents):
        mainfile = tmp_path / 'lobsterout'
        mainfile.write_text(contents)
        parser = lobster_parser.load()

        assert not parser.is_mainfile(
            str(mainfile), 'text/plain', contents.encode(), contents
        )
