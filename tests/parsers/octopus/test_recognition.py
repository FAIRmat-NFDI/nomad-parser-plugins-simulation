import pytest

from nomad_simulation_parsers.parsers import octopus_parser


@pytest.mark.unit
class TestOctopusRecognition:
    def test_recognizes_octopus_output(self, tmp_path):
        contents = 'header\n|0) ~ (0) |\nRunning octopus\n'
        mainfile = tmp_path / 'stdout.txt'
        mainfile.write_text(contents)
        parser = octopus_parser.load()

        assert (
            parser.is_mainfile(str(mainfile), 'text/plain', contents.encode(), contents)
            is True
        )

    @pytest.mark.parametrize(
        'contents', ['not an octopus output\n', 'Running octave\n']
    )
    def test_rejects_non_octopus_output(self, tmp_path, contents):
        mainfile = tmp_path / 'stdout.txt'
        mainfile.write_text(contents)
        parser = octopus_parser.load()

        assert not parser.is_mainfile(
            str(mainfile), 'text/plain', contents.encode(), contents
        )
