import pytest

from nomad_simulation_parsers.parsers import fhiaims_parser


@pytest.mark.unit
class TestFHIAimsRecognition:
    def test_recognizes_fhiaims_output(self, tmp_path):
        contents = 'header\n  Invoking FHI-aims ...\n'
        mainfile = tmp_path / 'aims.out'
        mainfile.write_text(contents)
        parser = fhiaims_parser.load()

        assert parser.is_mainfile(
            str(mainfile), 'text/plain', contents.encode(), contents
        ) is True

    @pytest.mark.parametrize('contents', ['not an FHI-aims output\n', 'aims started\n'])
    def test_rejects_non_fhiaims_output(self, tmp_path, contents):
        mainfile = tmp_path / 'aims.out'
        mainfile.write_text(contents)
        parser = fhiaims_parser.load()

        assert not parser.is_mainfile(
            str(mainfile), 'text/plain', contents.encode(), contents
        )
