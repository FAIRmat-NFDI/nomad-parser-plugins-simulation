import pytest

from nomad_simulation_parsers.parsers import gpaw_parser


@pytest.mark.unit
class TestGPAWRecognition:
    @pytest.mark.parametrize('filename', ['calculation.gpw', 'calculation.gpw2'])
    def test_recognizes_gpaw_extensions(self, tmp_path, filename):
        mainfile = tmp_path / filename
        mainfile.write_bytes(b'not parsed by recognition')
        parser = gpaw_parser.load()

        assert (
            parser.is_mainfile(
                str(mainfile),
                'application/octet-stream',
                mainfile.read_bytes(),
                None,
            )
            is True
        )

    @pytest.mark.parametrize(
        ('filename', 'mime'),
        [
            ('calculation.txt', 'application/octet-stream'),
            ('calculation.gpw', 'text/plain'),
        ],
    )
    def test_rejects_non_gpaw_mainfiles(self, tmp_path, filename, mime):
        mainfile = tmp_path / filename
        contents = b'not a GPAW file'
        mainfile.write_bytes(contents)
        parser = gpaw_parser.load()

        assert not parser.is_mainfile(str(mainfile), mime, contents, contents.decode())
