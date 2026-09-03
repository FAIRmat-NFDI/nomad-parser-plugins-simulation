import pytest

from nomad_simulation_parsers.parsers import abinit_parser


@pytest.mark.unit
class TestAbinitRecognition:
    def test_recognizes_abinit_header(self, tmp_path):
        contents = '\n.Version 9.10.4 of ABINIT\n'
        mainfile = tmp_path / 'calculation.out'
        mainfile.write_text(contents)
        parser = abinit_parser.load()

        assert (
            parser.is_mainfile(str(mainfile), 'text/plain', contents.encode(), contents)
            is True
        )


    @pytest.mark.parametrize(
        ('contents', 'compression'),
        [
            ('Version 9.10.4 of another program\n', None),
            ('\n.Version 9.10.4 of ABINIT\n', 'zip'),
        ],
    )
    def test_rejects_non_abinit_or_unsupported_compression(
        self, tmp_path, contents, compression
    ):
        mainfile = tmp_path / 'calculation.out'
        mainfile.write_text(contents)
        parser = abinit_parser.load()

        assert not parser.is_mainfile(
            str(mainfile),
            'text/plain',
            contents.encode(),
            contents,
            compression=compression,
        )
