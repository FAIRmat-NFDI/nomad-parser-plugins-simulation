import pytest

from nomad_simulation_parsers.parsers import orca_parser

ORCA_BANNER = (
    '                                 *****************\n'
    '                                 * O   R   C   A *\n'
    '                                 *****************\n'
)


@pytest.mark.unit
class TestOrcaRecognition:
    def test_recognizes_orca_output(self, tmp_path):
        contents = ORCA_BANNER + 'Program Version 5.0.4\n'
        mainfile = tmp_path / 'mol.out'
        mainfile.write_text(contents)
        parser = orca_parser.load()

        assert (
            parser.is_mainfile(str(mainfile), 'text/plain', contents.encode(), contents)
            is True
        )

    @pytest.mark.parametrize(
        'contents',
        [
            'not an orca output\n',
            '* G A U S S I A N *\n',  # a different code banner
            '* ORCA *\n',  # the letters must be spaced as in the real banner
        ],
    )
    def test_rejects_non_orca_output(self, tmp_path, contents):
        mainfile = tmp_path / 'mol.out'
        mainfile.write_text(contents)
        parser = orca_parser.load()

        assert not parser.is_mainfile(
            str(mainfile), 'text/plain', contents.encode(), contents
        )
