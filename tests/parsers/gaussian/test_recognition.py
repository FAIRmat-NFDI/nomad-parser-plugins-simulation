import pytest

from nomad_simulation_parsers.parsers import gaussian_parser

GAUSSIAN_CITATION = ' Cite this work as:\n Gaussian 09, Revision D.01,\n'


@pytest.mark.unit
class TestGaussianRecognition:
    def test_recognizes_gaussian_output(self, tmp_path):
        mainfile = tmp_path / 'calculation.out'
        mainfile.write_text(GAUSSIAN_CITATION)
        parser = gaussian_parser.load()

        assert (
            parser.is_mainfile(
                str(mainfile),
                'text/plain',
                GAUSSIAN_CITATION.encode(),
                GAUSSIAN_CITATION,
            )
            is True
        )

    @pytest.mark.parametrize(
        ('filename', 'contents'),
        [
            ('calculation.out', 'not a Gaussian output\n'),
            ('calculation.out', 'Gaussian 09, Revision D.01\n'),
            ('calculation.txt', GAUSSIAN_CITATION),
        ],
    )
    def test_rejects_non_gaussian_mainfile(self, tmp_path, filename, contents):
        mainfile = tmp_path / filename
        mainfile.write_text(contents)
        parser = gaussian_parser.load()

        assert not parser.is_mainfile(
            str(mainfile), 'text/plain', contents.encode(), contents
        )
