import pytest

from nomad_simulation_parsers.parsers import wannier90_parser
from nomad_simulation_parsers.parsers.wannier90.parser import Wannier90Parser


@pytest.mark.unit
class TestWannier90Recognition:
    def test_recognizes_wout(self, tmp_path):
        contents = '|                   WANNIER90                       |\n'
        mainfile = tmp_path / 'calculation.wout'
        mainfile.write_text(contents)
        parser = wannier90_parser.load()
        assert parser.is_mainfile(
            str(mainfile), 'text/plain', contents.encode(), contents
        )

    def test_returns_dft_child_when_present(self, tmp_path):
        contents = '|                   WANNIER90                       |\n'
        mainfile = tmp_path / 'lco.wout'
        mainfile.write_text(contents)
        dft = tmp_path / 'lco.OUTCAR'
        dft.write_text('')
        result = Wannier90Parser().is_mainfile(
            str(mainfile), 'text/plain', contents.encode(), contents
        )
        assert result == [str(dft)]

    @pytest.mark.parametrize(
        'contents', ['not a Wannier90 output\n', '| WANNIER90 without a calculation\n']
    )
    def test_rejects_non_wannier90_output(self, tmp_path, contents):
        mainfile = tmp_path / 'calculation.wout'
        mainfile.write_text(contents)
        parser = wannier90_parser.load()
        assert not parser.is_mainfile(
            str(mainfile), 'text/plain', contents.encode(), contents
        )
