import pytest
from nomad.parsing.parsers import run_parser

from nomad_simulation_parsers.parsers import vasp_parser
from nomad_simulation_parsers.parsers.vasp.parser import VASPParser


@pytest.mark.unit
class TestVASPRecognition:
    def test_recognizes_outcar(self, tmp_path):
        contents = ' vasp.6.4.2 complex executed on Linux date 2025.01.01 12:00:00\n'
        mainfile = tmp_path / 'OUTCAR'
        mainfile.write_text(contents)
        parser = vasp_parser.load()

        assert (
            parser.is_mainfile(str(mainfile), 'text/plain', contents.encode(), contents)
            is True
        )

    def test_recognizes_vasprun_xml(self, tmp_path):
        contents = (
            '<?xml version="1.0" encoding="ISO-8859-1"?>\n'
            '<modeling><generator><i name="program" type="string">vasp</i>'
            '</generator></modeling>\n'
        )
        mainfile = tmp_path / 'vasprun.xml'
        mainfile.write_text(contents)
        parser = vasp_parser.load()

        assert (
            parser.is_mainfile(str(mainfile), 'text/xml', contents.encode(), contents)
            is True
        )

    def test_recognizes_phonopy_workflow(self, tmp_path):
        (tmp_path / 'phonopy.yaml').write_text('phonopy:\n  calculator: vasp\n')
        displacement_dir = tmp_path / 'disp-001'
        displacement_dir.mkdir()
        contents = (
            '<?xml version="1.0" encoding="ISO-8859-1"?>\n'
            '<modeling><generator><i name="program" type="string">vasp</i>'
            '</generator></modeling>\n'
        )
        mainfile = displacement_dir / 'vasprun.xml'
        mainfile.write_text(contents)
        parser = vasp_parser.load()

        assert parser.is_mainfile(
            str(mainfile), 'text/xml', contents.encode(), contents
        ) == ['phonon_workflow:../phonopy.yaml']
        assert parser.creates_children

    def test_regular_vasp_matching_keeps_child_creation_enabled(self, tmp_path):
        contents = ' vasp.6.4.2 complex executed on Linux date 2025.01.01 12:00:00\n'
        mainfile = tmp_path / 'OUTCAR'
        mainfile.write_text(contents)
        parser = vasp_parser.load()

        assert (
            parser.is_mainfile(str(mainfile), 'text/plain', contents.encode(), contents)
            is True
        )
        assert parser.creates_children

    def test_phonopy_workflow_key_creates_child_archive(self, tmp_path, monkeypatch):
        (tmp_path / 'phonopy.yaml').write_text('phonopy:\n  calculator: vasp\n')
        displacement_dir = tmp_path / 'disp-001'
        displacement_dir.mkdir()
        contents = (
            '<?xml version="1.0" encoding="ISO-8859-1"?>\n'
            '<modeling><generator><i name="program" type="string">vasp</i>'
            '</generator></modeling>\n'
        )
        mainfile = displacement_dir / 'vasprun.xml'
        mainfile.write_text(contents)
        parser = vasp_parser.load()
        mainfile_keys = parser.is_mainfile(
            str(mainfile), 'text/xml', contents.encode(), contents
        )
        monkeypatch.setattr(VASPParser, 'parse', lambda *args, **kwargs: None)

        archives = run_parser(str(mainfile), parser, mainfile_keys)

        assert len(archives) == 2
        assert archives[1].metadata.mainfile_key == 'phonon_workflow:../phonopy.yaml'

    def test_does_not_recognize_outcar_with_vasprun_sibling(self, tmp_path):
        contents = ' vasp.6.4.2 complex executed on Linux date 2025.01.01 12:00:00\n'
        mainfile = tmp_path / 'OUTCAR'
        mainfile.write_text(contents)
        (tmp_path / 'vasprun.xml').write_text('')
        parser = vasp_parser.load()

        assert vasp_parser.mainfile_alternative is True
        assert not parser.is_mainfile(
            str(mainfile), 'text/plain', contents.encode(), contents
        )

    @pytest.mark.parametrize('compression', ['gz', 'bz2', 'xz'])
    def test_recognizes_supported_compression(self, tmp_path, compression):
        contents = (
            '<?xml version="1.0" encoding="ISO-8859-1"?>\n'
            '<modeling><generator><i name="program" type="string">vasp</i>'
            '</generator></modeling>\n'
        )
        mainfile = tmp_path / f'vasprun.xml.{compression}'
        mainfile.write_text(contents)
        parser = vasp_parser.load()

        assert (
            parser.is_mainfile(
                str(mainfile),
                'text/xml',
                contents.encode(),
                contents,
                compression=compression,
            )
            is True
        )

    def test_rejects_unsupported_compression(self, tmp_path):
        contents = (
            '<?xml version="1.0" encoding="ISO-8859-1"?>\n'
            '<modeling><generator><i name="program" type="string">vasp</i>'
            '</generator></modeling>\n'
        )
        mainfile = tmp_path / 'vasprun.xml.zip'
        mainfile.write_text(contents)
        parser = vasp_parser.load()

        assert not parser.is_mainfile(
            str(mainfile), 'text/xml', contents.encode(), contents, compression='zip'
        )

    @pytest.mark.parametrize(
        'contents',
        [
            'not a VASP output\n',
            ' vasp.6.4.2 started\n',
            '<?xml version="1.0"?><modeling><generator><i name="program">qe</i>',
        ],
    )
    def test_rejects_non_vasp_output(self, tmp_path, contents):
        mainfile = tmp_path / 'output'
        mainfile.write_text(contents)
        parser = vasp_parser.load()

        assert not parser.is_mainfile(
            str(mainfile), 'text/plain', contents.encode(), contents
        )
