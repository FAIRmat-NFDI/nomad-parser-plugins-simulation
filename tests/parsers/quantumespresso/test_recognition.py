import pytest

from nomad_simulation_parsers.parsers import quantumespresso_parser


@pytest.mark.unit
class TestQuantumEspressoRecognition:
    def test_recognizes_text_output(self, tmp_path):
        contents = 'Program PWSCF v.7.3 starts\n'
        mainfile = tmp_path / 'pw.out'
        mainfile.write_text(contents)
        parser = quantumespresso_parser.load()

        assert (
            parser.is_mainfile(str(mainfile), 'text/plain', contents.encode(), contents)
            is True
        )

    @pytest.mark.parametrize('compression', ['gz', 'bz2', 'xz'])
    def test_recognizes_supported_compressions(self, tmp_path, compression):
        contents = 'Program PWSCF v.7.3 starts\n'
        mainfile = tmp_path / f'pw.out.{compression}'
        mainfile.write_text(contents)
        parser = quantumespresso_parser.load()

        assert (
            parser.is_mainfile(
                str(mainfile),
                'text/plain',
                contents.encode(),
                contents,
                compression=compression,
            )
            is True
        )

    def test_rejects_unsupported_compression(self, tmp_path):
        contents = 'Program PWSCF v.7.3 starts\n'
        mainfile = tmp_path / 'pw.out.zip'
        mainfile.write_text(contents)
        parser = quantumespresso_parser.load()

        assert not parser.is_mainfile(
            str(mainfile),
            'text/plain',
            contents.encode(),
            contents,
            compression='zip',
        )

    def test_creates_children_for_related_qe_files(self, tmp_path):
        contents = 'Program PWSCF v.7.3 starts\n'
        mainfile = tmp_path / 'pw.out'
        mainfile.write_text(contents)
        (tmp_path / 'ph.out').write_text('Program PHONON v.7.3 starts\n')
        parser = quantumespresso_parser.load()

        children = parser.is_mainfile(
            str(mainfile), 'text/plain', contents.encode(), contents
        )

        assert children == ['workflow_generic']
        assert parser.creates_children is True

    def test_creates_child_workflow_when_auxiliary_xml_is_present(self, tmp_path):
        contents = 'Program PWSCF v.7.3 starts\n'
        mainfile = tmp_path / 'pw.out'
        mainfile.write_text(contents)
        (tmp_path / 'pw.in').write_text("&CONTROL\n calculation = 'scf'\n")
        save_dir = tmp_path / 'pw.save'
        save_dir.mkdir()
        (save_dir / 'data-file-schema.xml').write_text(
            '<?xml version="1.0" encoding="UTF-8"?>\n'
            '<qes:espresso xmlns:qes="http://www.quantum-espresso.org/ns/qes/qes-1.0">\n'
            '</qes:espresso>\n'
        )
        parser = quantumespresso_parser.load()

        children = parser.is_mainfile(
            str(mainfile), 'text/plain', contents.encode(), contents
        )

        assert children == ['workflow_generic']
        assert parser.creates_children is True

    def test_creates_children_for_multiple_programs_in_one_file(self, tmp_path):
        contents = 'Program PWSCF v.7.3 starts\nProgram PHONON v.7.3 starts\n'
        mainfile = tmp_path / 'combined.out'
        mainfile.write_text(contents)
        parser = quantumespresso_parser.load()

        children = parser.is_mainfile(
            str(mainfile), 'text/plain', contents.encode(), contents
        )

        assert children == ['workflow_modules', '1 phonon']
        assert parser.creates_children is True

    @pytest.mark.parametrize(
        'contents',
        ['not a quantum espresso output\n', 'Program VASP v.6.4 starts\n'],
    )
    def test_rejects_non_quantum_espresso_output(self, tmp_path, contents):
        mainfile = tmp_path / 'output.out'
        mainfile.write_text(contents)
        parser = quantumespresso_parser.load()

        assert not parser.is_mainfile(
            str(mainfile), 'text/plain', contents.encode(), contents
        )
