import pytest

from nomad_simulation_parsers.parsers.vasp.chgcar_parser import CHGCARFileParser
from nomad_simulation_parsers.parsers.vasp.doscar_parser import DOSCARFileParser
from nomad_simulation_parsers.parsers.vasp.outcar_parser import OutcarTextParser
from nomad_simulation_parsers.parsers.vasp.xml_parser import VasprunParser
from tests.parsers.common import assert_approx


@pytest.fixture
def outcar_parser():
    return OutcarTextParser()


@pytest.fixture
def vasprun_parser():
    return VasprunParser()


@pytest.fixture
def doscar_parser():
    return DOSCARFileParser()


@pytest.fixture
def chgcar_parser():
    return CHGCARFileParser()


@pytest.mark.unit
class TestOutcarReader:
    def test_reads_mock_header(self, tmp_path, outcar_parser):
        contents = (
            ' vasp.5.3.2 13Sep12 (build Apr 01 2013 09:32:17) complex\n\n'
            ' executed on             LinuxIFC date 2013.09.06 21:12:21\n'
        )
        mainfile = tmp_path / 'OUTCAR'
        mainfile.write_text(contents)
        parser = outcar_parser
        parser.mainfile = str(mainfile)
        source = parser.to_dict()

        assert source['header']['version'] == '5.3.2'
        assert source['header']['platform'] == 'LinuxIFC'

    def test_reads_mock_parameters_and_geometry(self, tmp_path, outcar_parser):
        contents = (
            ' vasp.5.3.2 13Sep12 (build Apr 01 2013 09:32:17) complex\n\n'
            ' executed on             LinuxIFC date 2013.09.06 21:12:21\n'
            'Startparameter for this run:\n'
            'ISPIN = 2\n'
            'ENCUT = 400\n'
            '----------------------------------------------------------------------------------------------------\n'
            'ions per type = 2 1\n'
            'NBANDS = 8\n'
            'position of ions in cartesian coordinates (Angst):\n'
            ' 1 0.100 0.200 0.300\n'
            ' 2 0.400 0.500 0.600\n\n'
        )
        mainfile = tmp_path / 'OUTCAR'
        mainfile.write_text(contents)
        outcar_parser.mainfile = str(mainfile)
        source = outcar_parser.to_dict()

        assert source['parameters'] == {'ISPIN': 2, 'ENCUT': 400}
        assert source['ions_per_type'].tolist() == [2, 1]
        assert source['nbands'] == 8
        assert source['positions'].tolist() == [
            [0.1, 0.2, 0.3],
            [0.4, 0.5, 0.6],
        ]

    @pytest.mark.parametrize(
        'kpoint_block',
        [
            (
                'Following reciprocal coordinates:\n'
                '           Coordinates               Weight\n'
                '  0.000000  0.000000  0.000000      1.000000\n'
            ),
        ],
    )
    def test_reads_kpoints_layouts(self, tmp_path, outcar_parser, kpoint_block):
        mainfile = tmp_path / 'OUTCAR'
        mainfile.write_text(kpoint_block)
        outcar_parser.mainfile = str(mainfile)

        assert_approx(outcar_parser.to_dict()['kpoints'], [0.0, 0.0, 0.0, 1.0])


@pytest.mark.unit
class TestDOSCARReader:
    def test_reads_mock_doscar(self, tmp_path, doscar_parser):
        contents = (
            'header 1\n'
            'header 2\n'
            'header 3\n'
            'header 4\n'
            'header 5\n'
            '-10.0 10.0 3 0.0 0.0\n'
            '-1.0 1.0 2.0 3.0 4.0\n'
            ' 0.0 5.0 6.0 7.0 8.0\n'
            ' 1.0 9.0 10.0 11.0 12.0\n'
        )
        mainfile = tmp_path / 'DOSCAR'
        mainfile.write_text(contents)
        doscar_parser.mainfile = str(mainfile)
        source = doscar_parser.to_dict()

        assert source['e_fermi'] == 0.0
        assert source['energies'].tolist() == [-1.0, 0.0, 1.0]
        assert source['total_dos'].tolist() == [
            [1.0, 5.0, 9.0],
            [2.0, 6.0, 10.0],
            [3.0, 7.0, 11.0],
            [4.0, 8.0, 12.0],
        ]


@pytest.mark.unit
class TestCHGCARReader:
    def test_reads_mock_chgcar(self, tmp_path, chgcar_parser):
        contents = (
            'CHGCAR mock\n'
            '1.0\n'
            '1 0 0\n'
            '0 1 0\n'
            '0 0 1\n'
            '1\n'
            'Direct\n'
            '0 0 0\n\n'
            '2 1 1\n'
            '1.0 2.0\n'
        )
        mainfile = tmp_path / 'CHGCAR'
        mainfile.write_text(contents)
        chgcar_parser.mainfile = str(mainfile)
        source = chgcar_parser.to_dict()

        assert len(source['values']) == 1
        assert_approx(source['values'][0], [[[1.0]], [[2.0]]])


@pytest.mark.unit
class TestVasprunReader:
    def test_reads_mock_xml_sections(self, tmp_path, vasprun_parser):
        contents = (
            '<?xml version="1.0" encoding="ISO-8859-1"?>\n'
            '<modeling><generator>'
            '<i name="program" type="string">vasp</i>'
            '<i name="version" type="string">6.4.2</i>'
            '</generator><atominfo><atoms>2</atoms></atominfo>'
            '<parameters><i name="NSW" type="int">3</i></parameters>'
            '<kpoints><varray name="kpointlist"><v>0 0 0</v></varray></kpoints>'
            '<structure/><calculation/><calculation/><calculation/>'
            '</modeling>\n'
        )
        mainfile = tmp_path / 'vasprun.xml'
        mainfile.write_text(contents)
        parser = vasprun_parser
        parser.filepath = str(mainfile)
        source = parser.to_dict()['modeling']

        generator = source['generator']['i']
        program = next(item for item in generator if item['@name'] == 'program')
        assert program['__value'].strip() == 'vasp'
        assert source['atominfo']['atoms'] == 2
        assert source['parameters']['i']['@name'] == 'NSW'
        assert source['kpoints']['varray']['@name'] == 'kpointlist'
        assert source['structure'] == {}
        assert len(source['calculation']) == 3
