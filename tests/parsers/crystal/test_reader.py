from pathlib import Path

import pytest

from nomad_simulation_parsers.parsers.crystal.file_parser import (
    F25Parser,
    OutputParser,
    to_float,
)
from nomad_simulation_parsers.parsers.crystal.parser import (
    CrystalF25Parser,
    CrystalOutputParser,
)

MINIMAL_OUTPUT = """
date ma 10.10.2016 13.59.45 +0300
hostname test
system Linux
user user
output data in output.out
input data in input.d12
CRYSTAL
DFT
EXCHANGE
PBE
END
\n * CRYSTAL14 *
 * version : 1
 EEEEEEEEEE STARTING  DATE 10 10 2016 TIME 13:59:45.000
 TEST TITLE

"""


def read_output_source(path):
    parser = OutputParser()
    parser.mainfile = str(path)
    return parser.to_dict()


@pytest.mark.unit
def test_extracts_crystal_numeric_notation():
    assert to_float('2**-3') == 0.125
    assert to_float('2**4') == 16


@pytest.mark.unit
def test_extracts_crystal_start_timestamp():
    timestamp = CrystalOutputParser().to_unix_time('10 10 2016 TIME 13:59:45.000')

    assert timestamp == pytest.approx(1476100785.0)


@pytest.mark.unit
def test_extracts_header_and_functional_from_minimal_output(tmp_path):
    mainfile = tmp_path / 'minimal.out'
    mainfile.write_text(MINIMAL_OUTPUT)

    source = read_output_source(mainfile)

    assert source['program_version'] == '14'
    assert source['distribution'] == 'version : 1'
    assert source['start_timestamp'] == '10 10 2016 TIME 13:59:45.000'
    assert source['title'] == 'TEST TITLE'
    assert source['dft']['exchange'] == 'PBE'


@pytest.mark.unit
def test_extracts_fort25_band_and_dos_sections(tmp_path):
    parser = F25Parser()
    parser.mainfile = str(
        Path(__file__).resolve().parents[2]
        / 'data'
        / 'crystal'
        / 'band_structure'
        / 'nacl_hf'
        / 'NaCl.f25'
    )
    source = parser.to_dict()

    assert len(source['segments']) == 4
    assert source['dos'] is not None
    crystal_parser = CrystalF25Parser()
    assert crystal_parser.get_band_structures(source['segments'])[0]['value'].shape == (
        20,
        18,
    )
    assert crystal_parser.get_dos(source['dos'])[0]['values'].shape == (302,)
