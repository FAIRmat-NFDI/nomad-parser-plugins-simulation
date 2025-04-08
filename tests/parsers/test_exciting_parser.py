from nomad.datamodel import EntryArchive

from nomad_simulation_parsers.parsers.exciting.parser import ExcitingParser


def test_parse_file():
    parser = ExcitingParser()
    archive = EntryArchive()
    parser.parse('tests/data/exciting/GaO_strucopt/INFO.OUT', archive)
