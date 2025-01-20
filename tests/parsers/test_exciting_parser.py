from nomad.datamodel import EntryArchive

from nomad_simulation_parsers.parsers.exciting.parser import ExcitingParser


def test_parse_file():
    parser = ExcitingParser()
    archive = EntryArchive()
    parser.parse('test/data/exciting/INFO.OUT', archive)
