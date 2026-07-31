import logging

from nomad.datamodel import EntryArchive

from nomad_parser_gsd.parsers.parser import GSDParser


def test_parse_file():
    parser = GSDParser()
    archive = EntryArchive()
    parser.parse('tests/data/trajectory.gsd', archive, logging.getLogger())

    assert archive.workflow2.name == 'test'
