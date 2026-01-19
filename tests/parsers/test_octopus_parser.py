from nomad.datamodel import EntryArchive
from nomad.utils import get_logger

from nomad_simulation_parsers.parsers.octopus.parser import OctopusParser

LOGGER = get_logger(__name__)


def test_parse_file():
    parser = OctopusParser()
    archive = EntryArchive()
    parser.parse('tests/data/octopus/Fe_spinpol/stdout.txt', archive, LOGGER)


def test_model_method():
    parser = OctopusParser()
    archive = EntryArchive()
    parser.parse('tests/data/octopus/Fe_spinpol/stdout.txt', archive, LOGGER)
    assert archive.data.model_method is not None
    assert len(archive.data.model_method) > 0
