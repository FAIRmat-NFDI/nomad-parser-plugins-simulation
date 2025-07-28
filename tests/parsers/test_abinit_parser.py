from nomad.datamodel import EntryArchive
from nomad.utils import get_logger

from nomad_simulation_parsers.parsers.abinit.parser import AbinitParser

LOGGER = get_logger(__name__)


# def test_parse_file():
#     parser = AbinitParser()
#     archive = EntryArchive()
#     parser.parse('tests/data/abinit/Fe/Fe.out', archive, LOGGER)
