from nomad.datamodel import EntryArchive
from nomad.utils import get_logger

from nomad_simulation_parsers.parsers.phonopy.parser import PhonopyParser

LOGGER = get_logger(__name__)


def test_parse_file():
    parser = PhonopyParser()
    archive = EntryArchive()
    parser.parse('tests/data/phonopy/vasp/phonopy.yaml', archive, LOGGER)
