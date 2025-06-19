from nomad.datamodel import EntryArchive
from nomad.utils import get_logger

from nomad_simulation_parsers.parsers.ams.parser import AMSParser

LOGGER = get_logger(__name__)


def test_parse_file():
    parser = AMSParser()
    archive = EntryArchive()
    parser.parse('tests/data/ams/scf/phenylrSmall-metagga.out', archive, LOGGER)
