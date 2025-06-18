from nomad.datamodel import EntryArchive
from nomad.utils import get_logger
from nomad_simulation_parsers.parsers.fhiaims.parser import FHIAimsParser

LOGGER = get_logger(__name__)


def test_parse_file():
    parser = FHIAimsParser()
    archive = EntryArchive()
    parser.parse('tests/data/fhiaims/Si_geomopt/out.out', archive, LOGGER)
