from nomad.datamodel import EntryArchive
from nomad.utils import get_logger
from nomad_simulation_parsers.parsers.crystal.parser import CrystalParser

LOGGER = get_logger(__name__)


def test_parse_file():
    parser = CrystalParser()
    archive = EntryArchive()
    parser.parse('tests/data/crystal/single_point/dft/output.out', archive, LOGGER)
