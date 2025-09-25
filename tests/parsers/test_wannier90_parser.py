from nomad.datamodel import EntryArchive
from nomad.utils import get_logger

from nomad_simulation_parsers.parsers.wannier90.parser import Wannier90Parser

LOGGER = get_logger(__name__)


def test_parse_file():
    parser = Wannier90Parser()
    archive = EntryArchive()
    parser.parse('tests/data/wannier90/lco_mlwf/lco.wout', archive, LOGGER)
