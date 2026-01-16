from nomad.datamodel import EntryArchive
from nomad.utils import get_logger

from nomad_simulation_parsers.parsers.exciting.parser import ExcitingParser

LOGGER = get_logger(__name__)


def test_parse_file():
    parser = ExcitingParser()
    archive = EntryArchive()
    parser.parse('tests/data/exciting/GaO_strucopt/INFO.OUT', archive, LOGGER)

    simulation = archive.data
    assert simulation is not None, 'No simulation data in archive'
    assert simulation.model_method, 'No model_method in simulation'
