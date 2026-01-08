from nomad.datamodel import EntryArchive
from nomad.utils import get_logger

from nomad_simulation_parsers.parsers.fhiaims.parser import FHIAimsParser

LOGGER = get_logger(__name__)


def test_parse_file():
    parser = FHIAimsParser()
    archive = EntryArchive()
    parser.parse('tests/data/fhiaims/Si_geomopt/out.out', archive, LOGGER)

    # Verify basic parsing succeeded
    simulation = archive.data
    assert simulation is not None, 'No simulation data in archive'
    assert simulation.model_method, 'No model_method in simulation'
