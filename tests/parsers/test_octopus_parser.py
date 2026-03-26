from nomad.datamodel import EntryArchive
from nomad.utils import get_logger

from nomad_simulation_parsers.parsers.octopus.parser import OctopusParser

LOGGER = get_logger(__name__)


def test_parse_file():
    parser = OctopusParser()
    archive = EntryArchive()
    parser.parse('tests/data/octopus/Fe_spinpol/stdout.txt', archive, LOGGER)

    outputs = archive.data.outputs
    assert outputs is not None
    assert len(outputs) > 0
    assert outputs[-1].total_energies is not None
    assert outputs[-1].total_forces is not None
    assert outputs[-1].electronic_eigenvalues is not None
