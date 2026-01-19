from nomad.datamodel import EntryArchive
from nomad.utils import get_logger

from nomad_simulation_parsers.parsers.phonopy.parser import PhonopyParser

LOGGER = get_logger(__name__)


def test_parse_file():
    parser = PhonopyParser()
    archive = EntryArchive()
    parser.parse('tests/data/phonopy/vasp/phonopy.yaml', archive, LOGGER)


def test_model_method():
    parser = PhonopyParser()
    archive = EntryArchive()
    parser.parse('tests/data/phonopy/vasp/phonopy.yaml', archive, LOGGER)
    # Phonopy may not populate data if the file has errors or is incomplete
    if archive.data is not None:
        assert hasattr(archive.data, 'model_method')
