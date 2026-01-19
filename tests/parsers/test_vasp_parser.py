from nomad.datamodel import EntryArchive
from nomad.utils import get_logger

from nomad_simulation_parsers.parsers.vasp.parser import VASPParser

LOGGER = get_logger(__name__)


def test_vasprun():
    parser = VASPParser()
    archive = EntryArchive()
    parser.parse('tests/data/vasp/AgAc_relax/vasprun.xml.relax', archive, LOGGER)


def test_outcar():
    parser = VASPParser()
    archive = EntryArchive()
    parser.parse('tests/data/vasp/AgAc_relax/OUTCAR', archive, LOGGER)


def test_model_method_vasprun():
    parser = VASPParser()
    archive = EntryArchive()
    parser.parse('tests/data/vasp/AgAc_relax/vasprun.xml.relax', archive, LOGGER)
    assert archive.data.model_method is not None
    assert len(archive.data.model_method) > 0


def test_model_method_outcar():
    parser = VASPParser()
    archive = EntryArchive()
    parser.parse('tests/data/vasp/AgAc_relax/OUTCAR', archive, LOGGER)
    assert archive.data.model_method is not None
    assert len(archive.data.model_method) > 0
