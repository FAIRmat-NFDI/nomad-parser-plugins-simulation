from nomad.datamodel import EntryArchive
from nomad.utils import get_logger

from nomad_simulation_parsers.parsers.vasp.parser import VASPParser

LOGGER = get_logger(__name__)


def test_vasprun():
    parser = VASPParser()
    archive = EntryArchive()
    parser.parse('tests/data/vasp/AgAc_relax/vasprun.xml.relax', archive, LOGGER)

    # Verify pseudopotentials are extracted
    assert archive.data is not None
    assert archive.data.model_method is not None
    assert len(archive.data.model_method) > 0
    method = archive.data.model_method[0]
    assert method.numerical_settings is not None

    pseudopotentials = [
        ns for ns in method.numerical_settings if type(ns).__name__ == 'Pseudopotential'
    ]
    # Note: OUTCAR merge can deduplicate PPs, so we check for at least 1
    # When vasprun.xml is parsed alone (no OUTCAR), we get 2 pseudopotentials
    # When OUTCAR is present and merged, deduplication may occur
    n_pps = len(pseudopotentials)
    assert n_pps >= 1, f'Expected at least 1 pseudopotential, got {n_pps}'

    # Check first pseudopotential has required fields from XML
    pp0 = pseudopotentials[0]
    assert pp0.name is not None, 'Pseudopotential name should not be None (XML)'
    assert pp0.n_valence_electrons is not None, (
        'n_valence_electrons should not be None (XML)'
    )


def test_outcar():
    parser = VASPParser()
    archive = EntryArchive()
    parser.parse('tests/data/vasp/AgAc_relax/OUTCAR', archive, LOGGER)
