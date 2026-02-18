from nomad.datamodel import EntryArchive
from nomad.utils import get_logger
from nomad_simulations.schema_packages.model_method import DFT

from nomad_simulation_parsers.parsers.orca.parser import OrcaParser

LOGGER = get_logger(__name__)


def test_parse_file():
    parser = OrcaParser()
    archive = EntryArchive()
    parser.parse('tests/data/orca/single-point-dft.out', archive, LOGGER)


def test_dft_xc_canonicalization():
    parser = OrcaParser()
    archive = EntryArchive()
    parser.parse('tests/data/orca/single-point-dft.out', archive, LOGGER)

    dft_methods = [m for m in (archive.data.model_method or []) if isinstance(m, DFT)]
    assert len(dft_methods) == 1

    dft = dft_methods[0]
    dft.normalize(archive, LOGGER)

    assert dft.jacobs_ladder == 'meta-GGA'
    assert dft.xc is not None
    assert dft.xc.functional_key == 'TPSS'
    assert dft.xc.global_exact_exchange == 0.1
    assert [c.canonical_label for c in dft.xc.components] == [
        'XC_MGGA_C_TPSS',
        'XC_MGGA_X_TPSS',
    ]
