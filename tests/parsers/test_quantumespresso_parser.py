from nomad.datamodel import EntryArchive
from nomad.utils import get_logger

from nomad_simulation_parsers.parsers.quantumespresso.parser import (
    QuantumEspressoParser,
)

LOGGER = get_logger(__name__)


def test_pwscf():
    parser = QuantumEspressoParser()
    archive = EntryArchive()
    parser.parse('tests/data/quantumespresso/pwscf/TiO2_opt/pw.out', archive, LOGGER)


def test_epw():
    parser = QuantumEspressoParser()
    archive = EntryArchive()
    parser.parse('tests/data/quantumespresso/epw/epw.out', archive, LOGGER)


def test_phonon():
    parser = QuantumEspressoParser()
    archive = EntryArchive()
    parser.parse('tests/data/quantumespresso/phonon/ph.out', archive, LOGGER)


def test_xspectra():
    parser = QuantumEspressoParser()
    archive = EntryArchive()
    parser.parse(
        'tests/data/quantumespresso/xspectra/ms-10734/Spectra-1-1-1/0/dipole1/xanes.out',
        archive,
        LOGGER,
    )


def test_model_method_pwscf():
    parser = QuantumEspressoParser()
    archive = EntryArchive()
    parser.parse('tests/data/quantumespresso/pwscf/TiO2_opt/pw.out', archive, LOGGER)
    assert archive.data is not None
    assert hasattr(archive.data, 'model_method')


def test_model_method_epw():
    parser = QuantumEspressoParser()
    archive = EntryArchive()
    parser.parse('tests/data/quantumespresso/epw/epw.out', archive, LOGGER)
    assert archive.data is not None
    assert hasattr(archive.data, 'model_method')


def test_model_method_phonon():
    parser = QuantumEspressoParser()
    archive = EntryArchive()
    parser.parse('tests/data/quantumespresso/phonon/ph.out', archive, LOGGER)
    assert archive.data is not None
    assert hasattr(archive.data, 'model_method')
