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

from devtools import debug
def test_gipaw_nmr():
    parser = QuantumEspressoParser()
    archive = EntryArchive()
    parser.parse(
        'tests/data/quantumespresso/gipaw/quartz-nmr/quartz-nmr.out',
        archive,
        LOGGER,
    )
    debug(archive.data.outputs[0])
    debug(archive.data.outputs[0].magnetic_shieldings)



def test_gipaw_scf():
    parser = QuantumEspressoParser()
    archive = EntryArchive()
    parser.parse(
        'tests/data/quantumespresso/gipaw/quartz-nmr/quartz-scf.out',
        archive,
        LOGGER,
    )
    debug(archive.data.outputs[0])
