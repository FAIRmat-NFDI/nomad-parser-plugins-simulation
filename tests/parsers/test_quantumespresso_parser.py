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



def test_pwscf_xml():
    parser = QuantumEspressoParser()
    archive = EntryArchive()
    parser.parse('tests/data/quantumespresso/pwscf/TiO2_opt/TiO2.save/data-file-schema.xml', archive, LOGGER)
    # parser.parse('tests/data/quantumespresso/gipaw/scf_xml_nmr_xml/quartz.xml', archive, LOGGER)
    from devtools import debug
    debug(archive.data)
    debug(archive.data.model_system)




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
def test_gipaw_nmr_text():
    parser = QuantumEspressoParser()
    archive = EntryArchive()
    parser.parse(
        'tests/data/quantumespresso/gipaw/scf_out_nmr_out_741/quartz-nmr.out',
        archive,
        LOGGER,
    )
    debug(archive)
    debug(archive.data)
    debug(archive.data.outputs)
    debug(archive.data.outputs[0])
    debug(archive.data.outputs[0].magnetic_shieldings)
    debug(archive.data.outputs[0].magnetic_susceptibilities)


def test_gipaw_nmr_text():
    parser = QuantumEspressoParser()
    archive = EntryArchive()
    parser.parse(
        'tests/data/quantumespresso/gipaw/scf_out_nmr_out_741/quartz-nmr.out',
        archive,
        LOGGER,
    )
    debug(archive)
    debug(archive.data)
    debug(archive.data.outputs)
    debug(archive.data.outputs[0])
    debug(archive.data.outputs[0].magnetic_shieldings)
    debug(archive.data.outputs[0].magnetic_susceptibilities)
    

def test_gipaw_nmr_xml():
    parser = QuantumEspressoParser()
    archive = EntryArchive()
    parser.parse(
        'tests/data/quantumespresso/gipaw/scf_xml_nmr_xml/quartz-nmr.xml',
        archive,
        LOGGER,
    )
    debug(archive.data.outputs)
    debug(archive.data.outputs[0].magnetic_shieldings)




