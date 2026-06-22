from collections.abc import Generator
from contextlib import nullcontext
from pathlib import Path

import numpy as np
import pytest
from nomad import files
from nomad.datamodel import EntryArchive, EntryMetadata
from nomad.datamodel.context import ServerContext
from nomad.utils import create_uuid, get_logger
from nomad_simulations.schema_packages.basis_set import BasisSetContainer
from nomad_simulations.schema_packages.model_method import (
    CC,
    DFT,
    HF,
    MultireferenceCI,
    MultireferencePT,
    OrbitalLocalization,
    PerturbationMethod,
    RelativityModel,
)
from nomad_simulations.schema_packages.properties.molecular_orbitals import (
    MolecularOrbitals,
)

from nomad_simulation_parsers.parsers.orca.parser import OrcaParser

LOGGER = get_logger(__name__)
DATA_DIR = Path(__file__).resolve().parents[1] / 'data' / 'orca'


class UploadForTest:
    def __init__(self, upload_id: str, upload_files: files.StagingUploadFiles) -> None:
        self.upload_id = upload_id
        self.upload_files = upload_files


def _parse_orca(filename: str) -> EntryArchive:
    parser = OrcaParser()
    archive = EntryArchive()
    parser.parse(str(DATA_DIR / filename), archive, LOGGER)
    return archive


@pytest.fixture
def archive_with_hdf5() -> Generator[EntryArchive, None, None]:
    upload_id = f'test_upload_orca_h5_{create_uuid()}'
    entry_id = 'test_entry_orca_h5'
    upload_files = files.StagingUploadFiles(upload_id, create=True)
    archive = EntryArchive(
        m_context=ServerContext(upload=UploadForTest(upload_id, upload_files)),
        metadata=EntryMetadata(upload_id=upload_id, entry_id=entry_id),
    )
    try:
        OrcaParser().parse(str(DATA_DIR / 'dft-print-MOs.out'), archive, LOGGER)
        yield archive
    finally:
        upload_files.delete()


def test_parse_file():
    _parse_orca('RI_MP2_water.out')


def test_mp2_method_and_basis_sets():
    archive = _parse_orca('RI_MP2_water.out')

    assert archive.data.program.name == 'ORCA'
    assert archive.data.program.version == '5.0.4'
    assert len(archive.data.model_method) == 1
    method = archive.data.model_method[0]
    assert isinstance(method, PerturbationMethod)
    assert method.type == 'MP'
    assert method.order == 2

    basis_sets = [
        settings
        for settings in method.numerical_settings
        if isinstance(settings, BasisSetContainer)
    ]
    assert len(basis_sets) == 1
    assert [
        component.basis_set for component in basis_sets[0].basis_set_components
    ] == ['def2-SVP', 'def2-SVP/C']


def test_multireference_methods_basis_sets_and_relativity():
    archive = _parse_orca('CoPc_CASCI_QD.out')

    casci = next(
        method
        for method in archive.data.model_method
        if isinstance(method, MultireferenceCI)
    )
    nevpt = next(
        method
        for method in archive.data.model_method
        if isinstance(method, MultireferencePT)
    )

    assert casci.type == 'CASCI'
    assert casci.active_space.n_active_electrons == 13
    assert casci.active_space.n_active_orbitals == 8
    assert list(casci.state_multiplicities) == [4, 2]
    assert list(casci.n_roots_per_multiplicity) == [40, 115]

    assert nevpt.type == 'NEVPT'
    assert nevpt.order == 2
    assert nevpt.name == 'QD-SC-NEVPT2'

    casci_basis = next(
        settings
        for settings in casci.numerical_settings
        if isinstance(settings, BasisSetContainer)
    )
    assert [
        (component.basis_set, component.role)
        for component in casci_basis.basis_set_components
    ] == [('cc-pVTZ-DK', 'orbital'), ('SARC/J', 'auxiliary_scf')]

    nevpt_basis = next(
        settings
        for settings in nevpt.numerical_settings
        if isinstance(settings, BasisSetContainer)
    )
    assert [
        (component.basis_set, component.role)
        for component in nevpt_basis.basis_set_components
    ] == [('cc-pVTZ-DK', 'orbital'), ('cc-pVTZ/C', 'auxiliary_post_hf')]

    for method in (casci, nevpt):
        relativity = next(
            contribution
            for contribution in method.contributions
            if isinstance(contribution, RelativityModel)
        )
        assert relativity.level == 'scalar'
        assert relativity.approximation == 'DKH'
        assert relativity.dkh_order == 2


def test_model_system_and_molecular_orbitals():
    archive = _parse_orca('RI_MP2_water.out')

    assert len(archive.data.model_system) == 1
    system = archive.data.model_system[0]
    assert system.is_representative
    assert [state.chemical_symbol for state in system.particle_states] == [
        'O',
        'H',
        'H',
    ]
    assert system.positions[0][2].to('angstrom').magnitude == pytest.approx(0.11779)
    assert system.total_charge == 0
    assert system.total_spin == 0

    assert len(archive.data.outputs) == 1
    output = archive.data.outputs[0]
    assert output.model_system_ref.m_proxy_value == '/data/model_system/0'
    assert len(output.molecular_orbitals) == 1

    molecular_orbitals = output.molecular_orbitals[0]
    assert isinstance(molecular_orbitals, MolecularOrbitals)
    assert molecular_orbitals.n_mo == 24
    assert molecular_orbitals.n_ao == 24
    assert molecular_orbitals.kind == 'canonical'
    assert molecular_orbitals.occupations[0] == pytest.approx(2.0)
    assert molecular_orbitals.occupations[-1] == pytest.approx(0.0)
    assert molecular_orbitals.value[0].to('electron_volt').magnitude == (
        pytest.approx(-559.0826)
    )


def test_molecular_orbital_coefficients(archive_with_hdf5):
    molecular_orbitals = archive_with_hdf5.data.outputs[0].molecular_orbitals[0]
    assert molecular_orbitals.n_mo == 77
    assert molecular_orbitals.n_ao == 77
    coefficients = molecular_orbitals.coefficients
    coefficient_context = (
        nullcontext(coefficients)
        if isinstance(coefficients, np.ndarray)
        else coefficients
    )
    with coefficient_context as coefficients:
        assert coefficients.shape == (77, 77)
        assert coefficients[0, 0] == pytest.approx(-0.000034)
        assert coefficients[5, 9] == pytest.approx(0.995146)
        assert coefficients[-1, -1] == pytest.approx(-0.006075)
    assert molecular_orbitals.occupations[28] == pytest.approx(0.0)
    assert molecular_orbitals.value[0].to('electron_volt').magnitude == (
        pytest.approx(-520.4853)
    )


def test_dft_method(archive_with_hdf5):
    dft = next(
        method
        for method in archive_with_hdf5.data.model_method
        if isinstance(method, DFT)
    )
    assert dft.reference_form == 'RKS'
    assert dft.xc.functional_key == 'B88+LYP'
    assert dft.xc.global_exact_exchange == pytest.approx(0.2)


def test_coupled_cluster_methods():
    archive = _parse_orca('dlpno-coupled-cluster.out')

    hf = next(method for method in archive.data.model_method if isinstance(method, HF))
    localization = next(
        method
        for method in archive.data.model_method
        if isinstance(method, OrbitalLocalization)
    )
    mp2 = next(
        method
        for method in archive.data.model_method
        if isinstance(method, PerturbationMethod)
    )
    cc = next(method for method in archive.data.model_method if isinstance(method, CC))

    assert hf.reference_form == 'RHF'
    assert localization.method == 'Foster-Boys'
    assert localization.n_localized_orbitals == 54
    assert mp2.type == 'MP'
    assert mp2.order == 2
    assert mp2.local_correlation.type == 'DLPNO'
    assert cc.type == 'CCSD'
    assert list(cc.excitation_order) == [1, 2]
    assert cc.perturbative_correction == '(T)'
    assert list(cc.perturbative_correction_order) == [3]
    assert cc.local_correlation.type == 'DLPNO'
