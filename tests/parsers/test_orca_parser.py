from collections.abc import Generator
from pathlib import Path

import pytest
from nomad import files
from nomad.datamodel import EntryArchive, EntryMetadata
from nomad.datamodel.context import ServerContext
from nomad.utils import create_uuid, get_logger
from nomad_simulations.schema_packages.basis_set import BasisSetContainer
from nomad_simulations.schema_packages.model_method import DFT, PerturbationMethod
from nomad_simulations.schema_packages.properties.molecular_orbitals import (
    MolecularOrbitals,
)

from nomad_simulation_parsers.parsers.orca.parser import OrcaParser, OutParser

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
        OrcaParser().parse(str(DATA_DIR / 'orca_orbitals.out'), archive, LOGGER)
        yield archive
    finally:
        upload_files.delete()


def test_parse_file():
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
    assert len(output.electronic_eigenvalues) == 1

    molecular_orbitals = output.electronic_eigenvalues[0]
    assert isinstance(molecular_orbitals, MolecularOrbitals)
    assert molecular_orbitals.n_mo == 24
    assert molecular_orbitals.n_ao == 24
    assert molecular_orbitals.mo_type == 'canonical'
    assert molecular_orbitals.mo_occupations[0] == pytest.approx(2.0)
    assert molecular_orbitals.mo_occupations[-1] == pytest.approx(0.0)
    assert molecular_orbitals.mo_energies[0].to('electron_volt').magnitude == (
        pytest.approx(-559.0826)
    )


def test_molecular_orbital_coefficients(archive_with_hdf5):
    molecular_orbitals = archive_with_hdf5.data.outputs[0].electronic_eigenvalues[0]
    assert molecular_orbitals.n_mo == 77
    assert molecular_orbitals.n_ao == 77
    with molecular_orbitals.mo_coefficients as coefficients:
        assert coefficients.shape == (77, 77)
        assert coefficients[0, 0] == pytest.approx(-0.000034)
        assert coefficients[5, 9] == pytest.approx(0.995146)
        assert coefficients[-1, -1] == pytest.approx(-0.006075)
    assert molecular_orbitals.mo_occupations[28] == pytest.approx(0.0)
    assert molecular_orbitals.mo_energies[0].to('electron_volt').magnitude == (
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


def test_coupled_cluster_method_data():
    methods = OutParser().get_coupled_cluster_methods(
        {
            'input_file': '! DLPNO-CCSD(T)-F12',
            'single_point': {
                'cc': {
                    'coupled_cluster_type': 'CCSD',
                    'perturbative_triple_excitations_on_off': 'ON',
                    'f12_correction_on_off': 'ON',
                    'tCutPNO': 1e-7,
                }
            },
        }
    )

    assert methods == [
        {
            'type': 'CCSD',
            'excitation_order': [1, 2],
            'perturbative_correction': '(T)',
            'perturbative_correction_order': [3],
            'explicit_correlation': 'F12',
            'local_correlation': {
                'type': 'DLPNO',
                'spaces': [
                    {
                        'space_kind': 'local_virtual_space',
                        'virtual_space_type': 'PNO',
                        'excitation_order': 2,
                    }
                ],
            },
            'numerical_settings': [
                {
                    'screening_thresholds': [
                        {
                            'name': 'TCutPNO',
                            'value': 1e-7,
                            'applies_to': 'local_virtual_space',
                        }
                    ]
                }
            ],
        }
    ]


def test_multireference_and_localization_method_data():
    parser = OutParser()
    source = {
        'input_file': '%casscf maxiter 1 end',
        'single_point': {
            'casscf': {
                'n_active_electrons': 6,
                'n_active_orbitals': 5,
                'block': [
                    {
                        'multiplicity': 1,
                        'n_roots': 2,
                        'root_weights': [0.5, 0.5],
                    }
                ],
            },
            'loc': [
                {'type': 'Pipek-Mezey', 'orbital_range': [2, 8]},
            ],
        },
    }

    method = parser.get_multireference_ci_methods(source)[0]
    assert method['type'] == 'CASCI'
    assert method['active_space'] == {
        'n_active_electrons': 6,
        'n_active_orbitals': 5,
        'orbital_space_type': 'CAS',
    }
    assert method['state_treatment'] == 'state_averaged'
    assert method['state_multiplicities'] == [1]
    assert method['n_roots_per_multiplicity'] == [2]
    assert method['state_weights'] == [0.5, 0.5]
    assert parser.get_orbital_localization_methods(source) == [
        {'method': 'Pipek-Mezey', 'n_localized_orbitals': 7}
    ]
