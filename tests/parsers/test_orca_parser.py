from pathlib import Path

import pytest
from nomad.datamodel import EntryArchive
from nomad.utils import get_logger
from nomad_simulations.schema_packages.basis_set import BasisSetContainer
from nomad_simulations.schema_packages.model_method import (
    CC,
    DFT,
    HF,
    MultireferenceCI,
    MultireferenceSCF,
    OrbitalLocalization,
    PerturbationMethod,
)
from nomad_simulations.schema_packages.numerical_settings import (
    LocalCorrelationSettings,
)
from nomad_simulations.schema_packages.workflow.general import SerialWorkflow

from nomad_simulation_parsers.parsers.orca.parser import OrcaParser

LOGGER = get_logger(__name__)
DATA_DIR = Path(__file__).resolve().parents[1] / 'data' / 'orca'


def _parse_orca(filename):
    parser = OrcaParser()
    archive = EntryArchive()
    parser.parse(str(DATA_DIR / filename), archive, LOGGER)
    return archive


def _methods_of_type(archive, section_cls):
    return [
        method
        for method in (archive.data.model_method or [])
        if isinstance(method, section_cls)
    ]


def _single_method(archive, section_cls):
    methods = _methods_of_type(archive, section_cls)
    assert len(methods) == 1
    return methods[0]


def _settings_of_type(method, section_cls):
    return [
        setting
        for setting in (method.numerical_settings or [])
        if isinstance(setting, section_cls)
    ]


def _basis_components_by_role(method):
    basis = _settings_of_type(method, BasisSetContainer)
    assert len(basis) == 1
    return {component.role: component for component in basis[0].basis_set_components}


def _thresholds_by_name(method):
    settings = _settings_of_type(method, LocalCorrelationSettings)
    assert len(settings) == 1
    return {threshold.name: threshold for threshold in settings[0].screening_thresholds}


def _assert_main_basis(components):
    assert components['orbital'].basis_set == 'cc-pVDZ'
    assert components['orbital'].n_total_basis_functions == 1320


def _assert_auxiliary_post_hf_basis(components):
    assert components['auxiliary_post_hf'].basis_set == 'cc-pVDZ/C'
    assert components['auxiliary_post_hf'].n_total_basis_functions == 3600


def _assert_local_correlation(method, localization):
    assert method.local_correlation is not None
    assert method.local_correlation.type == 'DLPNO'
    assert len(method.local_correlation.spaces) == 1
    assert method.local_correlation.spaces[0].kind == 'orbital'
    assert method.local_correlation.spaces[0].virtual_space_type == 'PNO'
    assert method.local_correlation.spaces[0].excitation_order == 2
    assert method.local_correlation.orbital_localization_ref is localization


def _assert_mp2_method(mp2, localization):
    assert mp2.type == 'MP'
    assert mp2.order == 2
    assert mp2.determinant == 'restricted'
    _assert_local_correlation(mp2, localization)

    components = _basis_components_by_role(mp2)
    _assert_main_basis(components)
    _assert_auxiliary_post_hf_basis(components)

    thresholds = _thresholds_by_name(mp2)
    assert thresholds['TCutPairs'].value == pytest.approx(1.0e-6)
    assert thresholds['TCutPNO'].value == pytest.approx(1.0e-8)
    assert thresholds['TCutPNOSingles'].value == pytest.approx(3.0e-10)
    assert thresholds['TCutMP2Pairs'].value == pytest.approx(1.0e-7)
    assert thresholds['TCutMKN'].value == pytest.approx(1.0e-3)
    assert thresholds['TCutPAO'].value == pytest.approx(1.0e-3)
    assert thresholds['TCutDO'].value == pytest.approx(1.0e-2)


def _assert_cc_method(cc, localization):
    assert cc.type == 'CCSD'
    assert list(cc.excitation_order) == [1, 2]
    assert cc.perturbative_correction == '(T)'
    assert list(cc.perturbative_correction_order) == [3]
    assert cc.determinant == 'restricted'
    _assert_local_correlation(cc, localization)

    components = _basis_components_by_role(cc)
    _assert_main_basis(components)
    _assert_auxiliary_post_hf_basis(components)

    thresholds = _thresholds_by_name(cc)
    assert thresholds['TCutPairs'].value == pytest.approx(1.0e-6)
    assert thresholds['TCutPairs'].applies_to == 'pair_screening'
    assert thresholds['TCutPNO'].value == pytest.approx(1.0e-8)
    assert thresholds['TCutPNO'].applies_to == 'orbital'
    assert thresholds['TCutPNOSingles'].value == pytest.approx(3.0e-10)
    assert thresholds['TCutPNOSingles'].applies_to == 'orbital'
    assert thresholds['TCutMP2Pairs'].value == pytest.approx(1.0e-7)
    assert thresholds['TCutMKN'].value == pytest.approx(1.0e-3)
    assert thresholds['TCutPAO'].value == pytest.approx(1.0e-3)
    assert thresholds['TCutEN'].value == pytest.approx(9.7e-1)
    assert thresholds['TCutDO'].value == pytest.approx(1.0e-2)
    assert thresholds['TCutTNO'].value == pytest.approx(1.0e-9)
    assert thresholds['TCutDOStrong'].value == pytest.approx(2.0e-3)
    assert thresholds['TCutMKNStrong'].value == pytest.approx(1.0e-2)
    assert thresholds['TCutMKNWeak'].value == pytest.approx(1.0e-1)
    assert thresholds['TCutDOWeak'].value == pytest.approx(4.0e-3)


def _assert_local_cc_workflow(archive, hf):
    assert isinstance(archive.workflow2, SerialWorkflow)
    assert archive.workflow2.method is not None
    assert archive.workflow2.method.initial_method is hf
    assert [task.name for task in archive.workflow2.tasks] == [
        'HF',
        'Orbital localization',
        'Local MP2',
        'Local CC',
    ]


def test_parse_file():
    parser = OrcaParser()
    archive = EntryArchive()
    parser.parse(str(DATA_DIR / 'single-point-dft.out'), archive, LOGGER)


def test_dft_xc_canonicalization():
    parser = OrcaParser()
    archive = EntryArchive()
    parser.parse(str(DATA_DIR / 'single-point-dft.out'), archive, LOGGER)

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


def test_local_cc_parsing():
    archive = _parse_orca('dlpno-coupled-cluster.out')

    system = archive.data.model_system[0]
    assert system.positions[0][0].to('angstrom').magnitude == pytest.approx(-1.487532)
    assert system.total_charge == 0
    assert system.total_spin == 0

    hf = _single_method(archive, HF)
    assert hf.type == 'RHF'
    _assert_main_basis(_basis_components_by_role(hf))

    localization = _single_method(archive, OrbitalLocalization)
    assert localization.method == 'Foster-Boys'
    assert localization.n_localized_orbitals == 54

    _assert_mp2_method(_single_method(archive, PerturbationMethod), localization)
    _assert_cc_method(_single_method(archive, CC), localization)

    assert [type(method).__name__ for method in (archive.data.model_method or [])] == [
        'HF',
        'OrbitalLocalization',
        'PerturbationMethod',
        'CC',
    ]

    _assert_local_cc_workflow(archive, hf)


@pytest.mark.parametrize(
    'case',
    [
        {
            'filename': 'CoPc_CASSCF_SS.out',
            'section_cls': MultireferenceSCF,
            'method_type': 'CASSCF',
            'active_space': (7, 5),
            'multiplicities': [4],
            'roots': [1],
        },
        {
            'filename': 'CoPc_CASSCF_SA.out',
            'section_cls': MultireferenceSCF,
            'method_type': 'CASSCF',
            'active_space': (7, 5),
            'multiplicities': [4, 2],
            'roots': [10, 40],
        },
        {
            'filename': 'CoPc_CASCI_QD.out',
            'section_cls': MultireferenceCI,
            'method_type': 'CASCI',
            'active_space': (13, 8),
            'multiplicities': [4, 2],
            'roots': [40, 115],
        },
    ],
)
def test_multireference_parsing(case):
    archive = _parse_orca(case['filename'])
    method = _single_method(archive, case['section_cls'])

    assert method.type == case['method_type']
    assert method.active_space is not None
    assert method.active_space.n_active_electrons == case['active_space'][0]
    assert method.active_space.n_active_orbitals == case['active_space'][1]
    assert method.active_space.orbital_space_type == 'CAS'
    assert list(method.state_multiplicities) == case['multiplicities']
    assert list(method.n_roots_per_multiplicity) == case['roots']
    assert len(method.state_weights) == sum(case['roots'])
