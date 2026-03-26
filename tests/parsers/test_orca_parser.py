from pathlib import Path

import pytest
from nomad.datamodel import EntryArchive
from nomad.utils import get_logger
from nomad_simulations.schema_packages.basis_set import BasisSetContainer
from nomad_simulations.schema_packages.model_method import (
    CC,
    DFT,
    HF,
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
    parser = OrcaParser()
    archive = EntryArchive()
    parser.parse(str(DATA_DIR / 'dlpno-coupled-cluster.out'), archive, LOGGER)

    hf_methods = [m for m in (archive.data.model_method or []) if isinstance(m, HF)]
    localization_methods = [
        m for m in (archive.data.model_method or []) if isinstance(m, OrbitalLocalization)
    ]
    mp2_methods = [
        m
        for m in (archive.data.model_method or [])
        if isinstance(m, PerturbationMethod)
    ]
    cc_methods = [m for m in (archive.data.model_method or []) if isinstance(m, CC)]
    assert len(hf_methods) == 1
    assert len(localization_methods) == 1
    assert len(mp2_methods) == 1
    assert len(cc_methods) == 1

    hf = hf_methods[0]
    assert hf.type == 'RHF'
    hf_basis = [
        setting for setting in (hf.numerical_settings or []) if isinstance(setting, BasisSetContainer)
    ]
    assert len(hf_basis) == 1
    hf_components = {component.role: component for component in hf_basis[0].basis_set_components}
    assert hf_components['orbital'].basis_set == 'cc-pVDZ'
    assert hf_components['orbital'].n_total_basis_functions == 1320

    localization = localization_methods[0]
    assert localization.method == 'Foster-Boys'
    assert localization.n_localized_orbitals == 54

    mp2 = mp2_methods[0]
    assert mp2.type == 'MP'
    assert mp2.order == 2
    assert mp2.determinant == 'restricted'
    assert mp2.local_correlation is not None
    assert mp2.local_correlation.type == 'DLPNO'
    assert len(mp2.local_correlation.spaces) == 1
    assert mp2.local_correlation.spaces[0].kind == 'orbital'
    assert mp2.local_correlation.spaces[0].virtual_space_type == 'PNO'
    assert mp2.local_correlation.spaces[0].excitation_order == 2
    assert mp2.local_correlation.orbital_localization_ref is localization
    mp2_local_settings = [
        setting
        for setting in (mp2.numerical_settings or [])
        if isinstance(setting, LocalCorrelationSettings)
    ]
    assert len(mp2_local_settings) == 1
    mp2_basis = [
        setting
        for setting in (mp2.numerical_settings or [])
        if isinstance(setting, BasisSetContainer)
    ]
    assert len(mp2_basis) == 1
    mp2_components = {
        component.role: component for component in mp2_basis[0].basis_set_components
    }
    assert mp2_components['orbital'].basis_set == 'cc-pVDZ'
    assert mp2_components['orbital'].n_total_basis_functions == 1320
    assert mp2_components['auxiliary_post_hf'].basis_set == 'cc-pVDZ/C'
    assert mp2_components['auxiliary_post_hf'].n_total_basis_functions == 3600

    mp2_thresholds = {
        threshold.name: threshold
        for threshold in mp2_local_settings[0].screening_thresholds
    }
    assert mp2_thresholds['TCutPairs'].value == pytest.approx(1.0e-6)
    assert mp2_thresholds['TCutPNO'].value == pytest.approx(1.0e-8)
    assert mp2_thresholds['TCutPNOSingles'].value == pytest.approx(3.0e-10)
    assert mp2_thresholds['TCutMP2Pairs'].value == pytest.approx(1.0e-7)

    cc = cc_methods[0]
    assert cc.type == 'CCSD'
    assert list(cc.excitation_order) == [1, 2]
    assert cc.perturbative_correction == '(T)'
    assert list(cc.perturbative_correction_order) == [3]
    assert cc.determinant == 'restricted'

    assert cc.local_correlation is not None
    assert cc.local_correlation.type == 'DLPNO'
    assert len(cc.local_correlation.spaces) == 1
    assert cc.local_correlation.spaces[0].kind == 'orbital'
    assert cc.local_correlation.spaces[0].virtual_space_type == 'PNO'
    assert cc.local_correlation.spaces[0].excitation_order == 2
    assert cc.local_correlation.orbital_localization_ref is localization

    assert [type(method).__name__ for method in (archive.data.model_method or [])] == [
        'HF',
        'OrbitalLocalization',
        'PerturbationMethod',
        'CC',
    ]

    assert isinstance(archive.workflow2, SerialWorkflow)
    assert archive.workflow2.method is not None
    assert archive.workflow2.method.initial_method is hf
    assert [task.name for task in archive.workflow2.tasks] == [
        'HF',
        'Orbital localization',
        'Local MP2',
        'Local CC',
    ]

    cc_local_settings = [
        setting
        for setting in (cc.numerical_settings or [])
        if isinstance(setting, LocalCorrelationSettings)
    ]
    assert len(cc_local_settings) == 1
    cc_basis = [
        setting for setting in (cc.numerical_settings or []) if isinstance(setting, BasisSetContainer)
    ]
    assert len(cc_basis) == 1
    cc_components = {
        component.role: component for component in cc_basis[0].basis_set_components
    }
    assert cc_components['orbital'].basis_set == 'cc-pVDZ'
    assert cc_components['orbital'].n_total_basis_functions == 1320
    assert cc_components['auxiliary_post_hf'].basis_set == 'cc-pVDZ/C'
    assert cc_components['auxiliary_post_hf'].n_total_basis_functions == 3600

    thresholds = {
        threshold.name: threshold for threshold in cc_local_settings[0].screening_thresholds
    }
    assert thresholds['TCutPairs'].value == pytest.approx(1.0e-6)
    assert thresholds['TCutPairs'].applies_to == 'pair_screening'
    assert thresholds['TCutPNO'].value == pytest.approx(1.0e-8)
    assert thresholds['TCutPNO'].applies_to == 'orbital'
    assert thresholds['TCutPNOSingles'].value == pytest.approx(3.0e-10)
    assert thresholds['TCutPNOSingles'].applies_to == 'orbital'
