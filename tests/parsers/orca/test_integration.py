from contextlib import nullcontext

import numpy as np
import pytest
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
from nomad_simulations.schema_packages.numerical_settings import (
    LocalCorrelationSettings,
)
from nomad_simulations.schema_packages.properties.molecular_orbitals import (
    MolecularOrbitals,
)

from tests.parsers.common import SimulationParserTestSuite, WorkflowTestSuite


class OrcaParserIntegrationSuite(SimulationParserTestSuite):
    """ORCA-specific configuration of the shared integration suite (molecular,
    so lattice vectors and periodic boundary conditions are not required)."""

    expected_program_name = 'ORCA'


class TestRIMP2Water(OrcaParserIntegrationSuite):
    archive_fixture = 'ri_mp2_water_archive'

    @pytest.mark.integration
    def test_mp2_method_and_basis_sets(self, archive):
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
        assert [
            (component.type, component.n_total_basis_functions)
            for component in basis_sets[0].basis_set_components
        ] == [('GTO', 24), ('GTO', 76)]

    @pytest.mark.integration
    def test_model_system_and_molecular_orbitals(self, archive):
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
        assert system.total_spin_multiplicity == 1

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


@pytest.mark.large_fixture
class TestCASCIQuantumDots(OrcaParserIntegrationSuite):
    # The CoPc CASCI/NEVPT2 mainfile is ~19 MB, too large for the default PR
    # gate; `large_fixture` moves this class to the nightly gate.
    archive_fixture = 'casci_qd_archive'
    # This CASCI/NEVPT2 run carries no parsed output block.
    required_simulation_sections = ('model_method', 'model_system')

    @pytest.mark.integration
    def test_multireference_methods_basis_sets_and_relativity(self, archive):
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
        assert casci.active_space.orbital_space_type == 'CAS'
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


class TestDLPNOCoupledCluster(OrcaParserIntegrationSuite, WorkflowTestSuite):
    archive_fixture = 'dlpno_cc_archive'
    workflow_name = 'SerialWorkflow'

    @pytest.mark.integration
    def test_coupled_cluster_methods(self, archive):
        hf = next(
            method for method in archive.data.model_method if isinstance(method, HF)
        )
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
        cc = next(
            method for method in archive.data.model_method if isinstance(method, CC)
        )

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

        mp2_local_settings = next(
            settings
            for settings in mp2.numerical_settings
            if isinstance(settings, LocalCorrelationSettings)
        )
        assert len(mp2_local_settings.screening_thresholds) == 8
        assert mp2_local_settings.screening_thresholds[0].name == 'TCutPairs'
        assert mp2_local_settings.screening_thresholds[0].value == pytest.approx(1e-6)

        cc_local_settings = next(
            settings
            for settings in cc.numerical_settings
            if isinstance(settings, LocalCorrelationSettings)
        )
        assert len(cc_local_settings.screening_thresholds) == 13
        assert cc_local_settings.screening_thresholds[-1].name == 'TCutDOWeak'
        assert cc_local_settings.screening_thresholds[-1].value == pytest.approx(4e-3)


class TestDFTPrintMOs(OrcaParserIntegrationSuite):
    archive_fixture = 'dft_mos_archive'

    @pytest.mark.integration
    def test_archive_serialization_round_trip(self, archive):
        # The MO coefficients are HDF5-backed here; skip the generic round-trip,
        # which detaches from the upload-backed context.
        pytest.skip('HDF5-backed coefficients are not round-tripped in isolation')

    @pytest.mark.integration
    def test_model_system_serialization_round_trip(self, archive):
        pytest.skip('HDF5-backed coefficients are not round-tripped in isolation')

    @pytest.mark.integration
    def test_dft_method(self, archive):
        dft = next(
            method for method in archive.data.model_method if isinstance(method, DFT)
        )
        assert dft.reference_form == 'RKS'
        assert dft.xc.functional_key == 'B88+LYP'
        assert dft.xc.global_exact_exchange == pytest.approx(0.2)

    @pytest.mark.integration
    def test_molecular_orbital_coefficients(self, archive):
        molecular_orbitals = archive.data.outputs[0].molecular_orbitals[0]
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
