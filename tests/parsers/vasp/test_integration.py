import numpy as np
import pytest
from nomad.datamodel import EntryArchive
from nomad.utils import get_logger

from tests.parsers.common import (
    SimulationParserTestSuite,
    WorkflowTestSuite,
    assert_approx,
)

LOGGER = get_logger(__name__)


class VASPIntegrationSuite(SimulationParserTestSuite, WorkflowTestSuite):
    expected_program_name = 'vasp'
    require_lattice_vectors = True
    require_periodic_boundary_conditions = True
    require_scf_steps = True
    require_convergence_targets = True

    @pytest.fixture(scope='class')
    @classmethod
    def archive(cls, request):
        source = request.getfixturevalue(cls.archive_fixture)
        return EntryArchive.m_from_dict(source.m_to_dict())

    @pytest.mark.integration
    def test_outputs_have_scf_steps(self, archive):
        if not self.require_scf_steps:
            pytest.skip('SCF-step contract is not required for this parser case')
        assert archive.data.outputs
        scf_outputs = [output for output in archive.data.outputs if output.scf_steps]
        assert scf_outputs
        assert all(
            output.scf_steps.energies_total is not None
            and len(output.scf_steps.energies_total) > 0
            for output in scf_outputs
        )

    @pytest.mark.integration
    def test_workflow_has_convergence_targets(self, archive):
        if not self.require_convergence_targets:
            pytest.skip(
                'convergence-target contract is not required for this parser case'
            )
        assert archive.workflow2 is not None
        assert archive.workflow2.method is not None
        assert archive.workflow2.method.convergence_targets



class TestAgacRelaxOutcarArchive(VASPIntegrationSuite):
    archive_fixture = 'agac_relax_outcar_archive'
    workflow_name = 'SinglePoint'

    @pytest.mark.integration
    def test_archive_contract(self, archive):
        simulation = archive.data
        assert len(simulation.model_system) == 1
        assert simulation.model_system[0].positions.shape == (2, 3)
        symbols = [
            state.chemical_symbol
            for state in simulation.model_system[0].particle_states
        ]
        assert symbols == [
            'Ac',
            'Ag',
        ]
        assert len(simulation.outputs) == 2
        assert len(simulation.outputs[0].scf_steps.energies_total) == 11
        assert simulation.outputs[-1].electronic_dos

    @pytest.mark.integration
    def test_electronic_outputs_and_convergence(self, archive):
        output = archive.data.outputs[0]
        assert output.electronic_eigenvalues
        assert output.electronic_eigenvalues[0].value is not None
        assert output.electronic_eigenvalues[0].occupation is not None
        if output.electronic_band_gaps:
            assert output.electronic_band_gaps[0].value is not None

        dos = archive.data.outputs[-1].electronic_dos[0]
        assert dos.value is not None
        assert dos.energies is not None
        assert dos.energies.points is not None
        assert dos.projected_dos is not None
        assert len(dos.projected_dos) == 32
        assert dos.projected_dos[0].value.to('1/eV').magnitude[218] == pytest.approx(
            0.07835
        )
        assert dos.projected_dos[15].value.to('1/eV').magnitude[277] == pytest.approx(
            0.20490
        )
        assert dos.projected_dos[16].value.to('1/eV').magnitude[238] == pytest.approx(
            0.33900
        )


    @pytest.mark.integration
    def test_workflow_convergence(self, archive):
        workflow = archive.workflow2
        target = workflow.method.convergence_targets[0]
        assert target.m_def.name == 'EnergyConvergenceTarget'
        assert target.threshold_type == 'absolute'
        assert target.threshold.to('eV').magnitude == pytest.approx(1.0e-4)
        scf_steps = archive.data.outputs[0].scf_steps
        assert len(scf_steps.delta_energies_total) == 10
        assert len(scf_steps.durations) == 11
        assert scf_steps.energies_total[-1].to('eV').magnitude == pytest.approx(
            -6.97148118
        )
        assert np.isfinite(scf_steps.delta_energies_total[-1].to('eV').magnitude)

    @pytest.mark.integration
    def test_xc_functional_and_jacobs_ladder(self, archive):
        dft = archive.data.model_method[0]
        assert dft.xc is not None
        assert dft.xc.functional_key == 'PBE'
        dft.normalize(archive, LOGGER)
        assert {component.canonical_label for component in dft.xc.components} == {
            'XC_GGA_X_PBE',
            'XC_GGA_C_PBE',
        }
        assert dft.jacobs_ladder == 'GGA'


class TestAgacRelaxVasprunArchive(VASPIntegrationSuite):
    archive_fixture = 'agac_relax_vasprun_archive'
    workflow_name = 'GeometryOptimization'

    @pytest.mark.integration
    def test_archive_contract(self, archive):
        simulation = archive.data
        assert len(simulation.model_system) == 3
        assert len(simulation.outputs) == 3
        assert simulation.model_system[0].positions.shape == (2, 3)
        assert simulation.outputs[0].scf_steps is not None
        assert all(output.scf_steps for output in simulation.outputs)
        assert np.isfinite(
            simulation.outputs[-1].scf_steps.energies_total[-1].to('eV').magnitude
        )

    @pytest.mark.integration
    def test_electronic_outputs_and_convergence(self, archive):
        output = next(
            output for output in archive.data.outputs if output.electronic_eigenvalues
        )
        assert output.electronic_eigenvalues[0].occupation is not None
        assert output.electronic_band_gaps
        assert output.electronic_dos
        assert output.electronic_dos[0].energies.points is not None


    @pytest.mark.integration
    def test_workflow_convergence(self, archive):
        workflow = archive.workflow2
        targets = workflow.method.convergence_targets
        assert len(targets) == 1
        assert targets[0].m_def.name == 'EnergyConvergenceTarget'
        assert targets[0].threshold.to('eV').magnitude == pytest.approx(1.0e-3)
        sp_targets = workflow.method.single_point_convergence_targets
        assert len(sp_targets) == 1
        assert sp_targets[0].threshold.to('eV').magnitude == pytest.approx(1.0e-4)
        assert [
            len(output.scf_steps.energies_total) for output in archive.data.outputs
        ] == [
            12,
            10,
            6,
        ]

    @pytest.mark.integration
    def test_xc_functional_and_jacobs_ladder(self, archive):
        dft = archive.data.model_method[0]
        assert dft.xc is not None
        assert dft.xc.functional_key == 'PBE'
        dft.normalize(archive, LOGGER)
        assert {component.canonical_label for component in dft.xc.components} == {
            'XC_GGA_X_PBE',
            'XC_GGA_C_PBE',
        }
        assert dft.jacobs_ladder == 'GGA'


class TestWithCHGCARArchive:
    @pytest.mark.integration
    def test_charge_density_is_mapped(self, chgcar_archive):
        outputs = chgcar_archive.data.outputs
        assert len(outputs) == 14
        charge_density = outputs[-1].charge_density
        assert len(charge_density) == 1
        assert charge_density[0].value_h5_dataset is not None
        with charge_density[0].value_h5_dataset as dataset:
            dataset_array = dataset[:]
            assert np.shape(dataset_array) == (10, 10, 10)
            assert dataset_array[9][9][9] == pytest.approx(0.18013097030e-05)


class TestMgStaticArchive:
    @pytest.mark.integration
    def test_mg_static_maps_electronic_data(self, mg_static_archive):
        simulation = mg_static_archive.data
        assert len(simulation.model_system) == 1
        assert len(simulation.model_method) == 1
        assert len(simulation.outputs) == 1
        output = simulation.outputs[0]
        assert output.electronic_eigenvalues
        assert output.electronic_dos
        assert output.scf_steps


class TestMgBandsArchive:
    @pytest.mark.integration
    def test_mg_bands_maps_band_calculation(self, mg_bands_archive):
        output = mg_bands_archive.data.outputs[0]
        assert output.electronic_eigenvalues
        assert output.electronic_band_gaps
        assert output.scf_steps


class TestSiliconDosArchive:
    @pytest.mark.integration
    def test_silicon_dos_maps_dos(self, silicon_dos_archive):
        output = silicon_dos_archive.data.outputs[0]
        assert output.electronic_dos
        dos = output.electronic_dos[0]
        assert dos.energies is not None
        assert dos.value is not None


class TestSiliconBandArchive:
    @pytest.mark.integration
    def test_silicon_band_maps_band_gap(self, silicon_band_archive):
        output = silicon_band_archive.data.outputs[0]
        assert output.electronic_eigenvalues
        assert output.electronic_band_gaps
        assert output.electronic_band_gaps[0].value is not None


class BasicVASPIntegrationSuite(VASPIntegrationSuite):
    required_simulation_sections = ('model_method', 'outputs')
    require_lattice_vectors = False
    require_periodic_boundary_conditions = False
    require_scf_steps = False
    require_convergence_targets = False
    workflow_name = 'SinglePoint'


class TestSiliconGWArchive(BasicVASPIntegrationSuite):
    archive_fixture = 'silicon_gw_archive'

    @pytest.mark.integration
    def test_kspace(self, silicon_gw_kspace_archive):
        method = silicon_gw_kspace_archive.data.model_method[0]
        k_space = next(
            setting
            for setting in method.numerical_settings
            if setting.m_def.name == 'KSpace'
        )
        assert len(k_space.k_mesh) == 1
        k_mesh = k_space.k_mesh[0]
        assert list(k_mesh.grid) == [6, 6, 6]
        assert_approx(k_mesh.points[0], [0.0, 0.0, 0.0])
        assert_approx(k_mesh.points[-1], [-0.33333333, 0.5, 0.16666667])
        assert_approx(
            k_mesh.weights,
            [
                0.00462963,
                0.03703704,
                0.03703704,
                0.01851852,
                0.02777778,
                0.11111111,
                0.11111111,
                0.11111111,
                0.05555556,
                0.02777778,
                0.11111111,
                0.05555556,
                0.01388889,
                0.11111111,
                0.11111111,
                0.05555556,
            ],
        )
        assert_approx(k_mesh.offset, [0.0, 0.0, 0.0])


class TestGammaOutcarArchive(BasicVASPIntegrationSuite):
    archive_fixture = 'gamma_outcar_archive'
    workflow_name = 'GeometryOptimization'

    @pytest.mark.integration
    def test_program_and_kspace(self, gamma_outcar_kmesh_archive):
        simulation = gamma_outcar_kmesh_archive.data
        assert simulation.program.version == (
            '5.4.1 05Feb16 gamma-only parallel IFC91_ompi'
        )
        k_space = next(
            setting
            for setting in simulation.model_method[0].numerical_settings
            if setting.m_def.name == 'KSpace'
        )
        assert len(k_space.k_mesh) == 1
        assert_approx(k_space.k_mesh[0].points, [[0.0, 0.0, 0.0]])
        assert_approx(k_space.k_mesh[0].weights, [1.0])


class TestAlNVasprunArchive(BasicVASPIntegrationSuite):
    archive_fixture = 'al_n_vasprun_archive'
    required_simulation_sections = ('model_method',)

    @pytest.mark.integration
    def test_potcar_metadata(self, archive):
        method = archive.data.model_method[0]
        assert method.m_def.name == 'DFT'
        assert method.xc.functional_key == 'PBE'
        assert archive.data.program.version == '5.4.4.18Apr17-6-g9f103f2a35'


class TestAlNOutcarArchive(BasicVASPIntegrationSuite):
    archive_fixture = 'al_n_outcar_archive'
    required_simulation_sections = ('model_method',)

    @pytest.mark.integration
    def test_potcar_metadata(self, archive):
        method = archive.data.model_method[0]
        assert method.m_def.name == 'DFT'
        assert method.xc.functional_key == 'PBE'
        assert archive.data.program.version == '5.4.4.18 Apr17 complex parallel LINUX'


class TestFGWOutcarArchive(BasicVASPIntegrationSuite):
    archive_fixture = 'f_gw_outcar_archive'
    required_simulation_sections = ('model_method',)

    @pytest.mark.integration
    def test_potcar_method_type(self, archive):
        assert archive.data.model_method[0].m_def.name == 'DFT'


class TestHybridArchive(BasicVASPIntegrationSuite):
    archive_fixture = 'hybrid_archive'

    @pytest.mark.integration
    def test_kspace_and_xc_functional(self, archive):
        method = archive.data.model_method[0]
        assert method.xc.functional_key == 'HSE06'
        k_space = next(
            setting
            for setting in method.numerical_settings
            if setting.m_def.name == 'KSpace'
        )
        assert len(k_space.k_mesh) == 1
        k_mesh = k_space.k_mesh[0]
        assert list(k_mesh.grid) == [1, 1, 1]
        assert_approx(k_mesh.points, [[0.0, 0.0, 0.0]])
        assert_approx(k_mesh.weights, [1.0])
        assert_approx(k_mesh.offset, [0.0, 0.0, 0.0])


class TestMetaggaArchive(BasicVASPIntegrationSuite):
    archive_fixture = 'metagga_archive'

    @pytest.mark.integration
    def test_kspace_and_xc_functional(self, archive):
        method = archive.data.model_method[0]
        assert method.xc.functional_key == 'PBE'
        k_space = next(
            setting
            for setting in method.numerical_settings
            if setting.m_def.name == 'KSpace'
        )
        assert len(k_space.k_mesh) == 1
        k_mesh = k_space.k_mesh[0]
        assert list(k_mesh.grid) == [1, 1, 1]
        assert_approx(k_mesh.points, [[0.0, 0.0, 0.0]])
        assert_approx(k_mesh.weights, [1.0])
        assert_approx(k_mesh.offset, [0.0, 0.0, 0.0])


class TestDftuMultiParameterArchive(BasicVASPIntegrationSuite):
    archive_fixture = 'dftu_multi_parameter_archive'


class TestDftuSingleParameterArchive(BasicVASPIntegrationSuite):
    archive_fixture = 'dftu_single_parameter_archive'


class TestDftuMultiParameterNoIncarArchive(BasicVASPIntegrationSuite):
    archive_fixture = 'dftu_multi_parameter_no_incar_archive'


class TestMalformedTimeArchive(BasicVASPIntegrationSuite):
    archive_fixture = 'malformed_time_archive'


class TestBrokenVasprunArchive(BasicVASPIntegrationSuite):
    archive_fixture = 'broken_vasprun_archive'
    required_simulation_sections = ()


class TestBooleanDottedArchive(BasicVASPIntegrationSuite):
    archive_fixture = 'boolean_dotted_archive'


class TestBooleanSingleCharArchive(BasicVASPIntegrationSuite):
    archive_fixture = 'boolean_single_char_archive'


class TestBooleanMixedCaseArchive(BasicVASPIntegrationSuite):
    archive_fixture = 'boolean_mixed_case_archive'


class TestBooleanLowercaseArchive(BasicVASPIntegrationSuite):
    archive_fixture = 'boolean_lowercase_archive'
