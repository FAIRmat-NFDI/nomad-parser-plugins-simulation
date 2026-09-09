import numpy as np
import pytest

from tests.parsers.common import (
    SimulationParserTestSuite,
    WorkflowTestSuite,
    approx,
    assert_approx,
)


class FHIAimsParserIntegrationSuite(SimulationParserTestSuite, WorkflowTestSuite):
    expected_program_name = 'FHI-aims'
    require_lattice_vectors = True
    require_periodic_boundary_conditions = True
    required_simulation_sections = None


class TestFeScfSpinPolarized(FHIAimsParserIntegrationSuite):
    archive_fixture = 'fe_scf_spinpol_archive'
    workflow_name = 'SinglePoint'

    @pytest.mark.integration
    def test_method_and_eigenvalue_mapping(self, archive):
        method = archive.data.model_method[0]
        assert method.m_def.name == 'DFT'
        assert list(method.numerical_settings[0].k_mesh[0].grid) == [16, 16, 16]
        assert method.numerical_settings[0].k_mesh[0].points.shape == (4, 3)
        np.testing.assert_allclose(
            method.numerical_settings[0].k_mesh[0].points,
            [[0, 0, 0], [0, 0, 0.0625], [0, 0, 0.125], [0, 0, 0.1875]],
        )
        assert method.xc.functional_key == 'LDA'

        system = archive.data.model_system[0]
        assert [state.chemical_symbol for state in system.particle_states] == ['Fe']
        np.testing.assert_allclose(
            system.lattice_vectors.to('angstrom').magnitude,
            [[-1.4, 1.4, 1.4], [1.4, -1.4, 1.4], [1.4, 1.4, -1.4]],
        )

        eigenvalue_sections = archive.data.outputs[0].electronic_eigenvalues
        assert len(eigenvalue_sections) == 2
        assert [section.spin_channel for section in eigenvalue_sections] == [0, 1]
        assert all(section.value.shape == (4, 19) for section in eigenvalue_sections)
        assert all(
            section.value.shape == section.occupation.shape
            for section in eigenvalue_sections
        )
        samples = [
            [
                ((0, 0), -258.846579),
                ((1, 1), -30.192224),
                ((2, 15), 0.509907),
                ((3, 18), 0.661682),
            ],
            [
                ((0, 0), -258.845593),
                ((1, 1), -30.148996),
                ((2, 15), 0.538254),
                ((3, 18), 0.697219),
            ],
        ]
        for section, section_samples in zip(eigenvalue_sections, samples, strict=True):
            values = section.value.to('hartree').magnitude
            for (kpoint, state), expected in section_samples:
                assert values[kpoint, state] == approx(expected, abs=1e-6)
        np.testing.assert_allclose(
            eigenvalue_sections[0].occupation[0][:5], [1, 1, 1, 1, 1]
        )


class TestFeBandSpinPolarized(FHIAimsParserIntegrationSuite):
    archive_fixture = 'fe_band_spinpol_archive'
    workflow_name = 'SinglePoint'

    @pytest.mark.integration
    def test_band_and_dos_mapping(self, archive):
        output = archive.data.outputs[0]
        assert len(output.electronic_band_structures) == 2
        assert [band.spin_channel for band in output.electronic_band_structures] == [
            0,
            1,
        ]
        assert all(
            band.value.shape == (1, 19)
            for band in output.electronic_band_structures
        )
        method = archive.data.model_method[0]
        assert_approx(
            method.numerical_settings[0].k_mesh[0].points, [[0, 0, 0]]
        )
        for band, expected_values, expected_occupations in zip(
            output.electronic_band_structures,
            [
                [-258.846579, -30.192224, -25.772599],
                [-258.845593, -30.148996, -25.738946],
            ],
            [[1, 1, 1], [1, 1, 1]],
            strict=True,
        ):
            assert_approx(
                band.value.to('hartree').magnitude[0, :3], expected_values, atol=1e-6
            )
            assert_approx(
                band.occupation[0, :3], expected_occupations
            )
        assert len(output.electronic_dos) == 2
        assert all(dos.value.shape == (50,) for dos in output.electronic_dos)
        assert all(dos.energies.points.shape == (50,) for dos in output.electronic_dos)
        for dos, expected_values in zip(
            output.electronic_dos,
            [
                [0.09655095, 0.20694765, 0.28868676],
                [0.10160034, 0.08330396, 0.06971712],
            ],
            strict=True,
        ):
            assert_approx(dos.value[:3].to('1 / eV').magnitude, expected_values, atol=1e-6)
            assert_approx(
                dos.energies.points[:3].to('eV').magnitude,
                [-15.0, -14.83673469, -14.67346939],
                atol=1e-6,
            )


class TestClNaDos(FHIAimsParserIntegrationSuite):
    archive_fixture = 'cl_na_dos_archive'
    workflow_name = 'SinglePoint'

    @pytest.mark.integration
    def test_dos_mapping(self, archive):
        method = archive.data.model_method[0]
        assert list(method.numerical_settings[0].k_mesh[0].grid) == [10, 10, 10]
        assert_approx(method.numerical_settings[0].k_mesh[0].points, [[0, 0, 0]])

        dos = archive.data.outputs[0].electronic_dos
        assert len(dos) == 1
        assert dos[0].value.shape == (50,)
        assert dos[0].energies.points.shape == (50,)
        assert_approx(
            dos[0].value[:3].to('1/eV').magnitude,
            [0.00233484, 0.24945582, 1.26529667],
            atol=1e-8,
        )
        assert_approx(
            dos[0].energies.points[:3].to('eV').magnitude,
            [-11.0, -10.71428571, -10.42857143],
            atol=1e-8,
        )


class TestSiGeometryOptimization(FHIAimsParserIntegrationSuite):
    archive_fixture = 'si_geomopt_archive'
    workflow_name = 'GeometryOptimization'
    required_simulation_sections = ('model_method', 'model_system', 'outputs')

    @pytest.mark.integration
    def test_archive_contract(self, archive):
        simulation = archive.data

        assert len(simulation.model_system) == 6
        assert len(simulation.outputs) == 5
        assert [
            state.chemical_symbol
            for state in simulation.model_system[0].particle_states
        ] == ['Si', 'Si']
        assert simulation.model_system[0].positions.shape == (2, 3)
        assert simulation.model_system[0].lattice_vectors.shape == (3, 3)
        np.testing.assert_allclose(
            simulation.model_system[0].lattice_vectors.to('angstrom').magnitude,
            [[2.7, 2.7, 0], [2.7, 0, 2.7], [0, 2.7, 2.7]],
        )
        np.testing.assert_allclose(
            simulation.model_system[0].positions.to('angstrom').magnitude,
            [[0, 0, 0], [1.35, 1.35, 1.35]],
        )
        assert archive.data.outputs[0].scf_steps.delta_energies_total[0].to(
            'eV'
        ).magnitude == approx(0.07631)
        assert archive.workflow2.method.optimization_method.strip() == 'Modified BFGS'
        assert archive.workflow2.method.convergence_targets[0].threshold.to(
            'eV/angstrom'
        ).magnitude == approx(0.01)

    @pytest.mark.integration
    def test_workflow_contract(self, archive):
        workflow = archive.workflow2

        assert workflow is not None
        assert workflow.m_def.name == 'GeometryOptimization'
        assert workflow.method is not None
        assert workflow.method.m_def.name == 'GeometryOptimizationMethod'
        assert workflow.method.optimization_method.strip() == 'Modified BFGS'

        force_targets = workflow.method.convergence_targets
        assert len(force_targets) == 1
        assert force_targets[0].m_def.name == 'ForceConvergenceTarget'
        assert force_targets[0].threshold_type == 'maximum'
        assert force_targets[0].threshold.to('eV/angstrom').magnitude == approx(
            0.01
        )

        energy_targets = workflow.method.single_point_convergence_targets
        assert len(energy_targets) == 1
        assert energy_targets[0].m_def.name == 'EnergyConvergenceTarget'
        assert energy_targets[0].threshold_type == 'absolute'
        assert energy_targets[0].threshold.to('eV').magnitude == approx(1e-6)

    @pytest.mark.integration
    def test_scf_steps_contract(self, archive):
        outputs = archive.data.outputs
        expected_n_scf = [11, 7, 8, 6, 7]

        assert len(outputs) == len(expected_n_scf)
        for output, n_scf in zip(outputs, expected_n_scf):
            scf_steps = output.scf_steps
            assert scf_steps is not None
            assert len(scf_steps.delta_energies_total) == n_scf
            assert len(scf_steps.delta_charge_abs) == n_scf
            assert len(scf_steps.durations) == n_scf

        first_steps = outputs[1].scf_steps
        assert first_steps.delta_energies_total[-1].to('eV').magnitude == approx(
            7.477e-09
        )
        assert (
            first_steps.delta_charge_abs[-1].to('coulomb').magnitude
            == approx(6.375e-08 * 1.602176634e-19)
        )

    @pytest.mark.integration
    def test_method_and_k_mesh_contract(self, archive):
        methods = archive.data.model_method

        assert len(methods) == 1
        method = methods[0]
        assert method.m_def.name == 'DFT'
        assert method.numerical_settings is not None

        k_spaces = [
            setting
            for setting in method.numerical_settings
            if setting.m_def.name == 'KSpace'
        ]
        assert len(k_spaces) == 1
        assert len(k_spaces[0].k_mesh) == 1

        k_mesh = k_spaces[0].k_mesh[0]
        assert k_mesh.m_def.name == 'KMesh'
        assert list(k_mesh.grid) == [8, 8, 8]
        assert list(k_mesh.offset) == approx([0.0, 0.0, 0.0])

    @pytest.mark.integration
    def test_scf_convergence_criteria_contract(self, archive):
        method = archive.data.model_method[0]
        criteria = [
            setting
            for setting in method.numerical_settings
            if setting.m_def.name == 'SelfConsistency'
        ]

        assert len(criteria) == 3
        by_name = {criterion.name: criterion for criterion in criteria}

        energy = by_name['total_energy_change']
        assert energy.threshold_change.to('eV').magnitude == approx(1e-6)

        density = by_name['charge_density_change']
        assert density.threshold_change == approx(1e-5)

        eigenvalues = by_name['sum_eigenvalues_change']
        assert eigenvalues.threshold_change.to('eV').magnitude == approx(1e-3)

    @pytest.mark.integration
    def test_eigenvalue_and_band_structure_contract(self, archive):
        expected_values = [
            [-65.168926, -65.168926, -5.071833, -5.071497, -3.510756],
            [-65.169070, -65.169070, -5.071843, -5.071511, -3.510774],
            [-65.169347, -65.169347, -5.071866, -5.071542, -3.510813],
            [-65.169364, -65.169364, -5.071867, -5.071544, -3.510815],
            [-65.169435, -65.169435, -5.071873, -5.071552, -3.510823],
        ]
        for output, expected in zip(archive.data.outputs, expected_values, strict=True):
            eigenvalues = output.electronic_eigenvalues
            band_structures = output.electronic_band_structures

            assert len(eigenvalues) == 1
            assert len(band_structures) == 1

            eigenvalue_section = eigenvalues[0]
            band_structure_section = band_structures[0]
            assert eigenvalue_section.spin_channel is None
            assert band_structure_section.spin_channel is None
            assert (
                eigenvalue_section.value.shape == eigenvalue_section.occupation.shape
            )
            assert (
                band_structure_section.value.shape
                == band_structure_section.occupation.shape
            )
            assert eigenvalue_section.value.shape[0] == 1
            assert_approx(
                eigenvalue_section.value.to('hartree').magnitude[0][:5],
                expected,
                atol=1e-6,
            )
            assert_approx(
                band_structure_section.value.to('hartree').magnitude[0][:5],
                expected,
                atol=1e-6,
            )
            assert eigenvalue_section.occupation[0][0] == approx(2.0)


class TestH2OMolecularDynamics(FHIAimsParserIntegrationSuite):
    archive_fixture = 'ho_md_archive'
    workflow_name = 'MolecularDynamics'

    @pytest.mark.integration
    def test_trajectory_mapping(self, archive):
        assert len(archive.data.model_system) == 6
        assert len(archive.data.outputs) == 5
        assert archive.data.model_system[2].velocities is not None
        assert archive.data.outputs[4].total_forces[0].value.shape == (3, 3)
        assert_approx(
            archive.data.model_system[0].positions.to('angstrom').magnitude,
            [
                [2.8, 2.6, 2.0],
                [1.2, 2.6, 2.0],
                [2.0, 2.0, 2.0],
            ],
        )
        assert_approx(
            archive.data.model_system[2].velocities.to(
                'angstrom / femtosecond'
            ).magnitude[0],
            [-28.17773990, 20.13038952, -13.61884175],
            atol=1e-8,
        )
        assert_approx(
            archive.data.outputs[1].total_energies[0].value.to('eV').magnitude,
            -2078.48223446601,
            atol=1e-6,
        )
        assert archive.workflow2.method.n_steps == 5
        assert archive.workflow2.method.integration_timestep.to('ps').magnitude == approx(
            0.001
        )
        assert archive.workflow2.results.n_steps == 5
        assert archive.workflow2.results.finished_normally is True


class TestGaAsHybrid(FHIAimsParserIntegrationSuite):
    archive_fixture = 'gaas_hse06_soc_archive'
    workflow_name = 'SinglePoint'

    @pytest.mark.integration
    def test_band_mapping(self, archive):
        method = archive.data.model_method[0]
        assert list(method.numerical_settings[0].k_mesh[0].grid) == [12, 12, 12]
        assert archive.data.outputs[0].electronic_band_structures[0].value.shape == (
            1,
            43,
        )


class TestCeO2DFTU(FHIAimsParserIntegrationSuite):
    archive_fixture = 'ceo2_dftu_archive'
    workflow_name = 'SinglePoint'

    @pytest.mark.integration
    def test_method_and_k_mesh_mapping(self, archive):
        method = archive.data.model_method[0]
        assert method.xc.functional_key == 'PBEsol'
        assert list(method.numerical_settings[0].k_mesh[0].grid) == [1, 1, 1]


class TestHeGW(FHIAimsParserIntegrationSuite):
    archive_fixture = 'he_gw_archive'
    workflow_name = 'SinglePoint'


class TestSiliconV071914(FHIAimsParserIntegrationSuite):
    archive_fixture = 'silicon_v071914_archive'
    workflow_name = 'SinglePoint'


class TestSiliconV171221(FHIAimsParserIntegrationSuite):
    archive_fixture = 'silicon_v171221_archive'
    workflow_name = 'SinglePoint'


class TestCHNGW(FHIAimsParserIntegrationSuite):
    archive_fixture = 'chn_gw_archive'
    workflow_name = 'SinglePoint'


class TestSiliconGWBands(FHIAimsParserIntegrationSuite):
    archive_fixture = 'si_gw_bands_archive'
    workflow_name = 'SinglePoint'

    @pytest.mark.integration
    def test_atoms_k_mesh_band_values_and_occupations(self, archive):
        system = archive.data.model_system[0]
        assert [state.chemical_symbol for state in system.particle_states] == [
            'Si',
            'Si',
        ]
        assert_approx(
            system.positions.to('angstrom').magnitude,
            [
                [0.0, 0.0, 0.0],
                [1.36203139, 1.36203140, 1.36206560],
            ],
            atol=1e-8,
        )
        assert_approx(
            system.lattice_vectors.to('angstrom').magnitude,
            [
                [0.00013684, 2.72413122, 2.72426803],
                [2.72413117, 0.00013679, 2.72426800],
                [2.72385756, 2.72385758, -0.00027363],
            ],
            atol=1e-8,
        )

        method = archive.data.model_method[0]
        k_mesh = method.numerical_settings[0].k_mesh[0]
        assert list(k_mesh.grid) == [8, 8, 8]
        assert_approx(k_mesh.points, [[0, 0, 0]])

        band = archive.data.outputs[0].electronic_band_structures[0]
        assert band.value.shape == (1, 78)
        assert_approx(
            band.value.to('hartree').magnitude[0, :6],
            [-65.782015, -65.782015, -5.139790, -5.139521, -3.523694, -3.523694],
            atol=1e-6,
        )
        assert_approx(band.occupation[0, :16], [2] * 14 + [0, 0])
        assert_approx(band.occupation[0, 14:], [0] * 64)
        assert len(archive.data.outputs[0].electronic_band_gaps) == 1
        assert archive.data.outputs[0].electronic_band_gaps[0].spin_channel is None
        assert archive.data.outputs[0].electronic_band_gaps[0].value.to(
            'hartree'
        ).magnitude == approx(0.093853)


class TestNativeTight(FHIAimsParserIntegrationSuite):
    archive_fixture = 'native_tight_archive'
    workflow_name = 'GeometryOptimization'


class TestNativeIntermediate(FHIAimsParserIntegrationSuite):
    archive_fixture = 'native_intermediate_archive'
    workflow_name = 'GeometryOptimization'


class TestNativeLightSpd(FHIAimsParserIntegrationSuite):
    archive_fixture = 'native_light_spd_archive'
    workflow_name = 'GeometryOptimization'
