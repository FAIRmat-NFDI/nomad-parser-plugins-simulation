"""Layer 3 mapping-contract tests.

Use small hand-built inputs to test the ORCA transformers directly.
"""

import pytest

from nomad_simulation_parsers.parsers.orca.parser import OutParser


@pytest.mark.unit
class TestDFTMapping:
    def test_maps_functional_key_and_reference_form(self):
        method = OutParser().get_dft(
            {
                'exchange_functional': 'B88',
                'correlation_functional': 'LYP',
                'fraction_hf_exchange': 0.2,
                'hf_type': 'RKS',
            }
        )
        assert method['xc']['functional_key'] == 'B88+LYP'
        assert method['xc']['global_exact_exchange'] == pytest.approx(0.2)
        assert method['reference_form'] == 'RKS'

    def test_absent_xc_yields_no_method(self):
        # absent (not empty-with-zeros): an SCF block with no functional is not DFT.
        assert OutParser().get_dft({}) == {}


@pytest.mark.unit
class TestCoupledClusterMapping:
    def test_maps_ccsd_t_with_dlpno(self):
        # The reader captures `coupled_cluster_type` as [A-Z]+ (regex stops at the
        # first paren), so a real CCSD(T) run yields base type `CCSD` plus the
        # perturbative-triples flag. The (T) is represented separately, not folded
        # into `type`/`excitation_order`.
        (method,) = OutParser().get_coupled_cluster_methods(
            {
                'coupled_cluster_type': 'CCSD',
                'perturbative_triple_excitations_on_off': 'ON',
                'kc_formation': 'DLPNO',
            },
            'input',
        )
        assert method['type'] == 'CCSD'
        assert method['excitation_order'] == [1, 2]  # ordering matters
        assert method['perturbative_correction'] == '(T)'
        assert method['perturbative_correction_order'] == [3]
        assert method['local_correlation']['type'] == 'DLPNO'

    def test_no_cc_type_yields_empty(self):
        assert OutParser().get_coupled_cluster_methods({}, '') == []


@pytest.mark.unit
class TestPerturbationMapping:
    def test_maps_scs_mp2_and_dlpno_from_input_hint(self):
        # local-correlation type is inferred from the input-file hint, not cc_data.
        (method,) = OutParser().get_perturbation_methods(
            {'spin_component_scaling': 'SCS'}, {}, 'some DLPNO-MP2 run'
        )
        assert method['type'] == 'MP'
        assert method['order'] == 2
        assert method['spin_component_scaling'] == 'SCS'
        assert method['local_correlation']['type'] == 'DLPNO'


@pytest.mark.unit
class TestMultireferenceMapping:
    def _casscf(self):
        return {
            'n_active_electrons': 13,
            'n_active_orbitals': 8,
            'block': [{'multiplicity': 4, 'root_weights': [1.0]}],
        }

    def test_casscf_active_space_and_states(self):
        (method,) = OutParser().get_multireference_scf_methods(self._casscf(), '')
        assert method['type'] == 'CASSCF'
        assert method['active_space'] == {
            'n_active_electrons': 13,
            'n_active_orbitals': 8,
            'orbital_space_type': 'CAS',
        }
        assert method['state_multiplicities'] == [4]

    def test_casci_discriminated_by_maxiter_marker(self):
        # `MAXITER 1` in the input file distinguishes CASCI from CASSCF.
        (method,) = OutParser().get_multireference_ci_methods(
            self._casscf(), 'nel MAXITER 1'
        )
        assert method['type'] == 'CASCI'

    def test_nevpt2_name_assembly(self):
        casscf = {**self._casscf(), 'qd_nevpt_type': 'QD', 'pt_method': 'SC_NEVPT2'}
        (method,) = OutParser().get_multireference_pt_methods(casscf, 'NEVPT2')
        assert method['type'] == 'NEVPT'
        assert method['order'] == 2
        assert method['name'] == 'QD-SC-NEVPT2'

    def test_empty_casscf_yields_empty(self):
        assert OutParser().get_multireference_scf_methods({}, '') == []


@pytest.mark.unit
class TestRelativityMapping:
    def test_maps_dkh_with_order(self):
        model = OutParser().get_relativity_model(
            {'relativistic_hamiltonian': {'method': 'DKH', 'dkh_order': 2}}
        )
        assert model == {'level': 'scalar', 'approximation': 'DKH', 'dkh_order': 2}

    def test_absent_relativity_yields_none(self):
        assert OutParser().get_relativity_model({}) is None


@pytest.mark.unit
class TestBasisSetMapping:
    _SOURCE = {
        'basis_set': {
            'basis_set_name': {
                'main_basis_set': 'def2-SVP',
                'auxc_basis_set': 'def2-SVP/C',
                'auxj_basis_set': 'def2/J',
            },
            'basis_set_total': {'main_basis_set': 24, 'auxc_basis_set': 76},
        }
    }

    def test_roles_depend_on_method_post_hf(self):
        # A post-HF method (CC) admits the correlation-fitting aux basis, not auxj.
        parser = OutParser()
        parser._method = 'CC'
        components = parser.get_basis_set_components(self._SOURCE)
        assert [(c['basis_set'], c['role']) for c in components] == [
            ('def2-SVP', 'orbital'),
            ('def2-SVP/C', 'auxiliary_post_hf'),
        ]

    def test_roles_depend_on_method_scf(self):
        # An SCF method (DFT) admits the Coulomb-fitting aux basis, not auxc.
        parser = OutParser()
        parser._method = 'DFT'
        components = parser.get_basis_set_components(self._SOURCE)
        assert [(c['basis_set'], c['role']) for c in components] == [
            ('def2-SVP', 'orbital'),
            ('def2/J', 'auxiliary_scf'),
        ]


@pytest.mark.unit
class TestSCFStepsMapping:
    @pytest.mark.parametrize(
        ('headers', 'rows', 'durations'),
        [
            pytest.param(
                [
                    'Energy Delta-E Max-DP RMS-DP [F,P] Damp',
                    'Energy Delta-E Grad Rot Max-DP RMS-DP',
                ],
                [
                    [-76.28, 0.0, 0.04, 0.005, 0.29, 0.85],
                    [-76.31, -0.03, 0.001, 0.002, 0.003, 0.0004],
                ],
                None,
                id='orca4',
            ),
            pytest.param(
                [
                    'Energy (Eh) Delta-E RMSDP MaxDP DIISErr Damp Time(sec)',
                    'Energy (Eh) Delta-E RMSDP MaxDP MaxGrad Time(sec)',
                ],
                [
                    [-76.28, 0.0, 0.005, 0.04, 0.29, 0.85, 0.0],
                    [-76.31, -0.03, 0.0004, 0.003, 0.001, 0.2],
                ],
                [0.0, 0.2],
                id='orca6',
            ),
        ],
    )
    def test_maps_scf_columns_across_iteration_blocks(self, headers, rows, durations):
        blocks = [
            {'header': header.split(), 'values': [row]}
            for header, row in zip(headers, rows, strict=True)
        ]
        steps = OutParser().get_scf_iterations(
            {'self_consistent': {'scf_iterations': blocks}}
        )
        assert [value.magnitude for value in steps['energy']] == pytest.approx(
            [-76.28, -76.31]
        )
        assert [value.magnitude for value in steps['deltae']] == pytest.approx(
            [0.0, -0.03]
        )
        assert steps['maxdp'] == pytest.approx([0.04, 0.003])
        assert steps['rmsdp'] == pytest.approx([0.005, 0.0004])
        if durations is None:
            assert 'time' not in steps
        else:
            assert [value.magnitude for value in steps['time']] == pytest.approx(
                durations
            )


@pytest.mark.unit
class TestGeometryOptimizationMapping:
    @pytest.mark.parametrize(
        ('banner', 'is_converged'),
        [
            ({'is_converged': ['THE', 'OPTIMIZATION', 'HAS', 'CONVERGED']}, True),
            (
                {'is_not_converged': ['The', 'optimization', 'did', 'not', 'converge']},
                False,
            ),
            ({}, None),
        ],
        ids=['converged', 'unconverged', 'absent'],
    )
    def test_maps_workflow_method_and_last_cycle_results(self, banner, is_converged):
        source = {
            'geometry_optimization': {
                'update_method': ['Update', 'BFGS'],
                'cycle': [
                    {
                        'geom_opt_max_gradient': 0.04,
                        'geom_opt_max_displacement': 0.15,
                    },
                    {
                        'geom_opt_max_gradient': 0.004,
                        'geom_opt_max_displacement': 0.01,
                    },
                ],
                **banner,
            }
        }
        parser = OutParser()
        assert parser.get_geometry_optimization_method(source) == {
            'optimization_type': 'atomic',
            'optimization_method': 'BFGS',
            'sampling_frequency': 1,
        }
        results = parser.get_geometry_optimization_results(source)
        assert results['is_converged'] is is_converged
        assert results['n_steps'] == 2
        assert results['steps'] == [0, 1]
        assert results['final_force_maximum'] == pytest.approx(0.004)
        assert results['final_displacement_maximum'] == pytest.approx(0.01)
