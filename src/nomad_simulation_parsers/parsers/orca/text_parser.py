import numpy as np
from nomad.parsing.file_parser import Quantity, TextParser
from nomad.units import ureg


class OutReader(TextParser):
    def init_quantities(self):
        re_float = r'[-+]?\d+\.?\d*(?:[Ee][-+]\d+)?'
        re_n = r'\r*\n'
        self._energy_mapping = {
            'Total Energy': 'energy_total',
            'Nuclear Repulsion': 'energy_nuclear_repulsion',
            'Electronic Energy': 'elec_energy',
            'One Electron Energy': 'one_elec_energy',
            'Two Electron Energy': 'two_elec_energy',
            'Potential Energy': 'potential_energy',
            'Kinetic Energy': 'energy_kinetic_electronic',
            r'E\(X\)': 'energy_exchange',
            r'E\(C\)': 'energy_correlation',
            r'E\(XC\)': 'energy_XC',
        }

        self._timing_mapping = {
            'Total time': 'final_time',
            'Sum of individual times': 'sum_individual_times',
            'Fock matrix formation': 'fock_matrix_formation',
            'Coulomb formation': 'coulomb_formation',
            r'Split\-RI-J': 'split_rj',
            'XC integration': 'xc_integration',
            r'Basis function eval\.': 'basis_fn_evaluation',
            r'Density eval\.': 'density_evaluation',
            r'XC\-Functional eval\.': 'xc_functional_evaluation',
            r'XC\-Potential eval\.': 'potential_evaluation',
            'Diagonalization': 'diagonalization',
            'Density matrix formation': 'density_matrix_formation',
            'Population analysis': 'population_analysis',
            'Initial guess': 'initial_guess',
            'Orbital Transformation': 'orbital_transformation',
            'Orbital Orthonormalization': 'orbital_orthonormalization',
            'DIIS solution': 'diis_solution',
            'Grid generation': 'grid_generation',
            'Total SCF gradient time': 'scf_gradient',
        }

        # Coupled cluster related quantities
        coupled_cluster_quantities = [
            Quantity(
                'coupled_cluster_type',
                r'Correlation treatment\s+\.\.\.\s+([A-Z]+)',
                repeats=False,
            ),
            Quantity(
                'cc_reference_wavefunction',
                r'Reference Wavefunction\s+\.\.\.\s+([A-Z]+)',
                repeats=False,
            ),
            Quantity(
                't1_diagnostic',
                r'T1 diagnostic\s+\.\.\.\s+([\d\.]+)',
                repeats=False,
            ),
            Quantity(
                'largest_t2_amplitudes',
                # r'\b\d+[ab]->\d+[ab]\b\s+\b\d+[ab]->\d+[ab]\b\s+([-+]?\d*\.\d+)\b',
                # r'\b\d+[ab]?->\d+[ab]?\s+\d+[ab]?->\d+[ab]?\s+([-+]?\d*\.\d+)',
                # r'(\d+[ab]?->\d+[ab]?)\s+(\d+[ab]?->\d+[ab]?)\s+([-+]?\d*\.\d+)',
                # r'\b\d+[ab]?->\d+[ab]?\s+\d+[ab]?->\d+[ab]?\s+([-+]?\d*\.\d+)',
                # r'\b(\d+[ab]?->\d+[ab]?)\s+(\d+[ab]?->\d+[ab]?)\s+([-+]?\d*\.\d+)\b',
                r'\b\d+[ab]?\s*->\s*\d+[ab]?\s+\d+[ab]?\s*->\s*\d+[ab]?\s+([-+]?\d*\.\d+)',
                repeats=True,
            ),  # TODO: find regex that is suitable for all 4 types
            Quantity(
                'reference_energy',
                r'E\(0\)\s*\.\.\.\s*(-?[\d\.]+)',
                repeats=False,
            ),
            Quantity(
                'corr_energy_strong',
                r'E\(CORR\)\(strong-pairs\)\s*\.\.\.\s*(-?[\d\.]+)',
                repeats=False,
            ),
            Quantity(
                'corr_energy_weak',
                r'E\(CORR\)\(weak-pairs\)\s*\.\.\.\s*(-?[\d\.]+)',
                repeats=False,
            ),
            Quantity(
                'tCutPairs',
                r'TCutPairs\s*=\s*([-+]?\d*\.\d+(?:[eE][-+]?\d+)?)',
                repeats=False,
            ),
            Quantity(
                'tCutPNO',
                r'TCutPNO\s*=\s*([-+]?\d*\.\d+(?:[eE][-+]?\d+)?)',
                repeats=False,
            ),
            Quantity(
                'tCutPNOSingles',
                r'TCutPNOSingles\s*=\s*([-+]?\d*\.\d+(?:[eE][-+]?\d+)?)',
                repeats=False,
            ),
            Quantity(
                'tCutMP2Pairs',
                r'TCutMP2Pairs\s*[:=]\s*([-+]?\d*\.\d+(?:[eE][-+]?\d+)?)',
                repeats=False,
            ),
            Quantity(
                'perturbative_triple_excitations_on_off',
                r'Perturbative triple excitations\s*\.+\s*(\S+)',
                convert=False,
            ),
            Quantity(
                'f12_correction_on_off',
                r'Calculation of F12 correction\s*\.+\s*(\S+)',
                convert=False,
            ),
            Quantity(
                'kc_formation',
                r'K\(C\) Formation\s*\.+\s*(\S+)',
                convert=False,
            ),
        ]

        # Basis set related quantities
        basis_set_quantities = [
            Quantity('basis_set_atom_labels', r'Type\s*(\w+)', repeats=True),
            Quantity('basis_set', r':\s*(\w+)\s*contracted\s*to', repeats=True),
            Quantity('basis_set_contracted', r'(\w+)\s*pattern', repeats=True),
        ]

        basis_set_naming_quantities = [
            Quantity(
                'main_basis_set',
                # r'Your calculation utilizes the basis:\s*(.*)',
                r'Your calculation utilizes the basis:\s*([^\s].*?)[\s]*\n',
                repeats=False,
                convert=False,
                flatten=False,
            ),
            Quantity(
                'auxj_basis_set',
                (
                    r'----- AuxJ basis set information -----\s*'
                    r'Your calculation utilizes the auxiliary basis:\s*'
                    r'([^\s].*?)(?:\s*\n|$)'
                ),
                repeats=False,
                convert=False,
            ),
            Quantity(
                'auxjk_basis_set',
                (
                    r'----- AuxJK basis set information -----\s*'
                    r'Your calculation utilizes the auxiliary basis:\s*'
                    r'([^\s].*?)(?:\s*\n|$)'
                ),
                repeats=False,
                convert=False,
            ),
            Quantity(
                'auxc_basis_set',
                (
                    r'----- AuxC basis set information -----\s*'
                    r'Your calculation utilizes the auxiliary basis:\s*'
                    r'([^\s].*?)(?:\s*\n|$)'
                ),
                repeats=False,
                convert=False,
                flatten=False,
            ),
            Quantity(
                'capped_ecp',
                # r'Group\s*\d+,\s*Type\s*([A-Z][a-z]*)\s*ECP\s*([A-Z]+\(\d+,[A-Z]+\))',
                r'Group \d+, Type (\w+)\s+ECP\s+(\w+\(\d+,[A-Z]+\))',
                repeats=True,
                convert=False,
            ),
        ]

        # Basis set statistics quantities
        basis_set_statistics_quantities = [
            Quantity(
                'nb_of_primitive_gaussian_shells',
                r'# of primitive gaussian shells\s*\.+\s*(\d+)',
                repeats=True,
                dtype=int,
            ),
            Quantity(
                'nb_of_primitive_gaussian_functions',
                r'# of primitive gaussian functions\s*\.+\s*(\d+)',
                repeats=True,
                dtype=int,
            ),
            Quantity(
                'nb_of_contracted_shells',
                r'# of contracted shells\s*\.+\s*(\d+)',
                repeats=True,
                dtype=int,
            ),
            Quantity(
                'nb_of_contracted_basis_functions',
                r'# of contracted (?:aux-)?basis functions\s*\.+\s*(\d+)',
                repeats=True,
                dtype=int,
            ),
            Quantity(
                'highest_angular_moment',
                r'Highest angular momentum\s*\.+\s*(\d+)',
                repeats=True,
                dtype=int,
            ),
            Quantity(
                'maximum_contraction_depth',
                r'Maximum contraction depth\s*\.+\s*(\d+)',
                repeats=True,
                dtype=int,
            ),
            Quantity(
                'main_basis_set',
                r'Number of basis functions\s*\.+\s*(\d+)',
                repeats=False,
                dtype=int,
            ),
            Quantity(
                'auxj_basis_set',
                r'   # of basis functions in Aux-J\s*\.+\s*(\d+)',
                repeats=False,
                dtype=int,
            ),
            Quantity(
                'auxjk_basis_set',
                r'   # of basis functions in Aux-JK\s*\.+\s*(\d+)',
                repeats=False,
                dtype=int,
            ),
            Quantity(
                'auxc_basis_set',
                r'   # of basis functions in Aux-C\s*\.+\s*(\d+)',
                repeats=False,
                dtype=int,
            ),
        ]

        # Grid related quantities
        grid_quantities = [
            Quantity(
                'gral_integ_accuracy',
                rf'General Integration Accuracy\s*IntAcc\s*\.+\s*({re_float})',
                dtype=float,
            ),
            Quantity(
                'radial_grid_type',
                r'Radial Grid Type\s*RadialGrid\s*\.+\s*(\S+)',
                convert=False,
            ),
            Quantity(
                'angular_grid',
                r'Angular Grid \(max\. acc\.\)\s*AngularGrid\s*\.+\s*(\S+)',
                convert=False,
            ),
            Quantity(
                'grid_pruning_method',
                r'Angular grid pruning method\s*GridPruning\s*\.+\s*(.+)',
                flatten=False,
                convert=False,
            ),
            Quantity(
                'weight_gener_scheme',
                r'Weight generation scheme\s*WeightScheme\s*\.+\s*(\w+)',
                convert=False,
            ),
            Quantity(
                'basis_fn_cutoff',
                rf'Basis function cutoff\s*BFCut\s*\.+\s*({re_float})',
                dtype=float,
            ),
            Quantity(
                'integr_weight_cutoff',
                rf'Integration weight cutoff\s*WCut\s*\.+\s*({re_float})',
                dtype=float,
            ),
            Quantity(
                'nb_grid_pts_after_initial_pruning',
                r'# of grid points \(after initial pruning\)\s*\.+\s*(\d+)',
                dtype=int,
            ),
            Quantity(
                'nb_grid_pts_after_weights_screening',
                r'# of grid points \(after weights\+screening\)\s*\.+\s*(\d+)',
                dtype=int,
            ),
            Quantity(
                'total_nb_grid_pts',
                r'Total number of grid points\s*\.+\s*(\d+)',
                dtype=int,
            ),
            Quantity(
                'total_nb_batches', r'Total number of batches\s*\.+\s*(\d+)', dtype=int
            ),
            Quantity(
                'avg_nb_points_per_batch',
                r'Average number of points per batch\s*\.+\s*(\d+)',
                dtype=int,
            ),
            Quantity(
                'avg_nb_grid_pts_per_atom',
                r'Average number of grid points per atom\s*\.+\s*(\d+)',
                dtype=int,
            ),
        ]

        # SCF convergence quantities
        scf_convergence_quantities = [
            Quantity(
                name.lower().replace(' ', '_').replace('-', '_'),
                rf'%s\s*\.+\s*({re_float})\s* Tolerance :\s*({re_float})' % name,
                dtype=float,
            )
            for name in [
                'Last Energy change',
                'Last MAX-Density change',
                'Last RMS-Density change',
            ]
        ]

        # Population analysis quantities
        population_quantities = [
            Quantity(
                'atomic_charges',
                r'[A-Z]+ ATOMIC CHARGES.*\n\-+([\s\S]+?)\-{10}',
                sub_parser=TextParser(
                    quantities=[
                        Quantity('species', r'\n *\d+\s*(\w+)', repeats=True),
                        Quantity(
                            'charge', rf':\s*({re_float})', repeats=True, dtype=float
                        ),
                        Quantity(
                            'total_charge',
                            rf'Sum of atomic charges\s*:\s*({re_float})',
                            dtype=float,
                        ),
                    ]
                ),
            ),
            Quantity(
                'orbital_charges',
                rf'[A-Z]+ REDUCED ORBITAL CHARGES.*\s*\-+([\s\S]+?{re_n}{re_n})',
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'atom',
                            r'([A-Z][a-z]?\s*[spdf][\s\S]+?)\n *(?:\d|\Z)',
                            repeats=True,
                            sub_parser=TextParser(
                                quantities=[
                                    Quantity(
                                        'species', r'([A-Z][a-z]?)', convert=False
                                    ),
                                    Quantity(
                                        'charge',
                                        rf'([spdf]\S*)\s*:\s*({re_float})',
                                        repeats=True,
                                    ),
                                ]
                            ),
                        )
                    ]
                ),
            ),
        ]

        # scf quantities
        self_consistent_quantities = [
            Quantity(
                'scf_settings',
                r'SCF SETTINGS\s*\-+([\s\S]+?)\-{12}',
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'exchange_functional',
                            r'Exchange Functional\s+Exchange\s*\.+\s+(\w+)',
                            convert=False,
                        ),
                        Quantity(
                            'correlation_functional',
                            r'Correlation Functional\s*Correlation\s*\.+\s+(\w+)',
                            convert=False,
                        ),
                        Quantity(
                            'n_max_iterations',
                            rf'MaxIter\s+\.*\s*({re_float})',
                            dtype=float,
                        ),
                        Quantity(
                            'XC_functional_type',
                            r'Ab initio Hamiltonian\s*Method\s*\.+\s*(\S+)',
                            convert=False,
                        ),
                        Quantity(
                            'XC_functional_type',
                            r'Density Functional\s*Method\s*\.+\s*(\S+)',
                            convert=False,
                        ),
                        Quantity(
                            'xalpha_param',
                            rf'X-Alpha parameter\s*XAlpha\s*\.+\s*({re_float})',
                            dtype=float,
                        ),
                        Quantity(
                            'beckes_beta_param',
                            rf'Becke\'s b parameter\s*XBeta\s*\.+\s*({re_float})',
                            dtype=float,
                        ),
                        Quantity(
                            'correl_functional',
                            r'Correlation Functional Correlation\s*\.+\s*(\S+)',
                            convert=False,
                        ),
                        Quantity(
                            'lda_part_of_gga_corr',
                            r'LDA part of GGA corr\.\s*LDAOpt\s*\.+\s*(\S+)',
                            convert=False,
                        ),
                        Quantity(
                            'fraction_hf_exchange',
                            r'   Fraction HF Exchange\s*ScalHFX\s*\.+\s*([\d\.]+)',
                            dtype=float,
                        ),
                        Quantity(
                            'scaling_exchange',
                            r'Scaling of DF-GGA-X\s*ScalDFX\s*\.+\s*([\d\.]+)',
                            dtype=float,
                        ),
                        Quantity(
                            'scaling_correlation',
                            r'Scaling of DF-GGA-C\s*ScalDFC\s*\.+\s*([\d\.]+)',
                            dtype=float,
                        ),
                        Quantity(
                            'scalar_relativistic_method',
                            r'Scalar relativistic method\s*\.+\s*(\w+)',
                            convert=False,
                        ),
                        Quantity(
                            'speed_of_light_used',
                            rf'Speed of light used\s*Velit\s*\.+\s*({re_float})',
                            dtype=float,
                        ),
                        Quantity(
                            'hf_type',
                            r'Hartree-Fock type\s*HFTyp\s*\.+\s*(\w+)',
                            convert=False,
                        ),
                        Quantity(
                            'total_charge',
                            rf'Total Charge\s*Charge\s*\.+\s*({re_float})',
                            dtype=float,
                        ),
                        Quantity(
                            'multiplicity',
                            rf'Multiplicity\s*Mult\s*\.+\s*({re_float})',
                            dtype=float,
                        ),
                        Quantity(
                            'nelectrons',
                            rf'Number of Electrons\s*NEL\s*\.+\s*({re_float})',
                            dtype=float,
                        ),
                        Quantity(
                            'basis_dimension',
                            rf'Basis Dimension\s*Dim\s*\.+\s*({re_float})',
                            dtype=float,
                        ),
                        Quantity(
                            'nuclear_repulsion',
                            rf'Nuclear Repulsion\s*ENuc\s*\.+\s*({re_float})',
                            dtype=float,
                            unit=ureg.hartree,
                        ),
                        Quantity(
                            'convergence_check_mode',
                            r'Convergence Check Mode ConvCheckMode\s*\.+\s*(\S+)',
                            convert=False,
                        ),
                        Quantity(
                            'energy_change_tolerance',
                            # rf'Energy Change\s*TolE\s*\.+\s*({re_float})',
                            rf'Energy\s+Change\s+TolE\s+\.*\s*({re_float})\s*(?:Eh)?',
                            dtype=float,
                            unit=ureg.hartree,
                        ),
                        Quantity(
                            '1_elect_energy_change',
                            rf'1\-El\. energy change\s*\.+\s*({re_float})',
                            dtype=float,
                        ),
                        Quantity(
                            'rij',
                            r'RI-approximation to the Coulomb term is turned (\w+)',
                            convert=False,
                        ),
                        Quantity(
                            'cosx',
                            (
                                r'RIJ-COSX \(HFX calculated with COS-X\)\)\s*'
                                r'\.\.\.\.\s*(\w+)'
                            ),
                            convert=False,
                        ),
                        Quantity(
                            'rijk',
                            r'RI-JK \(J\+K treated both via RI\)\s+\.\.\.\.\s+(\w+)',
                            convert=False,
                        ),
                    ]
                ),
            ),
            Quantity(
                'dft_grid_generation',
                r'DFT GRID GENERATION\s*\-+([\s\S]+?\-{10})',
                sub_parser=TextParser(quantities=grid_quantities),
            ),
            Quantity(
                'scf_iterations',
                r'SCF ITERATIONS\s*\-+([\s\S]+?)\*{10}',
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'energy',
                            rf'\n *\d+\s*({re_float})\s*{re_float}',
                            repeats=True,
                            dtype=float,
                            unit=ureg.hartree,
                        )
                    ]
                ),
            ),
            Quantity(
                'final_grid',
                r'Setting up the final grid:([\s\S]+?)\-{10}',
                sub_parser=TextParser(quantities=grid_quantities),
            ),
            Quantity(
                'total_scf_energy',
                r'TOTAL SCF ENERGY\s*\-+([\s\S]+?)\-{10}',
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            name,
                            rf'%s\s*:\s*({re_float})' % key,
                            dtype=float,
                            unit=ureg.hartree,
                        )
                        for key, name in self._energy_mapping.items()
                    ]
                    + [
                        Quantity(
                            'virial_ratio',
                            rf'Virial Ratio\s*:\s*({re_float})',
                            dtype=float,
                        ),
                        Quantity(
                            'nb_elect_alpha_channel',
                            rf'N\(Alpha\)\s*:\s*({re_float})',
                            dtype=float,
                        ),
                        Quantity(
                            'nb_elect_beta_channel',
                            rf'N\(Beta\)\s*:\s*({re_float})',
                            dtype=float,
                        ),
                        Quantity(
                            'nb_elect_total',
                            rf'N\(Total\)\s*:\s*({re_float})',
                            dtype=float,
                        ),
                    ]
                ),
            ),
            Quantity(
                'scf_convergence',
                r'SCF CONVERGENCE\s*\-+([\s\S]+?)\-{10}',
                sub_parser=TextParser(quantities=scf_convergence_quantities),
            ),
            Quantity(
                'orbital_energies',
                rf'NO\s*OCC\s*E\(Eh\)\s*E\(eV\)\s*([\s\S]+?){re_n}{re_n}',
                str_operation=lambda x: np.array(
                    [v.split()[:4] for v in x.split('\n')], dtype=float
                ),
                repeats=True,
            ),
            Quantity(
                'mulliken',
                r'MULLIKEN POPULATION ANALYSIS \*\s*\*+([\s\S]+?)\*{10}',
                sub_parser=TextParser(quantities=population_quantities),
            ),
            Quantity(
                'timings',
                r'\n *TIMINGS\s*\-+\s*([\s\S]+?)\-{10}',
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            name,
                            rf'%s\s*\.+\s*({re_float})' % key,
                            dtype=float,
                            unit=ureg.s,
                        )
                        for key, name in self._timing_mapping.items()
                    ]
                ),
            ),
            Quantity(
                'time_calculation',
                r'Total SCF time\: (\d+) days (\d+) hours (\d+) min (\d+) sec ',
                dtype=np.dtype(np.float64),
            ),
        ]

        # TODO parse more properties, add to metainfo
        tddft_quantities = [
            Quantity(
                'absorption_spectrum_electric',
                r'ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS\s*'
                r'\-+[\s\S]+?\-+\n([\s\S]+?)\-{10}',
                str_operation=lambda x: [v.split() for v in x.strip().split('\n')],
            )
        ]

        def str_to_iteration_energy(val_in):
            val = [v.split() for v in val_in.strip().split('\n')]
            keys = val[0]
            val = np.transpose(
                np.array([v for v in val[1:] if len(v) == len(keys)], dtype=float)
            )
            return {keys[i]: val[i] for i in range(len(keys))}

        ci_quantities = [
            Quantity(
                'electronic_structure_method',
                r'Correlation treatment\s*\.+\s*(\S+)',
                convert=False,
            ),
            Quantity(
                'single_excitations_on_off',
                r'Single excitations\s*\.+\s*(\S+)',
                convert=False,
            ),
            Quantity(
                'orbital_opt_on_off',
                r'Orbital optimization\s*\.+\s*(\S+)',
                convert=False,
            ),
            Quantity(
                'z_vector_calc_on_off',
                r'Calculation of Z vector\s*\.+\s*(\S+)',
                convert=False,
            ),
            Quantity(
                'Brueckner_orbitals_calc_on_off',
                r'Calculation of Brueckner orbitals\s*\.+\s*(\S+)',
                convert=False,
            ),
            Quantity(
                'perturbative_triple_excitations_on_off',
                r'Perturbative triple excitations\s*\.+\s*(\S+)',
                convert=False,
            ),
            Quantity(
                'f12_correction_on_off',
                r'Calculation of F12 correction\s*\.+\s*(\S+)',
                convert=False,
            ),
            Quantity(
                'frozen_core_treatment',
                r'Frozen core treatment\s*\.+\s*(.+)',
                flatten=False,
                convert=False,
            ),
            Quantity(
                'reference_wave_function',
                r'Reference Wavefunction\s*\.+\s*(.+)',
                flatten=False,
                convert=False,
            ),
            Quantity(
                'nb_of_atomic_orbitals', r'Number of AO\'s\s*\.+\s*(\d+)', dtype=int
            ),
            Quantity(
                'nb_of_electrons', r'Number of electrons\s*\.+\s*(\d+)', dtype=int
            ),
            Quantity(
                'nb_of_correlated_electrons',
                r'Number of correlated electrons\s*\.+\s*(\d+)',
                dtype=int,
            ),
            Quantity(
                'integral_transformation',
                r'Integral transformation\s*\.+\s*(.+)',
                flatten=False,
                convert=False,
            ),
            Quantity(
                'level_shift_amplitude_update',
                rf'Level shift for amplitude update\s*\.+\s*({re_float})',
                dtype=float,
            ),
            Quantity(
                'coulomb_transformation_type',
                r'Transformation type\s*\.+\s*(.+)',
                flatten=False,
                convert=False,
            ),
            Quantity(
                'coulomb_transformation_dimension_basis',
                r'Dimension of the basis\s*\.+\s*(\d+)',
                dtype=int,
            ),
            Quantity(
                'nb_internal_alpha_mol_orbitals',
                r'Number of internal alpha\-MOs\s*\.+\s*(\d+)',
                dtype=int,
            ),
            Quantity(
                'nb_internal_beta_mol_orbitals',
                r'Number of internal beta\-MOs\s*\.+\s*(\d+)',
                dtype=int,
            ),
            Quantity('pair_cutoff', rf'Pair cutoff\s*\.+\s*({re_float})', dtype=float),
            Quantity(
                'atomic_orbital_integral_source',
                r'AO\-integral source\s*\.+\s*(.+)',
                flatten=False,
                convert=False,
            ),
            Quantity(
                'integral_package_used',
                r'Integral package used\s*\.+\s*(.+)',
                flatten=False,
                convert=False,
            ),
            Quantity(
                'nb_alpha_pairs_included',
                r'Number of Alpha\-MO pairs included\s*\.+\s*(\d+)',
                dtype=int,
            ),
            Quantity(
                'nb_beta_pairs_included',
                r'Number of Beta\-MO pairs included\s*\.+\s*(\d+)',
                dtype=int,
            ),
            Quantity(
                'mp2_energy_spin_aa',
                rf'EMP2\(aa\)=\s*({re_float})',
                dtype=float,
                unit=ureg.hartree,
            ),
            Quantity(
                'mp2_energy_spin_bb',
                rf'EMP2\(bb\)=\s*({re_float})',
                dtype=float,
                unit=ureg.hartree,
            ),
            Quantity(
                'mp2_energy_spin_ab',
                rf'EMP2\(ab\)=\s*({re_float})',
                dtype=float,
                unit=ureg.hartree,
            ),
            Quantity(
                'mp2_initial_guess',
                rf'E\(0\)\s*\.+\s*({re_float})',
                dtype=float,
                unit=ureg.hartree,
            ),
            Quantity(
                'mp2_energy',
                rf'E\(MP2\)\s*\.+\s*({re_float})',
                dtype=float,
                unit=ureg.hartree,
            ),
            Quantity(
                'mp2_total_energy',
                rf'Initial E\(tot\)\s*\.+\s*({re_float})',
                dtype=float,
                unit=ureg.hartree,
            ),
            Quantity(
                'sl_mp2_correlation_energy',
                rf'E\(SL\-MP2\)\s*\.+\s*({re_float})',
                dtype=float,
                unit=ureg.hartree,
            ),
            Quantity(
                'T_and_T_energy',
                rf'<T\|T>\s*\.+\s*({re_float})',
                dtype=float,
                unit=ureg.hartree,
            ),
            Quantity(
                'local_mp2_treatment',
                r'Type of local MP2 treatment\s*\.+\s*(.+)',
                convert=False,
                flatten=False,
            ),
            Quantity(
                'pair_density_normalization',
                r'Pair density normalization\s*\.+\s*(.+)',
                convert=False,
                flatten=False,
            ),
            Quantity(
                'spin_component_scaling',
                r'Spin component scaling\s*\.+\s*(.+)',
                convert=False,
                flatten=False,
            ),
            Quantity(
                'tCutMP2Pairs',
                r'TCutMP2Pairs\s*[:=]\s*([-+]?\d*\.\d+(?:[eE][-+]?\d+)?)',
                dtype=float,
            ),
            Quantity(
                'total_nb_pairs_included',
                r'Number of pairs included\s*\.+\s*(\d+)',
                dtype=int,
            ),
            Quantity(
                'iteration_energy',
                r'(Iter\s*E\(tot\)[\s\S]+?)\-{3}',
                str_operation=str_to_iteration_energy,
                convert=False,
            ),
            Quantity(
                'ccsd_correlation_energy',
                rf'E\(CORR\)\s*\.+\s*({re_float})',
                dtype=float,
                unit=ureg.hartree,
            ),
            Quantity(
                'ccsd_total_energy',
                rf'E\(TOT\)\s*\.+\s*({re_float})',
                dtype=float,
                unit=ureg.hartree,
            ),
            Quantity(
                'single_norm_half_ss',
                rf'Singles Norm <S\|S>\*\*1/2\s*\.+\s*({re_float})',
                dtype=float,
                unit=ureg.hartree,
            ),
            Quantity(
                't1_diagnostic',
                rf'T1 diagnostic\s*\.+\s*({re_float})',
                dtype=float,
                unit=ureg.hartree,
            ),
            Quantity(
                'ccsdt_total_triples_correction',
                rf'Triples Correction \(T\)\s*\.+\s*({re_float})',
                dtype=float,
                unit=ureg.hartree,
            ),
            Quantity(
                'ccsdt_aaa_triples_contribution',
                rf'alpha\-alpha\-alpha\s*\.+\s*({re_float})',
                dtype=float,
                unit=ureg.hartree,
            ),
            Quantity(
                'ccsdt_aab_triples_contribution',
                rf'alpha\-alpha\-beta\s*\.+\s*({re_float})',
                dtype=float,
                unit=ureg.hartree,
            ),
            # typo in metainfo?
            Quantity(
                'ccsdt_aba_triples_contribution',
                rf'alpha\-beta\-beta\s*\.+\s*({re_float})',
                dtype=float,
                unit=ureg.hartree,
            ),
            Quantity(
                'ccsdt_bbb_triples_contribution',
                rf'beta\-beta\-beta\s*\.+\s*({re_float})',
                dtype=float,
                unit=ureg.hartree,
            ),
            Quantity(
                'ccsdt_final_corr_energy',
                rf'Final correlation energy\s*\.+\s*({re_float})',
                dtype=float,
                unit=ureg.hartree,
            ),
            Quantity(
                'ccsd_final_energy',
                rf'E\(CCSD\)\s*\.+\s*({re_float})',
                dtype=float,
                unit=ureg.hartree,
            ),
            Quantity(
                'energy_total',
                rf'E\(CCSD\(T\)\)\s*\.+\s*({re_float})',
                dtype=float,
                unit=ureg.hartree,
            ),
        ]

        mp2_quantities = [
            Quantity(
                'mp2_basis_dimension',
                r'Dimension of the basis\s*\.+\s*(\d+)',
                dtype=int,
            ),
            Quantity(
                'scaling_mp2_energy',
                rf'Overall scaling of the MP2 energy\s*\.+\s*({re_float})',
                dtype=float,
            ),
            Quantity(
                'mp2_aux_basis_dimension',
                r'Dimension of the aux\-basis\s*\.+\s*(\d+)',
                dtype=int,
            ),
            Quantity(
                'energy_method_current',
                rf'RI\-MP2 CORRELATION ENERGY:\s*({re_float})',
                dtype=float,
                unit=ureg.hartree,
            ),
            Quantity(
                'energy_total',
                rf'MP2 TOTAL ENERGY:\s*({re_float})',
                dtype=float,
                unit=ureg.hartree,
            ),
        ]

        localization_quantities = [
            Quantity(
                'type',
                r'Localization criterion\s*\.+\s*(\S+)',
                convert=False,
            ),
            Quantity(
                'n_max_iterations',
                rf'Max. number of iterations\s*\.+\s*({re_float})',
                dtype=float,
            ),
            Quantity(
                'energy_change_tolerance',
                rf'Convergence tolerance\s*\.+\s*({re_float})',
                dtype=float,
                unit=ureg.hartree,
            ),
            Quantity(
                'orbital_range',
                # r'Orbital range for localization\s*\.+\s*(\S+)',
                # r'Orbital range for localization\s*\.+\s*(\d+\s+to\s+\d+)',
                r'Orbital range for localization\s*\.+\s*(\d+)\s+to\s+(\d+)',
                convert=False,
            ),
        ]

        calculation_quantities = [
            Quantity(
                'cartesian_coordinates',
                rf'CARTESIAN COORDINATES \(ANGSTROEM\)\s*\-+\s*([\s\S]+?){re_n}{re_n}',
                # str_operation=str_to_cartesian_coordinates,
                repeats=False,
            ),
            Quantity(
                'basis_set',
                r'\n *BASIS SET INFORMATION\s*\-+([\s\S]+?)\-{10}',
                sub_parser=TextParser(quantities=basis_set_quantities),
            ),
            Quantity(
                'auxiliary_basis_set',
                r'\n *AUXILIARY BASIS SET INFORMATION\s*\-+([\s\S]+?)\-{10}',
                sub_parser=TextParser(quantities=basis_set_quantities),
            ),
            Quantity(
                'basis_set_statistics',
                r'BASIS SET STATISTICS AND STARTUP INFO([\s\S]+?)\-{10}',
                sub_parser=TextParser(quantities=basis_set_statistics_quantities),
            ),
            Quantity(
                'self_consistent',
                r'((?:ORCA SCF|DFT GRID GENERATION)\s*\-+[\s\S]+?(?:\-{70}|\Z))',
                sub_parser=TextParser(quantities=self_consistent_quantities),
            ),
            Quantity(
                'tddft',
                r'ORCA TD\-DFT(?:/TDA)* CALCULATION\s*\-+\s*([\s\S]+?E\(tot\).*)',
                sub_parser=TextParser(quantities=tddft_quantities),
            ),
            Quantity(
                'mp2',
                (
                    r'(ORCA MP2 CALCULATION|-{78}\s+ORCA\s+MP2\s+-{78})'
                    r'([\s\S]+?MP2 TOTAL ENERGY:.+)'
                ),
                sub_parser=TextParser(quantities=mp2_quantities),
            ),
            Quantity(
                'ci',
                r'ORCA\-MATRIX DRIVEN CI([\s\S]+?E\(CCSD\(T\)\).*)',
                sub_parser=TextParser(quantities=ci_quantities),
            ),
            # CASSCF / multireference block
            Quantity(
                'casscf',
                (
                    r'ORCA\-CASSCF\s*\-+\s*'
                    r'([\s\S]+?CAS\-?SCF ITERATIONS)'
                ),
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'n_active_electrons',
                            r'Number of active electrons\s*\.+\s*(\d+)',
                            dtype=int,
                        ),
                        Quantity(
                            'n_active_orbitals',
                            r'Number of active orbitals\s*\.+\s*(\d+)',
                            dtype=int,
                        ),
                        Quantity(
                            'n_total_electrons',
                            r'Total number of electrons\s*\.+\s*(\d+)',
                            dtype=int,
                        ),
                        Quantity(
                            'n_total_orbitals',
                            r'Total number of orbitals\s*\.+\s*(\d+)',
                            dtype=int,
                        ),
                        Quantity(
                            'n_blocks',
                            r'Number of multiplicity blocks\s*\.+\s*(\d+)',
                            dtype=int,
                        ),
                        Quantity(
                            'block',
                            r'(BLOCK\s+\d+\s+WEIGHT=\s*[-+]?\d*\.?\d+\s*[\s\S]+?)(?=(?:BLOCK|\Z))',
                            repeats=True,
                            sub_parser=TextParser(
                                quantities=[
                                    Quantity(
                                        'weight',
                                        r'BLOCK\s+\d+\s+WEIGHT=\s*([-+]?\d*\.?\d+)',
                                        dtype=float,
                                    ),
                                    Quantity(
                                        'multiplicity',
                                        r'Multiplicity\s*\.+\s*(\d+)',
                                        dtype=int,
                                    ),
                                    Quantity(
                                        'n_configurations',
                                        r'#\(Configurations\)\s*\.+\s*(\d+)',
                                        dtype=int,
                                    ),
                                    Quantity(
                                        'n_csfs',
                                        r'#\(CSFs\)\s*\.+\s*(\d+)',
                                        dtype=int,
                                    ),
                                    Quantity(
                                        'n_roots',
                                        r'#\(Roots\)\s*\.+\s*(\d+)',
                                        dtype=int,
                                    ),
                                    Quantity(
                                        'root_weights',
                                        r'ROOT=\d+\s+WEIGHT=\s*([-+]?\d*\.?\d+)',
                                        repeats=True,
                                        dtype=float,
                                    ),
                                    Quantity(
                                        'casci_marker',
                                        r'(MaxMacroIter\s+1\s+detected\s+>>>\s*CAS\-CI)',
                                        convert=False,
                                    ),
                                ]
                            ),
                        ),
                    ]
                ),
            ),
            Quantity(
                'loc',
                r'\n *ORCA ORBITAL LOCALIZATION\s*\-+([\s\S]+?)\-{10}',
                repeats=True,
                sub_parser=TextParser(quantities=localization_quantities),
            ),
            # FIX HERE
            Quantity(
                'cc',
                r'ORCA\-MATRIX DRIVEN CI([\s\S]+?E\(CCSD\(T\)\).*)',
                # optional regex: ECCSDT part is not a must to match
                # r'ORCA\-MATRIX DRIVEN CI([\s\S]+?)(E\(CCSD\(T\)\).*)?',
                sub_parser=TextParser(quantities=coupled_cluster_quantities),
            ),
        ]

        geometry_optimization_quantities = [
            Quantity(
                f'{key.lower().replace(" ", "_").replace(".", "")}_tol',
                rf'{key}\s*(\w+)\s*\.+\s*({re_float})',
                dtype=float,
            )
            for key in [
                'Energy Change',
                'Max. Gradient',
                'RMS Gradient',
                'Max. Displacement',
                'RMS Displacement',
            ]
        ]

        geometry_optimization_quantities += [
            Quantity('update_method', r'Update method\s*(\w+)\s*\.+\s*(.+)'),
            Quantity('coords_choice', r'Choice of coordinates\s*(\w+)\s*\.+\s*(.+)'),
            Quantity('initial_hessian', r'Initial Hessian\s*(\w+)\s*\.+\s*(.+)'),
        ]

        geometry_optimization_quantities += [
            Quantity(
                'cycle',
                (
                    r'OPTIMIZATION CYCLE\s*\d+\s*\*\s*\*+'
                    r'([\s\S]+?)(?:\*\s*GEOMETRY|OPTIMIZATION RUN DONE|\Z)'
                ),
                repeats=True,
                sub_parser=TextParser(quantities=calculation_quantities),
            ),
            Quantity(
                'final_energy_evaluation',
                (
                    r'FINAL ENERGY EVALUATION AT THE STATIONARY POINT'
                    r'([\s\S]+?FINAL SINGLE POINT ENERGY.*)'
                ),
                sub_parser=TextParser(quantities=calculation_quantities),
            ),
        ]

        self._quantities = [
            Quantity(
                'program_version',
                # r'Program Version\s*([\w_.].*)',
                r'Program Version\s*([\d\.]+)',
                convert=False,
                flatten=False,
            ),
            Quantity(
                'program_svn', r'\(SVN:\s*\$([^$]+)\$\)\s', convert=False, flatten=False
            ),
            Quantity(
                'program_compilation_date',
                r'\(\$Date\:\s*(\w.+?)\s*\$\)',
                convert=False,
                flatten=False,
            ),
            Quantity(
                'input_file',
                # r'INPUT FILE\s*\=+([\s\S]+?)END OF INPUT',
                r'INPUT FILE\s*=\s*([\s\S]+?)(?:(?:^=+$|END OF INPUT))',
                # sub_parser=TextParser(
                #    quantities=[
                #        Quantity('xc_functional', r'\d+>\s*!\s*(\S+)'),
                #        Quantity('tier', r'(\w+SCF)'),
                #    ]
                # ),
            ),
            Quantity(
                'basis_set_name',
                # r'----- Orbital basis set information -----([\s\S]+?)\={10}',
                r'----- Orbital basis set information -----([\s\S]+?)\={10}\s*',
                sub_parser=TextParser(quantities=basis_set_naming_quantities),
            ),
            Quantity(
                'basis_set_total',
                # r'----- Orbital basis set information -----([\s\S]+?)\={10}',
                r'ORCA GTO INTEGRAL CALCULATION([\s\S]+?)ORCA SCF',
                sub_parser=TextParser(quantities=basis_set_statistics_quantities),
            ),
            Quantity(
                'ecp_basis_set_name',
                r'ECP PARAMETER INFORMATION\s*\-+\s*([\s\S]+?)Atom\s*',
                sub_parser=TextParser(quantities=basis_set_naming_quantities),
            ),
            Quantity(
                'single_point',
                (
                    r'\* Single Point Calculation \*\s*\*+'
                    r'([\s\S]+?(?:FINAL SINGLE POINT ENERGY.*|\Z))'
                ),
                sub_parser=TextParser(quantities=calculation_quantities),
            ),
            Quantity(
                'geometry_optimization',
                (
                    r'\* Geometry Optimization Run \*\s*\*+'
                    r'([\s\S]+?(?:OPTIMIZATION RUN DONE|\Z))'
                ),
                sub_parser=TextParser(quantities=geometry_optimization_quantities),
            ),
        ]
