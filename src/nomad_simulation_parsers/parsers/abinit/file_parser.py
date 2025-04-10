import re

import numpy as np
from ase.data import chemical_symbols
from nomad.parsing.file_parser import Quantity, TextParser
from nomad.units import ureg

RE_FLOAT = r'[-+]?\d+\.\d*(?:[Ee][-+]\d+)?'
RE_N = r'[\n\r]'


class AbinitOutParser(TextParser):
    def __init__(self):
        self.energy_components = {
            'energy_kinetic_electronic': 'Kinetic energy',
            'energy_electronstatic': 'Hartree energy',
            'energy_XC': 'XC energy',
            'ewald': 'Ewald energy',
            'psp_core': 'PspCore energy',
            'psp_local': 'Loc. psp. energy',
            'psp_nonlocal': 'NL   psp  energy',
            'internal': '>>>>> Internal E',
            'energy_correction_entropy': r'\-kT*entropy',
            'energy_total': ' >>>>>>>>> Etotal',
            'energy_sum_eigenvalues': r'Band energy \(Ha\)',
        }
        self._input_vars = None
        self._dataset_numbers = None
        self._n_datasets = None
        super().__init__(None)

    def init_quantities(self):
        self._quantities = [
            Quantity(
                'program_version',
                r'\.Version ([\w\.]+) of ABINIT',
                repeats=False,
                convert=False,
                flatten=False,
            ),
            Quantity(
                'x_abinit_parallel_compilation',
                r'\.\((\w+) version,',
                repeats=False,
                convert=False,
                flatten=False,
            ),
            Quantity(
                'program_compilation_host',
                r'prepared for a ([\w\.]+) computer',
                repeats=False,
                convert=False,
                flatten=False,
            ),
            Quantity(
                'x_abinit_start_date',
                r'\.Starting date : ([\w ]+)\.',
                repeats=False,
                convert=False,
                flatten=False,
            ),
            Quantity(
                'x_abinit_start_time',
                r'\- \( at (\w+) \)',
                repeats=False,
                convert=False,
                flatten=False,
            ),
            Quantity(
                'x_abinit_input_file',
                r'\- input\s*file\s*\-\> ([\w\.]+)',
                repeats=False,
                convert=False,
                flatten=False,
            ),
            Quantity(
                'x_abinit_output_file',
                r'\- output\s*file\s*\-\> ([\w\.]+)',
                repeats=False,
                convert=False,
                flatten=False,
            ),
            Quantity(
                'x_abinit_input_files_root',
                r'\- root for input\s*files \-\> (\w+)',
                repeats=False,
                convert=False,
                flatten=False,
            ),
            Quantity(
                'x_abinit_output_files_root',
                r'\- root for output\s*files \-\> (\w+)',
                repeats=False,
                convert=False,
                flatten=False,
            ),
            Quantity(
                'x_abinit_total_cpu_time',
                r'\-\s*Total cpu\s*time\s*\(s,m,h\):\s*([\d\.]+)',
                dtype=float,
            ),
            Quantity(
                'x_abinit_total_wallclock_time',
                r'\-\s*Total wall clock time\s*\(s,m,h\):\s*([\d\.]+)',
                dtype=float,
            ),
            Quantity(
                'run_clean_end',
                r'(Calculation completed)',
                repeats=False,
                convert=False,
                flatten=False,
            ),
        ]

        self._quantities.append(
            Quantity(
                'input_variables',
                r'\-outvars: echo values of preprocessed input variables '
                r'\-+([\s\S]+?)\={10}',
                repeats=False,
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'key_value',
                            r'([a-zA-Z\d]+)\s*([\d\.\+\-E\s]+)',
                            repeats=True,
                        )
                    ]
                ),
            )
        )

        def str_to_array(val_in):
            val = val_in.strip().split('\n')
            val = [v.split()[-3:] for v in val]
            return np.array(val, dtype=float)

        def str_to_stress_tensor(val_in):
            val = np.array(val_in.split(), dtype=float)
            stress_tensor = np.zeros((3, 3))
            stress_tensor[0][0] = val[0]
            stress_tensor[2][1] = stress_tensor[1][2] = val[1]
            stress_tensor[1][1] = val[2]
            stress_tensor[2][0] = stress_tensor[0][2] = val[3]
            stress_tensor[2][2] = val[4]
            stress_tensor[1][0] = stress_tensor[0][1] = val[5]
            return stress_tensor * (ureg.hartree / ureg.bohr**3)

        def str_to_eigenvalues(val_in):
            return [float(v) for v in val_in.split() if v[-1].isdecimal()]

        self_consistent = [
            Quantity(
                'energy_total_scf_iteration',
                r'ETOT\s*\d+\s*([\+\-\d\.Ee ]+)\n',
                repeats=True,
                dtype=float,
            ),
            Quantity(
                'convergence',
                r'At SCF step\s*([0-9]+)\s*'
                r', etot|, forces|vres2\s*=\s*[\-\+\d\.Ee]+?\s*'
                r'<\s*tolvrs=\s*[\-\+\d\.Ee]+?\s*=>\s*([\w ]+)',
                repeats=False,
            ),
        ]

        relaxation = [
            Quantity(
                'stress_tensor',
                r'Cartesian components of stress tensor \(hartree/bohr\^3\)\s*'
                r'sigma\(\d \d\)\s*=\s*([\-\+\d\.Ee]+)\s*'
                r'sigma\(\d \d\)\s*=\s*([\-\+\d\.Ee]+)\s*'
                r'sigma\(\d \d\)\s*=\s*([\-\+\d\.Ee]+)\s*'
                r'sigma\(\d \d\)\s*=\s*([\-\+\d\.Ee]+)\s*'
                r'sigma\(\d \d\)\s*=\s*([\-\+\d\.Ee]+)\s*'
                r'sigma\(\d \d\)\s*=\s*([\-\+\d\.Ee]+)\s*',
                repeats=False,
                str_operation=str_to_stress_tensor,
                convert=False,
            ),
            Quantity(
                'self_consistent',
                r'SELF\-CONSISTENT\-FIELD CONVERGENCE\-*\s*(iter\s*'
                r'Etot[\s\S]+?)\-*OUTPUT',
                sub_parser=TextParser(quantities=self_consistent),
                repeats=False,
            ),
            Quantity(
                'cartesian_coordinates',
                r'Cartesian coordinates \(xcart\) \[bohr\]\s*([\s\d\.Ee\+\-]+)',
                repeats=False,
                str_operation=str_to_array,
                convert=False,
                unit=ureg.bohr,
            ),
            Quantity(
                'cartesian_forces',
                r'Cartesian forces \(fcart\)[\s\S]+?\(free atoms\)\s*([\s\d\.Ee\+\-]+)',
                repeats=False,
                str_operation=str_to_array,
                convert=False,
                unit=ureg.hartree / ureg.bohr,
            ),
            Quantity(
                'energy_total',
                r'Total energy \(etotal\) \[Ha\]=\s*([\-\+\d\.Ee]+)',
                repeats=False,
                dtype=float,
                unit=ureg.hartree,
            ),
        ]

        energy_components = [
            Quantity(
                key,
                rf'\s*{val}\s*=\s*([\d\.E\-\+]+)',
                repeats=False,
                dtype=float,
                unit=ureg.hartree,
            )
            for key, val in self.energy_components.items()
        ]

        results = [
            Quantity(
                'cartesian_coordinates',
                r'\s*cartesian coordinates \(angstrom\) at end:\s*([\s\d\.Ee\+\-]+)',
                repeats=False,
                str_operation=str_to_array,
                convert=False,
                unit=ureg.angstrom,
            ),
            Quantity(
                'cartesian_forces',
                r'\s*cartesian forces \(hartree/bohr\) at end\:\s*([\s\d\.Ee\+\-]+)',
                repeats=False,
                str_operation=str_to_array,
                convert=False,
                unit=ureg.hartree / ureg.bohr,
            ),
            Quantity(
                'x_abinit_eig_filename',
                r'\s*prteigrs : about to open file\s*(\S+)',
                repeats=False,
                convert=False,
            ),
            Quantity(
                'fermi_energy',
                r'\s*Fermi \(or HOMO\) energy \(hartree\) =\s*([-+0-9.]+)',
                repeats=False,
                dtype=float,
                unit=ureg.hartree,
            ),
            Quantity(
                'x_abinit_magnetisation',
                r'\s*Magnetisation \(Bohr magneton\)=\s*([-+0-9.eEdD]*)\s*',
                repeats=False,
                dtype=float,
                unit=ureg.bohr_magneton,
            ),
            Quantity(
                'eigenvalues',
                r'\s*kpt#\s*(\d+), nband=\s*(\d+), wtk=\s*([\d\.]+)\s*, '
                r'kpt=([\d\.\-\+Ee ]+)\s*\(reduced coord\)\s*([\d\.\-\+Ee\s]+)',
                repeats=True,
                str_operation=str_to_eigenvalues,
                convert=False,
            ),
            Quantity(
                'occupation_numbers',
                r'\s*occupation numbers for kpt#\s*\d+\s*([\d\.\-\+Ee ]+)',
                repeats=True,
                dtype=float,
            ),
            Quantity(
                'energy_total',
                r'Total energy\(eV\)=\s*([\d\.\-\+Ee]+)\s*;',
                repeats=False,
                dtype=float,
                unit=ureg.eV,
            ),
            Quantity(
                'stress_tensor',
                r'Cartesian components of stress tensor \(hartree/bohr\^3\)\s*'
                r'sigma\(\d \d\)\s*=\s*([\-\+\d\.Ee]+)\s*'
                r'sigma\(\d \d\)\s*=\s*([\-\+\d\.Ee]+)\s*'
                r'sigma\(\d \d\)\s*=\s*([\-\+\d\.Ee]+)\s*'
                r'sigma\(\d \d\)\s*=\s*([\-\+\d\.Ee]+)\s*'
                r'sigma\(\d \d\)\s*=\s*([\-\+\d\.Ee]+)\s*'
                r'sigma\(\d \d\)\s*=\s*([\-\+\d\.Ee]+)\s*',
                repeats=False,
                str_operation=str_to_stress_tensor,
                convert=False,
            ),
        ] + energy_components

        self._quantities.append(
            Quantity(
                'dataset',
                r'==\s*(DATASET\s*\d+\s*==[\s\S]+?)(?:\-Cartesian components of '
                r'stress tensor \(GPa\)|\s*prteigrs\s*\:\s*prtvol\=[\w\s\,\-\.\n]*==='
                r'|\s*Average fulfillment[\s\w\[\]\-\:\.\%]*===|== END DATASET)',
                repeats=True,
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'x_abinit_dataset_number',
                            r'DATASET\s*(\d+)',
                            dtype=int,
                            repeats=False,
                        ),
                        Quantity(
                            'x_abinit_var_ixc',
                            r'\-\s*ixc\s*=\s*([\d\-]+)',
                            dtype=int,
                            repeats=False,
                        ),
                        Quantity(
                            'x_abinit_vprim',
                            r'R\(\d\)=\s*(\-*[\d\.]+)\s*(\-*[\d\.]+)\s*'
                            r'(\-*[\d\.]+)\s*G\(\d\)=',
                            repeats=True,
                            dtype=float,
                        ),
                        Quantity(
                            'x_abinit_unit_cell_volume',
                            r'Unit cell volume ucvol\s*=\s*([\d\.\+Ee]+)',
                            dtype=float,
                            unit=ureg.bohr**3,
                        ),
                        Quantity(
                            'self_consistent',
                            r'\={10}\s*(iter\s*Etot[\s\S]+?)\={10}',
                            repeats=False,
                            sub_parser=TextParser(quantities=self_consistent),
                        ),
                        Quantity(
                            'relaxation',
                            r'\-\-\- Iteration: \( \d+\/\d+\)([\s\S]+?Total energy '
                            r'\(etotal\) \[Ha\]=\s*[\-\+\d\.Ee]+)',
                            repeats=True,
                            sub_parser=TextParser(quantities=relaxation),
                        ),
                        Quantity(
                            'results',
                            r'\-+iterations are completed or convergence reached'
                            r'\-+([\s\S]+)'
                            r'(?:(sigma\(\d \d\)\s*=\s*[\+\-\d\.E]+)|={80})',
                            repeats=False,
                            sub_parser=TextParser(quantities=results),
                        ),
                    ]
                ),
            )
        )

        rpa_quantities = [
            Quantity(
                'precision_algorithm',
                rf'{RE_N}\.Using[a-zA-Z\s\;]+\=\s*(\d+)',
                repeats=False,
            ),
            Quantity(
                'kmesh',
                rf'{RE_N}(\s*====\s*K-mesh[a-zA-Z\s]*wavefunctions[\s\S]+?)'
                rf'(?:\s*====\s*Q-mesh)',
                repeats=False,
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'n_mesh',
                            rf'{RE_N}\s*Number of points in the irreducible '
                            rf'wedge\s*\:\s*(\d+)',
                            repeats=False,
                        ),
                        Quantity(
                            'mesh',
                            rf'{RE_N}[\d\)\s]+{RE_FLOAT}\s*{RE_FLOAT}\s*{RE_FLOAT}\s*'
                            rf'{RE_FLOAT}',
                            repeats=True,
                        ),
                    ]
                ),
            ),
            Quantity(
                'qmesh',
                rf'{RE_N}(\s*====\s*Q-mesh[a-zA-Z\s]*screening function[\s\S]+?)'
                rf'(?:\s*setmesh\:)',
                repeats=False,
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'n_mesh',
                            rf'{RE_N}\s*Number of points in the irreducible wedge'
                            rf'\s*\:\s*(\d+)',
                            repeats=False,
                        ),
                        Quantity(
                            'mesh',
                            rf'{RE_N}[\d\)\s]+{RE_FLOAT}\s*{RE_FLOAT}\s*{RE_FLOAT}\s*'
                            rf'{RE_FLOAT}',
                            repeats=True,
                        ),
                    ]
                ),
            ),
            Quantity(
                'fftmesh',
                r'FFT mesh size selected[\s\=]*(\d+)x\s*(\d+)x\s*(\d+)',
                repeats=False,
            ),
            Quantity('n_fftmesh', r'total number of points[\s\=]*(\d+)', repeats=False),
            Quantity(
                'symm_screening',
                rf'{RE_N}\- screening\:([a-zA-Z\-\s]+){RE_N}',
                repeats=False,
            ),
            Quantity(
                'max_band_occ',
                rf'{RE_N}\- Maximum band index for partially occupied states'
                rf'[a-zA-Z\s]+\=\s*(\d+)',
                repeats=False,
            ),
            Quantity(
                'n_bands_per_proc',
                rf'{RE_N}\- Remaining bands to be divided among processors'
                rf'[a-zA-Z\s]+\=\s*(\d+)',
                repeats=False,
            ),
            Quantity(
                'n_bands_per_node',
                rf'{RE_N}\- Number of bands treated by each node[a-zA-Z\s]+\~(\d+)',
                repeats=False,
            ),
            Quantity(
                'n_electrons',
                rf'{RE_N}\s*Number of electrons calculated from density'
                rf'\s*\=\s*{RE_FLOAT}'
                rf'\;\s*Expected\s*\=\s*{RE_FLOAT}'
                rf'{RE_N}\s*average of density\,\s*n\s*\=\s*{RE_FLOAT}',
                repeats=False,
            ),
            Quantity(
                'wigner_seitz_radius', rf'{RE_N}\s*r_s\s*\=\s*{RE_FLOAT}', repeats=False
            ),
            Quantity(
                'omega_plasma',
                rf'{RE_N}\s*omega_plasma\s*\=\s*{RE_FLOAT}\s*\[(?P<__unit>\w+)\]',
                repeats=False,
            ),
        ]

        screening_quantities = rpa_quantities + [
            Quantity(
                'frequencies',
                rf'{RE_N}\s*(calculating chi0 at frequencies \[[a-zA-Z]+\]\s*\:'
                rf'[\s\S]+?)(?:\-*{RE_N}\s*q-point number\s*1)',
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'values',
                            rf'{RE_N}\s*\d*\s*{RE_FLOAT}\s*{RE_FLOAT}',
                            repeats=True,
                        )
                    ]
                ),
            ),
            Quantity(
                'static_diel_const',
                rf'{RE_N}\s*dielectric constant\s*\=\s*{RE_FLOAT}',
                repeats=False,
            ),
            Quantity(
                'static_diel_const_nofields',
                rf'{RE_N}\s*dielectric constant without local fields\s*\=\s*{RE_FLOAT}',
                repeats=False,
            ),
            Quantity(
                'chi_q',
                rf'{RE_N}(\s*q-point number\s*1[\s\S]+?)(?:\s*Average fulfillment'
                rf'[\s\w\[\]\-\:\.\%]*===)',
                repeats=False,
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'q_point',
                            rf'\s*q-point number\s*(\d+)\s*q =\s*\(\s*{RE_FLOAT}\,\s*'
                            rf'{RE_FLOAT}\,\s*{RE_FLOAT}\)\s*\[r\.l\.u\.\]',
                            repeats=True,
                        ),
                        Quantity(
                            'av_fulfillment',
                            rf'{RE_N}\s*Average fulfillment[a-zA-Z\s\[\]\-\d]*\:'
                            rf'\s*([\d\.]+)',
                            repeats=True,
                        ),
                        # TODO regex chi0(G, G') for each qpoint and frequency
                    ]
                ),
            ),
        ]

        self._quantities.append(
            Quantity(
                'screening_dataset',
                rf'{RE_N}(==\s*DATASET\s*3[\s\S]+?)(?:==\s*DATASET\s*4|== END DATASET)',
                repeats=False,
                sub_parser=TextParser(quantities=screening_quantities),
            )
        )

        def params_to_pairs(val_in):
            key = '_'.join(val_in.split()[:-1]).replace('-', '_')
            value = np.int(val_in.split()[-1])
            return [key, value]

        gw_quantities = rpa_quantities + [
            Quantity(
                'ks_band_gaps',
                rf'{RE_N}\s*(\>\>\>\>\s*For spin\s*1[\s\S]+?)(?:\s*SIGMA fundamental'
                rf' parameters)',
                repeats=False,
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'min_direct_gap',
                            rf'\s*Minimum direct gap\s*=\s*{RE_FLOAT}\s*\[\w*\]\,\s*'
                            rf'located at k-point\s*\:\s*{RE_FLOAT}\s*{RE_FLOAT}\s*'
                            rf'{RE_FLOAT}',
                            repeats=False,
                        ),
                        Quantity(
                            'fundamental_gap',
                            rf'\s*Fundamental gap\s*=\s*{RE_FLOAT}\s*'
                            rf'\[(?P<__unit>\w+)\]',
                            repeats=False,
                        ),
                        Quantity(
                            'k_top_valence_band',
                            rf'\s*Top of valence bands at\s*\:\s*{RE_FLOAT}\s*'
                            rf'{RE_FLOAT}\s*{RE_FLOAT}',
                            repeats=False,
                        ),
                        Quantity(
                            'k_bottom_conduction_band',
                            rf'\s*Bottom of conduction at\s*\:\s*{RE_FLOAT}\s*'
                            rf'{RE_FLOAT}\s*{RE_FLOAT}',
                            repeats=False,
                        ),
                    ]
                ),
            ),
            Quantity(
                'sigma_parameters',
                rf'{RE_N}(\s*SIGMA fundamental parameters[\s\S]+?)(?:\s*EPSILON\^\-1)',
                repeats=False,
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'model', r'([a-zA-Z\s]+) MODEL\s*(\d+)', repeats=False
                        ),
                        Quantity(
                            'params',
                            r'\s*number of ([a-zA-Z\s\-\/]+)(\d+)',
                            repeats=True,
                            str_operation=params_to_pairs,
                        ),
                        Quantity(
                            'freq_step',
                            rf'{RE_N}\s*frequency step[a-zA-Z\s\/]*\[eV\]\s*{RE_FLOAT}',
                            repeats=False,
                        ),
                        Quantity(
                            'max_omega_sigma',
                            rf'{RE_N}\s*max omega for Sigma[a-zA-Z\s]*\[eV\]\s*'
                            rf'{RE_FLOAT}',
                            repeats=False,
                        ),
                        Quantity(
                            'zcut_avoid',
                            rf'{RE_N}\s*zcut for avoiding poles\s*\[eV\]\s*{RE_FLOAT}',
                            repeats=False,
                        ),
                    ]
                ),
            ),
            Quantity(
                'epsilon_inv',
                rf'{RE_N}(\s*EPSILON\^\-1[\s\S]+?)(?:\s*Perturbative Calculation)',
                repeats=False,
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'dimensions',
                            r'\s*dimension of the eps\^\-1 matrix([a-zA-Z\s]+)(\d+)',
                            repeats=True,
                            str_operation=params_to_pairs,
                        ),
                        Quantity(
                            'params',
                            r'\s*number of ([a-zA-Z\s\-\/]+)(\d+)',
                            repeats=True,
                            str_operation=params_to_pairs,
                        ),
                    ]
                ),
            ),
            Quantity(
                'self_energy_ee',
                rf'{RE_N}(\-\-\-\s*\!SelfEnergy\_ee[\s\S]+?)(?:\.\.\.)',
                repeats=True,
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'kpoint',
                            rf'{RE_N}kpoint\s*\:[\s\[]*{RE_FLOAT}\, *{RE_FLOAT}\, *'
                            rf'{RE_FLOAT}',
                            repeats=False,
                        ),
                        Quantity(
                            'params',
                            rf'{RE_N}([a-zA-Z\_\s]+)\:\s*([\d\.\-\+]+)',
                            repeats=True,
                        ),
                        Quantity(
                            'data',
                            rf'{RE_N} *(\d+) *{RE_FLOAT} *{RE_FLOAT} *{RE_FLOAT} *'
                            rf'{RE_FLOAT} *{RE_FLOAT} *{RE_FLOAT} *{RE_FLOAT} *'
                            rf'{RE_FLOAT} *{RE_FLOAT}',
                            repeats=True,
                        ),
                    ]
                ),
            ),
        ]

        self._quantities.append(
            Quantity(
                'gw_dataset',
                rf'{RE_N}(==\s*DATASET\s*4[\s\S]+?)(?:==\s*DATASET\s*5|==\s*'
                rf'END DATASET\s*\(S\))',
                repeats=False,
                sub_parser=TextParser(quantities=gw_quantities),
            )
        )

    @property
    def dataset_numbers(self):
        if self._dataset_numbers is None:
            self._dataset_numbers = [
                d.get('x_abinit_dataset_number', 1) for d in self.get('dataset', [])
            ]
        return self._dataset_numbers

    @property
    def n_datasets(self):
        if self._n_datasets is None:
            self._n_datasets = max(self.dataset_numbers) if self.dataset_numbers else 1
        return self._n_datasets

    @property
    def input_vars(self):
        if self._input_vars is None:
            # set defaults
            self._input_vars = {
                key: [1] + [None] * (self.n_datasets - 1)
                for key in [
                    'ntypat',
                    'npsp',
                    'nshiftk',
                    'natrd',
                    'nsppol',
                    'nspden',
                    'nkpt',
                    'occopt',
                    'ixc',
                ]
            }

            for key_val in self.get('input_variables', {}).get('key_value', []):
                key, n_dataset = re.search(r'(\D+)?(\d*)', key_val[0]).groups()
                self._input_vars.setdefault(key, [None] * self.n_datasets)

                # m_quantity = x_abinit_section_input.m_def.all_quantities.get(
                #     f'x_abinit_var_{key}'
                # )
                # if m_quantity is None:
                #     continue

                val = key_val
                if '-' in key_val:  # exception when the next line starts with -
                    val = key_val[:-1]

                val = val[1:]
                # val = np.array(
                #     key_val[1:],
                #     dtype=m_quantity.type.standard_type()
                #     if hasattr(m_quantity.type, 'standard_type')
                #     else m_quantity.type,
                # )
                # if not m_quantity.shape:
                #     val = val[0]
                if n_dataset:
                    self._input_vars[key][int(n_dataset) - 1] = val
                else:
                    self._input_vars[key] = [val] * self.n_datasets
        return self._input_vars

    def get_input_var(self, key, n_dataset, default=None):
        val = self.input_vars.get(key)
        if val is None or val[n_dataset - 1] is None:
            val = [default] * n_dataset
        return val[n_dataset - 1]

    def get_atom_labels(self):
        znucl = self.get_input_var('znucl', 1)
        typat = self.get_input_var('typat', 1)
        if znucl is None or typat is None:
            return
        return [chemical_symbols[int(znucl[n_at - 1])] for n_at in typat]
