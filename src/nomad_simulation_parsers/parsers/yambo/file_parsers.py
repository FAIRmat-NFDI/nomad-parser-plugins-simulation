import numpy as np
from netCDF4 import Dataset  # pylint: disable=no-name-in-module
from nomad.parsing.file_parser import FileParser
from nomad.parsing.file_parser.text_parser import DataTextParser, Quantity, TextParser
from nomad.units import ureg

RE_FLOAT = r'[-+]*\d*\.\d+[Ee]*[-+]*\d*'


class SpectraParser(DataTextParser):  # EM: moved sp_type here (from MainfileParser),  Jul 1st, 2026
    def init_quantities(self) -> None:
        self._quantities = [        
            Quantity(
                'sp_type',
                r'(EELS|Polarizability|Absorption)',
                repeats=False,
            ),
            #HB, Aug 3rd, 2026
            Quantity(
                'n_energies',
                r'(?:BEnSteps=|ETStpsXd=)\s*(\d+)',
                dtype=np.int32,
            ),
            end HB
            # EM, Jul 23rd, 2026
          #  Quantity(
           #     'n_energies',
            #    r'BEnSteps=|ETStpsXd=|\s+(\d+)',
             #   dtype=np.int32,
            # ),
            # end EM
        ]    
#    pass


class InputParser(TextParser):
    def init_quantities(self) -> None:
        def str_to_key_block(val_in: str) -> tuple[str, list[np.ndarray]]:
            val = val_in.strip().split('\n')
            return val[0].strip(), [
                np.array(v.split('#')[0].split('|')[:-1], dtype=np.float64)
                for v in val[1:]
            ]

        self._quantities = [
            Quantity('key_value', r'\n *(\w+) *= *(.+) *#*', repeats=True),
            Quantity(
                'key_block',
                r'\n *\% *(\w+)([\s\S]+?)\%',
                repeats=True,
                str_operation=str_to_key_block,
            ),
        ]


class NetCDFParser(FileParser):
    def init_parameters(self) -> None:
        self._keys = []

    @property
    def netcdf_file(self):
        if self._file_handler is None:
            try:
                self._file_handler = Dataset(self.mainfile)
            except Exception:
                self.logger.warning('Error loading file.')
                raise

        return self._file_handler

    def parse(self, key=None) -> None:
        self._results = dict() if self._results is None else self._results
        if self.netcdf_file is None:
            return

        self._keys = list(self.netcdf_file.variables.keys())
        for netcdf_variable in self._keys:
            self._results[netcdf_variable] = self.netcdf_file.variables[
                netcdf_variable
            ][:].data


class MainfileParser(TextParser):
    def init_quantities(self) -> None:
        io_quantities = [
            Quantity(
                'key_value',
                r'([A-Z\d].+?)(?:\(.+\)|\[.+\]| |)(:.+?)(?:\[|\n)',
                str_operation=lambda x: [v.strip() for v in x.split(':')],
                repeats=True,
            ),
            Quantity('file', r'\[(?:RD|WR)(.+?)\]', dtype=str),
            Quantity('sn', r'- S/N *(\d+)', dtype=str),
        ]

        energies_quantities = [
            Quantity(
                'fermi',
                rf'Fermi Level.+?: +({RE_FLOAT})',
                dtype=np.float64,
                unit=ureg.eV,
            ),
            Quantity(
                'conduction',
                rf'Conduction Band Min +: +({RE_FLOAT})',
                dtype=np.float64,
                unit=ureg.eV,
            ),
            Quantity(
                'valence',
                rf'Valence Band Max +: +({RE_FLOAT})',
                dtype=np.float64,
                unit=ureg.eV,
            ),
            Quantity(
                'valence_conduction',
                rf'VBM / CBm +\[ev\]: +({RE_FLOAT}) +({RE_FLOAT})',
                dtype=np.dtype(np.float64),
                unit=ureg.eV,
            ),
            Quantity(
                'x_yambo_filled_bands',
                r'Filled Bands +: +(\d+)',
                dtype=np.int32,
                str_operation=lambda x: [1, int(x)],
            ),
            Quantity(
                'x_yambo_empty_bands',
                r'Empty Bands +: +([\d ]+)',
                dtype=np.dtype(np.int32),
            ),
            Quantity(
                'x_yambo_electronic_temperature',
                rf'Electronic Temp.+?: +{RE_FLOAT} +({RE_FLOAT})',
                dtype=np.float64,
                unit=ureg.kelvin,
            ),
            Quantity(
                'x_yambo_bosonic_temperature',
                rf'Bosonic +Temp.+?: +{RE_FLOAT} +({RE_FLOAT})',
                dtype=np.float64,
                unit=ureg.kelvin,
            ),
            Quantity(
                'x_yambo_finite_temperature_mode',
                r'Finite Temperature mode: +(\S+)',
                str_operation=lambda x: x == 'yes',
            ),
            Quantity(
                'x_yambo_electronic_density',
                r'El\. density.+?: +(.+?)(?:\[|\n)',
                str_operation=lambda x: x.strip().split()[-1],
                dtype=np.float64,
            ),
            Quantity(
                'states_summary',
                r'States summary +: Full +Metallic +Empty\s+(.+)',
                str_operation=lambda x: [v.split('-') for v in x.strip().split()],
            ),
            Quantity(
                'x_yambo_indirect_gaps',
                rf'Indirect Gaps.+?: +({RE_FLOAT}) +({RE_FLOAT})',
                dtype=np.dtype(np.float64),
                unit=ureg.eV,
            ),
            Quantity(
                'x_yambo_direct_gaps',
                rf'Direct Gaps.+?: +({RE_FLOAT}) +({RE_FLOAT})',
                dtype=np.dtype(np.float64),
                unit=ureg.eV,
            ),
            Quantity(
                'x_yambo_indirect_gap',
                rf'Indirect Gap.+?: +({RE_FLOAT})',
                dtype=np.float64,
                unit=ureg.eV,
            ),
            Quantity(
                'x_yambo_direct_gap',
                rf'Direct Gap.+?: +({RE_FLOAT})',
                dtype=np.float64,
                unit=ureg.eV,
            ),
            Quantity(
                'x_yambo_direct_gap_kpoint',
                r'Direct Gap localized at k-point.+?: +(\d+)',
                dtype=np.int32,
            ),
            Quantity(
                'x_yambo_indirect_gap_kpoints',
                r'Indirect Gap between k-points.+?: +(\d+) +(\d+)',
                dtype=np.int32,
            ),
        ]

        qp_properties_quantity = Quantity(
            'qp_properties',
            r'QP properties and I/O([\s\S]+? S/N \d+.+)',
            sub_parser=TextParser(
                quantities=[
                    Quantity(
                        'qp_energy',
                        r'(QP \[eV\] @ K[\s\S]+?)\n *\n',
                        repeats=True,
                        sub_parser=TextParser(
                            quantities=[
                                Quantity(
                                    'band',
                                    rf'B= *(\d+) Eo= *({RE_FLOAT}) E= *({RE_FLOAT}) '
                                    rf'E-Eo= *({RE_FLOAT}) '
                                    rf'Re\(Z\)= *({RE_FLOAT}) Im\(Z\)= *({RE_FLOAT}) '
                                    rf'nlXC= *({RE_FLOAT}) lXC= *({RE_FLOAT}) '
                                    rf'So= *({RE_FLOAT})',
                                    repeats=True,
                                    dtype=np.dtype(np.float64),
                                ),
                                Quantity(
                                    'kpoint',
                                    r'K *\[\d+\].+?\: *(.+)',
                                    dtype=np.dtype(np.float64),
                                ),
                            ]
                        ),
                    ),
                    Quantity(
                        'output',
                        r'(\[WR.+?\.QP\][\s\S]+?- S/N \d+.+)',
                        repeats=True,
                        sub_parser=TextParser(quantities=io_quantities),
                    ),
                ]
            ),
        )

        module_quantities = [
            Quantity(
                'dipoles',
                r'Dipoles *\n([\s\S]+?)\n *\[\d+\]',
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'input',
                            r'(\[RD.+?\][\s\S]+?- S/N \d+.+)',
                            repeats=True,
                            sub_parser=TextParser(quantities=io_quantities),
                        ),
                        Quantity(
                            'output',
                            r'(\[WR.+?\.dipoles\][\s\S]+?- S/N \d+.+)',
                            repeats=True,
                            sub_parser=TextParser(quantities=io_quantities),
                        ),
                    ]
                ),
            ),
            Quantity(
                'local_xc_nonlocal_fock',
                r'Local Exchange-Correlation \+ '
                r'Non-Local Fock([\s\S]+?(?:\n *\[\d+\]|\Z))',
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'output',
                            r'(\[WR.+?\.HF_and_locXC\][\s\S]+?- S/N \d+.+)',
                            repeats=True,
                            sub_parser=TextParser(quantities=io_quantities),
                        ),
                        Quantity(
                            'x_yambo_plane_waves_vxc',
                            r'\[VXC\] Plane waves : *(\d+)',
                            dtype=np.int32,
                        ),
                        Quantity(
                            'x_yambo_plane_waves_exs',
                            r'\[EXS\] Plane waves : *(\d+)',
                            dtype=np.int32,
                        ),
                        Quantity(
                            'x_yambo_mesh_size',
                            r'Mesh size: *(\d+) *(\d+) *(\d+)',
                            dtype=np.dtype(np.int32),
                        ),
                        Quantity(
                            'energy_xc',
                            rf'E_xc *: *({RE_FLOAT}) \[Ha\]',
                            dtype=np.float64,
                            unit='hartree',
                        ),
                        Quantity(
                            'corrections',
                            r'Corrections @ K \[\d+\] *: *\[eV\]([\s\S]+?)\n *\n',
                            repeats=True,
                            sub_parser=TextParser(
                                quantities=[
                                    Quantity(
                                        'band',
                                        rf'\<\d+\|nlXC\|\d+\> *= *({RE_FLOAT}) *'
                                        rf'{RE_FLOAT} \<\d+\|lXC\|\d+\> *= *'
                                        rf'({RE_FLOAT}) *{RE_FLOAT}',
                                        repeats=True,
                                        dtype=np.dtype(np.float64),
                                    ),
                                    Quantity(
                                        'band_sp',
                                        rf'\<\d+\((?:up|dn)\)\|nlXC\|\d+\((?:up|dn)\)\>'
                                        rf' *= *({RE_FLOAT}) *{RE_FLOAT} '
                                        rf'\<\d+\((?:up|dn)\)\|lXC\|\d+\((?:up|dn)\)\>'
                                        rf' *= *({RE_FLOAT}) *{RE_FLOAT}',
                                        repeats=True,
                                        dtype=np.dtype(np.float64),
                                    ),
                                ]
                            ),
                        ),
                        Quantity(
                            'hf_occupations',
                            r'Hartree-Fock occupations report([\s\S]+?)'
                            r'(?:\n *\[\d+|\Z)',
                            sub_parser=TextParser(quantities=energies_quantities),
                        ),
                    ]
                ),
            ),
            # TODO add support for em1d
            Quantity(
                'dynamic_dielectric_matrix',
                r'Dynamic.+?Dielectric Matrix([\s\S]+?(?:\n *\[\d+\]|\Z))',
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'output',
                            r'(\[WR.+?(?:pp|em1d|dip_iR_and_P)\][\s\S]+?- S/N \d+.+)',
                            repeats=True,
                            sub_parser=TextParser(quantities=io_quantities),
                        ),
                        Quantity(
                            'x_yambo_mesh_size',
                            r'Mesh size: *(\d+) *(\d+) *(\d+)',
                            dtype=np.dtype(np.int32),
                        ),
                    ]
                ),
            ),
            Quantity(
                'bare_xc',
                r'Bare local and non-local Exchange-Correlation'
                r'([\s\S]+?(?:\n *\[\d+\]|\Z))',
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'output',
                            r'(\[WR.+?\.HF_and_locXC\][\s\S]+?- S/N \d+.+)',
                            repeats=True,
                            sub_parser=TextParser(quantities=io_quantities),
                        ),
                        Quantity(
                            'xc_hf_dft',
                            r'XC HF and DFT \[eV\]([\s\S]+?)\n *\n',
                            repeats=True,
                            sub_parser=TextParser(
                                quantities=[
                                    Quantity(
                                        'band',
                                        rf'\<\d+\|HF\|\d+\> *= *({RE_FLOAT}) *'
                                        rf'{RE_FLOAT} *\<\d+\|DFT\|\d+\> *= *'
                                        rf'({RE_FLOAT}) *{RE_FLOAT}',
                                        repeats=True,
                                        dtype=np.dtype(np.float64),
                                    )
                                ]
                            ),
                        ),
                        Quantity(
                            'hf_occupations',
                            r'HF occupations report([\s\S]+?Direct Gaps.+)',
                            sub_parser=TextParser(quantities=energies_quantities),
                        ),
                    ]
                ),
            ),
            Quantity(
                'dyson',
                r'Dyson equation: Newton solver([\s\S]+?(?:\n *\[\d+\]|\Z))',
                sub_parser=TextParser(
                    quantities=[
                        qp_properties_quantity,
                        Quantity(
                            'g0w0',
                            r'G0W0([\s\S]+?\n *\[\d+\.\d+\])',
                            sub_parser=TextParser(
                                quantities=[
                                    Quantity(
                                        'x_yambo_bands_range',
                                        r'Bands range *: *(\d+) *(\d+)',
                                        dtype=np.dtype(np.int32),
                                    ),
                                    Quantity(
                                        'x_yambo_g_damping',
                                        rf'G damping.+?: *({RE_FLOAT})',
                                        dtype=np.float64,
                                    ),
                                    Quantity(
                                        'x_yambo_mesh_size',
                                        r'Mesh size: *(\d+) *(\d+) *(\d+)',
                                        dtype=np.dtype(np.int32),
                                    ),
                                    Quantity(
                                        'input',
                                        r'(\[RD.+?\.pp\][\s\S]+?- S/N \d+.+)',
                                        repeats=True,
                                        sub_parser=TextParser(quantities=io_quantities),
                                    ),
                                ]
                            ),
                        ),
                    ]
                ),
            ),
        ]

        self._quantities = [
            Quantity(
                'version', r'Version ([\d.]+ Revision \d+)', flatten=False, dtype=str
            ),
            Quantity('hash', r'Hash (\S+)', dtype=str),
            Quantity('build', r'(\S+) Build', dtype=str),
            Quantity(
                'date_start',
                r' (\d\d/\d\d/\d\d\d\d) at (\d\d:\d\d) YAMBO @ .+',
                flatten=False,
                dtype=str,
            ),
            Quantity(
                'cpu_files_io',
                r'((?:Cores |CPU structure)[\s\S]+?)\n *\[\d+\]',
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'parameters',
                            r'([A-Z][\w/ ]+).*?(?: is | in |: ) *(\S+)',
                            repeats=True,
                            str_operation=lambda x: [
                                v.strip() for v in x.rsplit(' ', 1)
                            ],
                        ),
                        Quantity(
                            'input',
                            r'( \[RD.+[\s\S]+?- S/N \d+.+)',
                            repeats=False,
                            sub_parser=TextParser(quantities=io_quantities),
                        ),
                    ]
                ),
            ),
            Quantity(
                'core_variables_setup',
                r'(CORE Variables Setup[\s\S]+?)\n *\[\d+\]',
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'energies_occupations',
                            r'Energies.+?& Occupations([\s\S]+?)(?:\[0|\Z)',
                            sub_parser=TextParser(
                                quantities=energies_quantities
                                + [
                                    Quantity(
                                        'eigenenergies',
                                        rf'Energy unit is electronVolt \[eV\]([\s\S]+'
                                        rf'E *{RE_FLOAT} *{RE_FLOAT} *{RE_FLOAT}.+)',
                                        sub_parser=TextParser(
                                            quantities=[
                                                Quantity(
                                                    'energies',
                                                    rf'\n *E *({RE_FLOAT} .+)',
                                                    repeats=True,
                                                    dtype=np.dtype(np.float64),
                                                ),
                                                Quantity(
                                                    'kpoints',
                                                    rf'({RE_FLOAT} *{RE_FLOAT} *'
                                                    rf'{RE_FLOAT}) \(rlu\)',
                                                    repeats=True,
                                                    dtype=np.dtype(np.float64),
                                                ),
                                                Quantity(
                                                    'kpoints_weights',
                                                    rf'weight +({RE_FLOAT})',
                                                    repeats=True,
                                                    dtype=np.float64,
                                                ),
                                            ]
                                        ),
                                    )
                                ]
                            ),
                        )
                    ]
                ),
            ),
#            Quantity(
#                'sp_type',
#                r'(EELS|Polarizability|Absorption)',
#                repeats=False,
#            ),
            Quantity(
                'transferred_momenta',
                r'Transferred momenta grid([\s\S]+?)\n *\[\d+\]',
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'input',
                            r'( \[RD.+[\s\S]+?- S/N \d+.+)',
                            repeats=True,
                            sub_parser=TextParser(quantities=io_quantities),
                        ),
                        Quantity(
                            'qpoints',
                            rf'Q \[\d+\] *: *({RE_FLOAT}) *({RE_FLOAT}) *({RE_FLOAT})'
                            rf' *\(iku\) \* weight *({RE_FLOAT})',
                            repeats=True,
                            dtype=np.dtype(np.float64),
                        ),
                        Quantity(
                            'module',
                            r'((?:Dipoles *\n|Dynamic Dielectric|Dyson|Bare local|'
                            r'Local Exchange)[\s\S]+?\n *\[\d+\.\d+\])',
                            repeats=True,
                            sub_parser=TextParser(quantities=module_quantities),
                        ),
                        qp_properties_quantity,
                    ]
                ),
            ),
            Quantity(
                'module',
                r'((?:Dipoles *\n|Dynamic.+?Dielectric|Dyson|Bare local|'
                r'Local Exchange)[\s\S]+?\n *\[\d+\])',
                repeats=True,
                sub_parser=TextParser(quantities=module_quantities),
            ),
        ]
