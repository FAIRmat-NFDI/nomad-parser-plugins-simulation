import re

from nomad_file_parser.text_parser import Quantity, TextParser

RE_N = r'[\n\r]'


class WInParser(TextParser):
    def init_quantities(self):
        def str_proj_to_list(val_in: str) -> list[str]:
            # To avoid inconsistent regex that can contain or not spaces
            val_n = [re.sub(r'\s.*', '', x) for x in val_in.split('\n') if x]
            return [v.strip('[]').replace(' ', '').split(':') for v in val_n]

        self._quantities = [
            Quantity(
                'energy_fermi', r'\n\rfermi_energy\s*=\s*([\d\.\-]+)', repeats=False
            ),
            Quantity(
                'projections',
                r'[bB]egin [pP]rojections([\s\S]+?)(?:[eE]nd [pP]rojections)',
                repeats=False,
                str_operation=str_proj_to_list,
            ),
        ]


class HrParser(TextParser):
    def init_quantities(self):
        self._quantities = [
            Quantity('degeneracy_factors', r'\s*written on[\s\w]*:\d*:\d*\s*([\d\s]+)'),
            Quantity('hoppings', r'\s*([-\d\s.]+)', repeats=False),
        ]


class WOutParser(TextParser):
    def __init__(self):
        super().__init__(None)

    def init_quantities(self):
        kmesh_quantities = [
            Quantity('n_points', r'Total points[\s=]*(\d+)', dtype=int, repeats=False),
            Quantity(
                'grid', r'Grid size *\= *(\d+) *x *(\d+) *x *(\d+)', repeats=False
            ),
            Quantity('k_points', r'\|[\s\d]*(-*\d.[^\|]+)', repeats=True, dtype=float),
        ]

        klinepath_quantities = [
            Quantity(
                'high_symm_name',
                r'\| *From\: *([a-zA-Z]+) [\d\.\-\s]*To\: *([a-zA-Z]+)',
                repeats=True,
            ),
            Quantity(
                'high_symm_value',
                r'\| *From\: *[a-zA-Z]* *([\d\.\-\s]+)To\: *[a-zA-Z]* *([\d\.\-\s]+)\|',
                repeats=True,
            ),
        ]

        disentangle_quantities = [
            Quantity(
                'outer',
                r'\|\s*Outer:\s*([-\d.]+)\s*\w*\s*([-\d.]+)\s*\((?P<__unit>\w+)\)',
                dtype=float,
                repeats=False,
            ),
            Quantity(
                'inner',
                r'\|\s*Inner:\s*([-\d.]+)\s*\w*\s*([-\d.]+)\s*\((?P<__unit>\w+)\)',
                dtype=float,
                repeats=False,
            ),
        ]

        structure_quantities = [
            Quantity('labels', r'\|\s*([A-Z][a-z]*)', repeats=True),
            Quantity(
                'positions',
                r'\|\s*([\-\d\.]+)\s*([\-\d\.]+)\s*([\-\d\.]+)',
                repeats=True,
                dtype=float,
            ),
        ]

        self._quantities = [
            # Program quantities
            Quantity('version', r'\s*\|\s*Release\:\s*([\d\.]+)\s*', repeats=False),
            # System quantities
            Quantity('lattice_vectors', r'\s*a_\d\s*([\d\-\s\.]+)', repeats=True),
            Quantity(
                'reciprocal_lattice_vectors', r'\s*b_\d\s*([\d\-\s\.]+)', repeats=True
            ),
            Quantity(
                'structure',
                rf'(\s*Fractional Coordinate[\s\S]+?)(?:{RE_N}\s*'
                rf'(PROJECTIONS|K-POINT GRID))',
                repeats=False,
                sub_parser=TextParser(quantities=structure_quantities),
            ),
            # Method quantities
            Quantity(
                'k_mesh',
                r'\s*(K-POINT GRID[\s\S]+?)(?:-\s*MAIN)',
                repeats=False,
                sub_parser=TextParser(quantities=kmesh_quantities),
            ),
            Quantity(
                'k_line_path',
                r'\s*(K-space path sections\:[\s\S]+?)(?:\*-------)',
                repeats=False,
                sub_parser=TextParser(quantities=klinepath_quantities),
            ),
            Quantity(
                'Nwannier',
                r'\|\s*Number of Wannier Functions\s*\:\s*(\d+)',
                repeats=False,
            ),
            Quantity(
                'Nband',
                r'\|\s*Number of input Bloch states\s*\:\s*(\d+)',
                repeats=False,
            ),
            Quantity(
                'Niter', r'\|\s*Total number of iterations\s*\:\s*(\d+)', repeats=False
            ),
            Quantity(
                'conv_tol',
                r'\|\s*Convergence tolerence\s*\:\s*([\d.eE-]+)',
                repeats=False,
            ),
            Quantity(
                'energy_windows',
                r'(\|\s*Energy\s*Windows\s*\|[\s\S]+?)'
                r'(?:Number of target bands to extract:)',
                repeats=False,
                sub_parser=TextParser(quantities=disentangle_quantities),
            ),
            # Band related quantities
            Quantity(
                'n_k_segments',
                r'\|\s*Number of K-path sections\s*\:\s*(\d+)',
                repeats=False,
            ),
            Quantity(
                'div_first_k_segment',
                r'\|\s*Divisions along first K-path section\s*\:\s*(\d+)',
                repeats=False,
            ),
            Quantity(
                'band_segments_points',
                r'\|\s*From\:\s*\w+([\d\s\-\.]+)To\:\s*\w+([\d\s\-\.]+)',
                repeats=True,
            ),
        ]
