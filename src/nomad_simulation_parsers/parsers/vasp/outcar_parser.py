from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    pass

import re

import numpy as np
from nomad.parsing.file_parser import ArchiveWriter, Quantity, TextParser
from nomad.parsing.file_parser.mapping_parser import MetainfoParser, Path
from nomad.parsing.file_parser.mapping_parser import TextParser as MappingTextParser
from nomad.units import ureg
from nomad.utils import get_logger
from nomad_simulations.schema_packages.general import Simulation
from nomad_simulations.schema_packages.model_method import XCFunctional
from nomad_simulations.schema_packages.numerical_settings import PPCutoff

from nomad_simulation_parsers.schema_packages import vasp

RE_N = r'[\n\r]'
LOGGER = get_logger(__name__)


def get_key_values(val_in):
    val = [v for v in val_in.split('\n') if '=' in v]
    data = {}
    pattern = re.compile(r'([A-Z_]+)\s*=\s*(\.?[a-zA-Z]*[\d\-\.\+\sE]*\.?)')

    def convert(v):
        if isinstance(v, list):
            v = [convert(vi) for vi in v]
        elif isinstance(v, str):
            try:
                v = float(v) if '.' in v else int(v)
            except Exception:
                pass
        else:
            pass
        return v

    for v in val:
        res = pattern.findall(v)
        for resi in res:
            vi = resi[1].split()
            vi = vi[0] if len(vi) == 1 else vi
            if isinstance(vi, str):
                vi = vi.strip()
                vi_upper = vi.upper()
                if vi_upper in ['T', '.TRUE.', 'TRUE']:
                    vi = True
                elif vi_upper in ['F', '.FALSE.', 'FALSE']:
                    vi = False
            data[resi[0]] = convert(vi)
    return data


# TODO temporary fix for structlog unable to propagate logger
class VASPMetainfoParser(MetainfoParser):
    @property
    def logger(self):
        return LOGGER


class OutcarTextParser(TextParser):
    def __init__(self):
        self._chemical_symbols = None

        super().__init__(None)

    def init_quantities(self):  # noqa: PLR0915
        def str_to_array(val_in):
            val = [
                re.findall(r'(\-?\d+\.[\dEe]+)', v)
                for v in val_in.strip().split('\n')
                if '--' not in v
            ]
            return [
                np.array([v[0:3] for v in val], float),
                np.array([v[3:6] for v in val], float),
            ]

        def str_to_stress(val_in):
            val = [float(v) for v in val_in.strip().split()]
            stress = np.zeros((3, 3))
            stress[0][0] = val[0]
            stress[1][1] = val[1]
            stress[2][2] = val[2]
            stress[0][1] = stress[1][0] = val[3]
            stress[1][2] = stress[2][1] = val[4]
            stress[0][2] = stress[2][0] = val[5]
            return stress

        def str_to_header(val_in):
            (
                version,
                build_date,
                build_type,
                platform,
                date,
                time,
                parallel,
            ) = val_in.split()
            parallel = 'parallel' if parallel == 'running' else parallel
            subversion = ' '.join([build_date, build_type, parallel])
            date = date.replace('.', ' ')
            return dict(
                version=version,
                subversion=subversion,
                platform=platform,
                date=date,
                time=time,
            )

        def str_to_positions(val_in):
            re_position = re.compile(
                r'\d*\s*(\-*\d+\.\d+)\s*(\-*\d+\.\d+)\s*(\-*\d+\.\d+)'
            )
            positions = []
            for val in val_in.strip().split('\n'):
                position = re_position.search(val)
                if position:
                    positions.append(position.groups())
            return np.array(positions, dtype=float)

        def str_to_eigenvalues(val_in):
            val = []
            for line in val_in.strip().splitlines():
                val.extend(['nan' if '*' in v else v for v in line.split()])
            return np.array(val, np.float64)

        scf_iteration = [
            Quantity(
                'energy_total',
                r'free energy\s*TOTEN\s*=\s*([\d\.\-]+)\s*eV',
                repeats=False,
                dtype=float,
            ),
            Quantity(
                'energy_entropy0',
                r'energy without entropy\s*=\s*([\d\.\-]+)',
                repeats=False,
                dtype=float,
            ),
            Quantity(
                'energy_T0',
                r'energy\(sigma\->0\)\s*=\s*([\d\.\-]+)',
                repeats=False,
                dtype=float,
            ),
            Quantity(
                'energy_components',
                r'Free energy of the ion-electron system \(eV\)\s*\-+([\s\S]+?)\-{10}',
                str_operation=get_key_values,
                convert=False,
            ),
            Quantity(
                'time',
                r'LOOP\: +cpu time +([\d\.]+)\: +real time +([\d\.]+)',
                dtype=np.dtype(np.float64),
            ),
        ]

        calculation_quantities = [
            Quantity(
                'scf_iteration',
                r'Iteration\s*\d+\(\s*\d+\s*\)([\s\S]+?energy\(sigma\->0\)\s*=\s*.+)',
                repeats=True,
                sub_parser=TextParser(quantities=scf_iteration),
            ),
            Quantity(
                'energies',
                r'FREE ENERGIE OF THE ION-ELECTRON SYSTEM \(eV\)'
                r'\s*\-+\s*([\s\S]+?)\-{10}',
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'energy_total',
                            r'free\s*energy\s*TOTEN\s*=\s*([\-\d\.]+)',
                            repeats=False,
                            dtype=float,
                        ),
                        Quantity(
                            'energy_entropy0',
                            r'energy\s*without\s*entropy\s*=\s*([\-\d\.]+)',
                            repeats=False,
                            dtype=float,
                        ),
                        Quantity(
                            'energy_T0',
                            r'energy\(sigma\->0\)\s*=\s*([\-\d\.]+)',
                            repeats=False,
                            dtype=float,
                        ),
                    ]
                ),
            ),
            Quantity(
                'stress',
                r'in kB\s*(\-?\d+\.\d+)\s*(\-?\d+\.\d+)\s*(\-?\d+\.\d+)\s*'
                r'(\-?\d+\.\d+\s*)(\-?\d+\.\d+\s*)(\-?\d+\.\d+\s*)',
                str_operation=str_to_stress,
                convert=False,
            ),
            Quantity(
                'positions_forces',
                r'POSITION\s*TOTAL\-FORCE \(eV/Angst\)\s*\-+\s*([\d\.\s\-E]+)',
                str_operation=str_to_array,
                convert=False,
            ),
            Quantity(
                'lattice_vectors',
                r'direct lattice vectors\s*reciprocal lattice vectors\s*'
                r'(\-?\d+\.\d+\s*)(\-?\d+\.\d+\s*)(\-?\d+\.\d+\s*)(\-?\d+\.\d+\s*)(\-?\d+\.\d+\s*)(\-?\d+\.\d+\s*)'
                r'(\-?\d+\.\d+\s*)(\-?\d+\.\d+\s*)(\-?\d+\.\d+\s*)(\-?\d+\.\d+\s*)(\-?\d+\.\d+\s*)(\-?\d+\.\d+\s*)'
                r'(\-?\d+\.\d+\s*)(\-?\d+\.\d+\s*)(\-?\d+\.\d+\s*)(\-?\d+\.\d+\s*)(\-?\d+\.\d+\s*)(\-?\d+\.\d+\s*)',
                str_operation=str_to_array,
                convert=False,
            ),
            Quantity(
                'converged',
                r'aborting loop because (EDIFF is reached)',
                repeats=False,
                dtype=str,
                convert=False,
            ),
            Quantity(
                'fermi_energy', r'E\-fermi :\s*([\d\.]+)', dtype=str, repeats=False
            ),
            Quantity(
                'eigenvalues',
                r'band No\.\s*band energies\s*occupation\s*([\d\.\s\-\*]+?)'
                r'(?:k\-point|spin|\-{10})',
                repeats=True,
                dtype=float,
                str_operation=str_to_eigenvalues,
            ),
            Quantity('convergence', r'(aborting loop because EDIFF is reached)'),
            Quantity(
                'time',
                r'LOOP\+\: +cpu time +([\d\.]+)\: +real time +([\d\.]+)',
                dtype=np.dtype(np.float64),
            ),
        ]

        self._quantities = [
            Quantity(
                'calculation',
                r'(\-\-\s*Iteration\s*\d+\(\s*1\s*\)\s*[\s\S]+?)'
                r'((?:FREE ENERGIE OF THE ION\-ELECTRON SYSTEM \(eV\)[\s\S]+?'
                r'LOOP\+.+)|\Z)',
                repeats=True,
                sub_parser=TextParser(quantities=calculation_quantities),
            ),
            Quantity(
                'header',
                r'vasp\.([\d\.]+)\s*(\w+)\s*[\s\S]+?\)\s*(\w+)\s*'
                r'executed on\s*(\w+)\s*date\s*([\d\.]+)\s*([\d\:]+)\s*(\w+)',
                repeats=False,
                str_operation=str_to_header,
                convert=False,
            ),
            Quantity(
                'parameters',
                r'Startparameter for this run:([\s\S]+?)\-{100}',
                str_operation=get_key_values,
                repeats=False,
                convert=False,
            ),
            Quantity(
                'ions_per_type', r'ions per type =\s*([ \d]+)', dtype=int, repeats=False
            ),
            Quantity(
                'species',
                r'(\w+) +([A-Z][a-z]*).+?:\s*energy of atom +\d+',
                dtype=str,
                repeats=True,
            ),  # TODO: deprecate
            Quantity(
                'pseudopotentials',
                r'(POTCAR:\s*.+?[\s\S]+?)(?=POTCAR:|end of INCAR parameters|\Z)',
                repeats=True,
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'titel',
                            r'TITEL\s*=\s*(.+)',
                            dtype=str,
                            repeats=False,
                        ),
                        Quantity(
                            'vrhfin',
                            r'VRHFIN\s*=\s*(.+)',
                            str_operation=lambda x: (
                                x.split(':', 1)[1].strip() if ':' in x else x.strip()
                            ),
                            repeats=False,
                        ),
                        Quantity(
                            'lexch', r'LEXCH\s*=\s*(\w+)', dtype=str, repeats=False
                        ),
                        Quantity(
                            'zval',
                            r'ZVAL\s*=\s*([\d\.]+)',
                            dtype=float,
                            repeats=False,
                        ),
                        Quantity(
                            'rcore',
                            r'RCORE\s*=\s*([\d\.]+)',
                            dtype=float,
                            repeats=False,
                        ),
                        Quantity(
                            'enmax',
                            r'ENMAX\s*=\s*([\d\.]+)',
                            dtype=float,
                            repeats=False,
                        ),
                        Quantity(
                            'enmin',
                            r'ENMIN\s*=\s*([\d\.]+)',
                            dtype=float,
                            repeats=False,
                        ),
                        Quantity(
                            'lpaw',
                            r'LPAW\s*=\s*([TF])',
                            str_operation=lambda x: x == 'T',
                            repeats=False,
                        ),
                        Quantity(
                            'lultra',
                            r'LULTRA\s*=\s*([TF])',
                            str_operation=lambda x: x == 'T',
                            repeats=False,
                        ),
                        Quantity(
                            'sha256',
                            r'SHA256\s*=\s*([a-f0-9]{64})',
                            dtype=str,
                            repeats=False,
                        ),
                        Quantity(
                            'lmax',
                            r'number of l-projection\s+operators is LMAX\s*=\s*(\d+)',
                            dtype=int,
                            repeats=False,
                        ),
                        Quantity(
                            'lmmax',
                            r'number of lm-projection operators is LMMAX\s*=\s*(\d+)',
                            dtype=int,
                            repeats=False,
                        ),
                    ]
                ),
                convert=False,
            ),
            Quantity(
                'kpoints',
                r'Following reciprocal coordinates:[\s\S]+?\n([\d\.\s\-]+)',
                repeats=False,
                dtype=float,
            ),
            Quantity('nbands', r'NBANDS\s*=\s*(\d+)', dtype=int, repeats=False),
            Quantity(
                'lattice_vectors',
                r'direct lattice vectors\s*reciprocal lattice vectors\s*'
                r'(\-?\d+\.\d+\s*)(\-?\d+\.\d+\s*)(\-?\d+\.\d+\s*)(\-?\d+\.\d+\s*)(\-?\d+\.\d+\s*)(\-?\d+\.\d+\s*)'
                r'(\-?\d+\.\d+\s*)(\-?\d+\.\d+\s*)(\-?\d+\.\d+\s*)(\-?\d+\.\d+\s*)(\-?\d+\.\d+\s*)(\-?\d+\.\d+\s*)'
                r'(\-?\d+\.\d+\s*)(\-?\d+\.\d+\s*)(\-?\d+\.\d+\s*)(\-?\d+\.\d+\s*)(\-?\d+\.\d+\s*)(\-?\d+\.\d+\s*)',
                str_operation=str_to_array,
                convert=False,
            ),
            Quantity(
                'positions',
                r'ion\s*position\s*nearest neighbor table([\s\S]+?)LATTYP',
                str_operation=str_to_positions,
                convert=False,
            ),
            # alternative format
            Quantity(
                'positions',
                r'position of ions in cartesian coordinates\s*\(Angst\):'
                r'([\s\S]+?)\n *\n',
                str_operation=str_to_positions,
                convert=False,
            ),
            Quantity(
                'response_functions',
                r'\s*Response functions by sum over occupied states\:'
                r'([\s\S]+?)(?:\-\-\-\-\-\-)',
                repeats=False,
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'input_parameters',
                            rf'{RE_N}* *(\w+) *\= *([\w\.\-]+) *.*',
                            repeats=True,
                        )
                    ]
                ),
            ),
        ]


class OutcarParser(MappingTextParser):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.xc_functional_mapping = {
            '--': ['GGA_X_PBE', 'GGA_C_PBE'],
            'HL': ['LDA_C_HL'],
            'WI': ['LDA_C_WIGNER'],
            'PZ': ['LDA_C_PZ'],
            '91': ['GGA_X_PW91', 'GGA_C_PW91'],
            'PE': ['GGA_X_PBE', 'GGA_C_PBE'],
            'PBE': ['GGA_X_PBE', 'GGA_C_PBE'],
            'RE': ['GGA_X_PBE_R'],
            'VW': ['LDA_C_VWN'],
            'RP': ['GGA_X_RPBE', 'GGA_C_PBE'],
            'PS': ['GGA_C_PBE_SOL', 'GGA_X_PBE_SOL'],
            'AM': ['GGA_X_AM05', 'GGA_C_AM05'],
            'B3': ['HYB_GGA_XC_B3LYP3'],
            'B5': ['HYB_GGA_XC_B3LYP5'],
            'BF': ['GGA_X_BEEFVDW', 'GGA_XC_BEEFVDW'],
            'CO': [],  # TODO check if this is ever used
            'OR': ['GGA_X_OPTPBE_VDW'],
            'BO': ['GGA_X_OPTB88_VDW'],
            'MK': ['GGA_X_OPTB86B_VDW'],
            'ML': ['VDW_XC_DF2'],
            'CX': ['VDW_XC_DF_CX'],
            'TPSS': ['MGGA_X_TPSS', 'MGGA_C_TPSS'],
            'RTPSS': ['MGGA_X_RTPSS'],
            'M06L': ['MGGA_C_M06_L'],
            'MS0': ['MGGA_X_MS0'],
            'MS1': ['MGGA_X_MS1'],
            'MS2': ['MGGA_X_MS2'],
            'SCAN': ['MGGA_X_SCAN'],
            'RSCAN': ['MGGA_X_RSCAN', 'MGGA_C_RSCAN'],
            'R2SCAN': ['MGGA_X_R2SCAN', 'MGGA_C_R2SCAN'],
            'SCANL': ['MGGA_X_SCANL', 'MGGA_C_SCANL'],
            'RSCANL': [],  # not in LibXC, nor any paper, just deorbitalized SCANL
            'R2SCANL': ['MGGA_X_R2SCANL', 'MGGA_C_R2SCANL'],
            'OFR2': [],
            'MBJ': ['MGGA_X_BJ06'],
            'LBMJ': [],  # TODO ask Miguel Marquez
            'HLE17': ['MGGA_XC_HLE17'],  # TODO check if this is ever used
            'RA': ['LDA_C_PW_RPA'],  # TODO check if this is ever used
        }

    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER

    def get_version(self, source: dict[str, Any]) -> str:
        return ' '.join(
            [
                source[key]
                for key in ['version', 'subversion', 'platform']
                if source.get(key)
            ]
        )

    def get_data(self, source: Any, **kwargs) -> Any:
        if isinstance(source, dict) and source.get('value') is not None:
            return source['value']
        path = kwargs.get('path')
        if path is None:
            return
        parser = Path(path=path)
        return parser.get_data(source)

    def get_forces(self, source: Any) -> dict[str, Any]:
        value = self.get_data(source, path='.positions_forces | [1]')
        if value is None:
            return {}
        return dict(forces=value, npoints=len(value), rank=[3])

    def get_energy_contributions(
        self, source: dict[str, Any], **kwargs
    ) -> list[dict[str, Any]]:
        exclude = kwargs.get('exclude', [])
        return [
            {'name': key, 'value': val}
            for key, val in source.items()
            if key not in exclude
        ]

    def get_eigenvalues(
        self, eigenvalues: np.ndarray, parameters: dict[str, Any]
    ) -> list[dict[str, Any]]:
        ispin = parameters.get('ISPIN', 1)
        n_kpts = len(eigenvalues) // ispin
        n_bands = len(eigenvalues[0]) // 3
        eigenvalues = np.reshape(eigenvalues, (ispin, n_kpts, n_bands, 3))
        data = []
        for nspin in range(ispin):
            eigs, occs = eigenvalues[nspin].T[1:3]
            data.append(
                dict(
                    eigenvalues=eigs.T,
                    occupations=occs.T,
                    n_bands=n_bands,
                    npoints=n_kpts,
                )
            )
        return data

    def get_pseudopotentials(self, pseudopotentials: list[Any]) -> list[dict[str, Any]]:
        """
        Extract and filter pseudopotential data from sub-parser results.

        The OUTCAR parser extracts 4 POTCAR entries but the first 2 are just
        short header lines without full data. This transformer converts
        TextParser objects to dictionaries and filters out empty ones.
        """
        result = []
        # List of quantity names we expect from the sub-parser
        quantity_names = [
            'titel',
            'vrhfin',
            'lexch',
            'zval',
            'rcore',
            'enmax',
            'enmin',
            'lpaw',
            'lultra',
            'sha256',
            'lmax',
            'lmmax',
        ]

        for pp in pseudopotentials:
            pp_dict = {}

            # Extract quantities from TextParser object using get() method
            if hasattr(pp, 'get'):
                for qty_name in quantity_names:
                    value = pp.get(qty_name)
                    if value is not None:
                        # Convert numpy arrays to strings for titel
                        if (
                            qty_name == 'titel'
                            and hasattr(value, '__iter__')
                            and not isinstance(value, str)
                        ):
                            value = ' '.join(str(v) for v in value)
                        pp_dict[qty_name] = value
            elif isinstance(pp, dict):
                pp_dict = pp

            # Only include non-empty dictionaries
            if pp_dict:
                result.append(pp_dict)

        return result

    def get_xc_functionals(self, parameters: dict[str, Any]) -> list[dict[str, Any]]:
        xc_functionals = []
        if parameters.get('LHFCALC', False):
            hfscreen_06, hfscreen_03 = 0.2, 0.3
            aexx_b3, aggax_b3, aggac_b3, aldac_b3 = 0.2, 0.72, 0.81, 0.19
            xc_functional = {}
            gga = parameters.get('GGA', 'PE')
            aexx = parameters.get('AEXX', 0.0)
            aggax = parameters.get('AGGAX', 1.0)
            aggac = parameters.get('AGGAC', 1.0)
            aldac = parameters.get('ALDAC', 1.0)
            hfscreen = parameters.get('HFSCREEN', 0.0)

            if hfscreen == hfscreen_06:
                xc_functional['name'] = 'HYB_GGA_XC_HSE06'
            elif hfscreen == hfscreen_03:
                xc_functional['name'] = 'HYB_GGA_XC_HSE03'
            elif (
                gga == 'B3'
                and aexx == aexx_b3
                and aggax == aggax_b3
                and aggac == aggac_b3
                and aldac == aldac_b3
            ):
                xc_functional['name'] = 'HYB_GGA_XC_B3LYP3'
            elif aexx == 1.0 and aldac == 0.0 and aggac == 0.0:
                xc_functional['name'] = 'HF_X'
            elif gga == 'PE':
                xc_functional['name'] = 'HYB_GGA_XC_PBEH'
            else:
                xc_functional['name'] = f'HYB_GGA_XC_{gga}'
            xc_functionals.append(xc_functional)
        else:
            metagga = parameters.get('METAGGA')
            if metagga:
                functionals = self.xc_functional_mapping.get(metagga, [metagga])
            else:
                functionals = self.xc_functional_mapping.get(parameters.get('GGA'), [])
            for functional in functionals:
                xc_functionals.append({'name': functional})
        return xc_functionals


class OutcarArchiveWriter(ArchiveWriter):
    def _process_pseudopotentials(  # noqa: PLR0912, PLR0915
        self,
        archive_data: Simulation,
        parser_data: dict[str, Any],
    ) -> None:
        """
        Post-process Pseudopotential instances created by mapping annotations.

        This method handles complex transformations that cannot be expressed
        in declarative mapping annotations:
        - PPCutoff subsection creation from ENMAX and ENMIN values
        - Type determination (PAW, NC-PAW, US, NC) from LPAW/LULTRA flags
        - is_norm_conserving flag based on type
        - GW optimization detection
        - XC functional subsection creation and normalization
        - Linking pseudopotentials to AtomsState

        Args:
            archive_data: The simulation archive being populated
            parser_data: Raw parsed data containing lpaw, lultra, lexch flags
        """
        # Ensure model_method exists
        if not archive_data.model_method or len(archive_data.model_method) == 0:
            return

        model_method = archive_data.model_method[0]

        # Check if any pseudopotentials were created by mapping annotations
        if not model_method.numerical_settings:
            return

        # Get pseudopotential data from parser (filter out empty dicts)
        # The OUTCAR parser extracts 4 POTCAR entries but the first 2 are just
        # short header lines without full data
        raw_pseudopotentials = [
            pp for pp in parser_data.get('pseudopotentials', []) if pp
        ]

        pp_index = (
            0  # Track actual pseudopotential index (skipping non-PP numerical_settings)
        )
        for pp in model_method.numerical_settings:
            # Skip non-pseudopotential numerical settings (if any exist)
            if not isinstance(pp, vasp.Pseudopotential):
                continue

            # Create PPCutoff subsections from ENMAX and ENMIN
            raw_pp = (
                raw_pseudopotentials[pp_index]
                if pp_index < len(raw_pseudopotentials)
                else {}
            )

            # ENMAX: recommended cutoff for standard precision
            if enmax_value := raw_pp.get('enmax'):
                enmax_cutoff = PPCutoff(
                    cutoff_kind='wavefunction',
                    cutoff_role='recommended',
                    value=enmax_value * ureg.eV,
                )
                if not pp.cutoffs:
                    pp.cutoffs = []
                pp.cutoffs.append(enmax_cutoff)

            # ENMIN: minimum recommended cutoff for fast calculations
            if enmin_value := raw_pp.get('enmin'):
                enmin_cutoff = PPCutoff(
                    cutoff_kind='wavefunction',
                    cutoff_role='recommended_min',
                    value=enmin_value * ureg.eV,
                )
                if not pp.cutoffs:
                    pp.cutoffs = []
                pp.cutoffs.append(enmin_cutoff)

            # Get lpaw, lultra, lexch from raw parser data (not stored in schema)
            # These flags are used only to derive type and is_norm_conserving
            lpaw = raw_pp.get('lpaw', False)
            lultra = raw_pp.get('lultra', False)
            lexch = raw_pp.get('lexch')
            pp_index += 1

            # Determine type from POTCAR flags and name
            # VASP uses LPAW and LULTRA flags in POTCAR (parsed from OUTCAR):
            #   - LPAW=T: PAW potential
            #     - Title contains '_GW' → NC-PAW-GW (norm-conserving)
            #     - Title contains 'nc' and 'paw' → NC-PAW (norm-conserving)
            #     - Otherwise → PAW (NOT norm-conserving)
            #   - LULTRA=T: Ultrasoft (US, NOT norm-conserving)
            #   - Both LPAW=F and LULTRA=F: Fully norm-conserving (NC)
            # VASP does NOT have an explicit "LNORMCONS" flag - absence of both
            # LPAW and LULTRA indicates standard norm-conserving pseudopotential.
            is_gw = '_GW' in pp.name if pp.name else False

            if lpaw:
                # PAW potential
                if is_gw:
                    pp.type = 'NC-PAW-GW'
                    pp.is_norm_conserving = True  # GW potentials have NC partial waves
                else:
                    # Check for NC-PAW variant (rare without GW, but possible)
                    titel_lower = pp.name.lower() if pp.name else ''
                    if 'nc' in titel_lower and 'paw' in titel_lower:
                        pp.type = 'NC-PAW'
                        pp.is_norm_conserving = True
                    else:
                        pp.type = 'PAW'
                        pp.is_norm_conserving = False
            elif lultra:
                # Ultrasoft (always Vanderbilt formalism)
                pp.type = 'US'
                pp.is_norm_conserving = False
            else:
                # Fully norm-conserving (non-PAW, non-ultrasoft)
                pp.type = 'NC'
                pp.is_norm_conserving = True

            # GW optimization detection
            if is_gw:
                pp.gw_optimized = True

            # XC functional - use lexch from raw parser data
            if lexch:
                xc = XCFunctional()
                lexch_code = lexch

                # Map VASP LEXCH codes to standard functional names
                vasp_xc_map = {
                    'PE': 'PBE',
                    'PS': 'PBEsol',
                    'CA': 'LDA',
                    '91': 'PW91',
                    'AM': 'AM05',
                    'RP': 'RPBE',
                    'PW': 'PW91',
                }

                functional_name = vasp_xc_map.get(lexch_code, lexch_code)
                xc.functional_key = functional_name
                pp.xc_functional = xc

        # Link pseudopotentials to AtomsState
        # VASP lists pseudopotentials in POSCAR species order
        if archive_data.model_system and len(archive_data.model_system) > 0:
            model_system = archive_data.model_system[0]
            if (
                hasattr(model_system, 'particle_states')
                and model_system.particle_states
            ):
                # Get only the Pseudopotential objects
                pseudopotentials = [
                    ns
                    for ns in model_method.numerical_settings
                    if isinstance(ns, vasp.Pseudopotential)
                ]

                # Link each AtomsState to its corresponding Pseudopotential
                for atoms_state, pp in zip(
                    model_system.particle_states, pseudopotentials
                ):
                    atoms_state.pseudopotential = pp

    def write_to_archive(self) -> None:
        # set up archive parser
        archive_data_parser = VASPMetainfoParser()
        archive_data = Simulation()
        archive_data_parser.data_object = archive_data
        archive_data_parser.annotation_key = vasp.OUTCAR_KEY

        # set up outcar parser
        source_parser = OutcarParser()
        source_parser.text_parser = OutcarTextParser()
        source_parser.filepath = self.mainfile

        # convert
        source_parser.convert(archive_data_parser)

        # Post-process pseudopotentials created by mapping annotations
        # This handles complex transformations (type, XC functional, linking)
        # Pass parser data to access lpaw, lultra, lexch flags (not stored in schema)
        self._process_pseudopotentials(archive_data, source_parser.data)

        # assign simulation section to archive data
        self.archive.data = archive_data_parser.data_object

        # close file handles
        archive_data_parser.close()
        source_parser.close()
