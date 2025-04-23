import os
from typing import Any

import numpy as np
from ase.io import read
from nomad.datamodel import EntryArchive
from nomad.parsing.parser import MatchingParser
from nomad.units import ureg
from nomad.utils import get_logger
from nomad_file_parser import ArchiveWriter
from nomad_file_parser.mapping_parser import MetainfoParser, TextParser
from nomad_simulations.schema_packages.general import Program, Simulation
from structlog.stdlib import BoundLogger

from nomad_simulation_parsers.parsers.utils.general import (
    calculate_band_gap_from_occupations,
)
from nomad_simulation_parsers.schema_packages import octopus

from .file_parser import EigenvalueParser, InfoParser, InpParser, LogParser, OutParser

LOGGER = get_logger(__name__)


class OctopusMainfileParser(TextParser):
    log_parser = LogParser()
    inp_parser = InpParser()
    # List was generated from https://github.com/qsnake/octopus/blob/master/share/varinfo
    # In parser.log, XCFunctional is an integer given as a sum of exchange
    # and correlation. In stdout.txt, XCFunctional can be read as a string
    # but it is not always consistent with the mapping that we have below (see for
    # example Exchange which is Slater Exchange in stdout.txt). Therefore, in
    # resolving the nomad xc_functional name, we first try to get it from the string
    # output of stdout.txt and if it does not work, we get it from the integer value
    # in parser.log.
    # TODO resolve solely from the integer value
    _xc_functionals = {
        'Exchange': ('lda_x', 1),
        'Slater exchange': ('lda_x', 1),
        'Wigner parametrization': ('lda_c_wigner', 2000),
        'Random Phase Approximation': ('lda_c_rpa', 3000),
        'Hedin & Lundqvist': ('lda_c_hl', 4000),
        'Gunnarson & Lundqvist': ('lda_c_gl', 5000),
        'Slater Xalpha': ('lda_c_xalpha', 6000),
        'Vosko, Wilk, & Nussair': ('lda_c_vwn', 7000),
        'Vosko, Wilk, & Nussair (RPA)': ('lda_c_vwn_rpa', 8000),
        'Perdew & Zunger': ('lda_c_pz', 9000),
        'Perdew & Zunger (Modified)': ('lda_c_pz_mod', 10000),
        'Ortiz & Ballone (PZ)': ('lda_c_ob_pz', 11000),
        'Perdew & Wang': ('lda_c_pw', 12000),
        'Perdew & Wang (Modified)': ('lda_c_pw_mod', 13000),
        'Ortiz & Ballone (PW)': ('lda_c_ob_pw', 14000),
        'Attacalite et al': ('lda_c_2d_amgb', 15000),
        'Pittalis, Rasanen & Marques correlation in 2D': ('lda_c_2d_prm', 16000),
        'von Barth & Hedin': ('lda_c_vbh', 17000),
        'Casula, Sorella, and Senatore 1D correlation': ('lda_c_1d_csc', 18000),
        'Exchange in 2D': ('lda_x_2d', 19),
        'Teter 93 parametrization': ('lda_xc_teter93', 20000),
        'Exchange in 1D': ('lda_x_1d', 21),
        'Modified LSD (version 1) of Proynov and Salahub': ('lda_c_ml1', 22000),
        'Modified LSD (version 2) of Proynov and Salahub': ('lda_c_ml2', 23000),
        'Gombas parametrization': ('lda_c_gombas', 24000),
        'Thomas-Fermi kinetic energy functional': ('lda_k_tf', 50),
        'Lee and Parr Gaussian ansatz': ('lda_k_lp', 51),
        'Perdew, Burke & Ernzerhof exchange': ('gga_x_pbe', 101),
        'Perdew, Burke & Ernzerhof exchange (revised)': ('gga_x_pbe_r', 102),
        'Becke 86 Xalfa,beta,gamma': ('gga_x_b86', 103),
        'Herman et al original GGA': ('gga_x_herman', 104),
        'Becke 86 Xalfa,beta,gamma (with mod. grad. correction)': (
            'gga_x_b86_mgc',
            105,
        ),
        'Becke 88': ('gga_x_b88', 106),
        'Gill 96': ('gga_x_g96', 107),
        'Perdew & Wang 86': ('gga_x_pw86', 108),
        'Perdew & Wang 91': ('gga_c_pw91', 134000),
        'Handy & Cohen OPTX 01': ('gga_x_optx', 110),
        'dePristo & Kress 87 (version R1)': ('gga_x_dk87_r1', 111),
        'dePristo & Kress 87 (version R2)': ('gga_x_dk87_r2', 112),
        'Lacks & Gordon 93': ('gga_x_lg93', 113),
        'Filatov & Thiel 97 (version A)': ('gga_x_ft97_a', 114),
        'Filatov & Thiel 97 (version B)': ('gga_x_ft97_b', 115),
        'Perdew, Burke & Ernzerhof exchange (solids)': ('gga_x_pbe_sol', 116),
        'Hammer, Hansen & Norskov (PBE-like)': ('gga_x_rpbe', 117),
        'Wu & Cohen': ('gga_x_wc', 118),
        'Modified form of PW91 by Adamo & Barone': ('gga_x_mpw91', 119),
        'Armiento & Mattsson 05 exchange': ('gga_x_am05', 120),
        'Madsen (PBE-like)': ('gga_x_pbea', 121),
        'Adamo & Barone modification to PBE': ('gga_x_mpbe', 122),
        'xPBE reparametrization by Xu & Goddard': ('gga_c_xpbe', 136000),
        'Becke 86 MGC for 2D systems': ('gga_x_2d_b86_mgc', 124),
        'Bayesian best fit for the enhancement factor': ('gga_x_bayesian', 125),
        'JSJR reparametrization by Pedroza, Silva & Capelle': (
            'gga_x_pbe_jsjr',
            126,
        ),
        'Becke 88 in 2D': ('gga_x_2d_b88', 127),
        'Perdew, Burke & Ernzerhof exchange in 2D': ('gga_x_2d_pbe', 129),
        'Perdew, Burke & Ernzerhof correlation': ('gga_c_pbe', 130000),
        'Lee, Yang & Parr': ('gga_c_lyp', 131000),
        'Perdew 86': ('gga_c_p86', 132000),
        'Perdew, Burke & Ernzerhof correlation SOL': ('gga_c_pbe_sol', 133000),
        'Armiento & Mattsson 05 correlation': ('gga_c_am05', 135000),
        'Langreth and Mehl correlation': ('gga_c_lm', 137000),
        'JRGX reparametrization by Pedroza, Silva & Capelle': (
            'gga_c_pbe_jrgx',
            138000,
        ),
        'Becke 88 reoptimized to be used with vdW functional of Dion et al': (
            'gga_x_optb88_vdw',
            139,
        ),
        'PBE reparametrization for vdW': ('gga_x_optpbe_vdw', 141),
        'Regularized PBE': ('gga_c_rge2', 143000),
        'refitted Perdew & Wang 86': ('gga_x_rpw86', 144),
        'Keal and Tozer version 1': ('gga_x_kt1', 145),
        'Keal and Tozer version 2': ('gga_xc_kt2', 146000),
        'Wilson & Levy': ('gga_c_wl', 147000),
        'Wilson & Ivanov': ('gga_c_wi', 148000),
        'van Leeuwen & Baerends': ('gga_x_lb', 160),
        'HCTH functional fitted to 93 molecules': ('gga_xc_hcth_93', 161000),
        'HCTH functional fitted to 120 molecules': ('gga_xc_hcth_120', 162000),
        'HCTH functional fitted to 147 molecules': ('gga_xc_hcth_407', 164000),
        'Empirical functionals from Adamson, Gill, and Pople': (
            'gga_xc_edf1',
            165000,
        ),
        'XLYP functional': ('gga_xc_xlyp', 166000),
        'Becke 97': ('hyb_gga_xc_b97', 407000),
        'Becke 97-1': ('hyb_gga_xc_b97_1', 408000),
        'Becke 97-2': ('hyb_gga_xc_b97_2', 410000),
        'Grimme functional to be used with C6 vdW term': ('gga_xc_b97_d', 170000),
        'Boese-Martin for Kinetics': ('hyb_gga_xc_b97_k', 413000),
        'Becke 97-3': ('hyb_gga_xc_b97_3', 414000),
        'Functionals fitted for water': ('gga_xc_pbelyp1w', 175000),
        'Schmider-Becke 98 parameterization 1a': ('hyb_gga_xc_sb98_1a', 420000),
        'Schmider-Becke 98 parameterization 1b': ('hyb_gga_xc_sb98_1b', 421000),
        'Schmider-Becke 98 parameterization 1c': ('hyb_gga_xc_sb98_1c', 422000),
        'Schmider-Becke 98 parameterization 2a': ('hyb_gga_xc_sb98_2a', 423000),
        'Schmider-Becke 98 parameterization 2b': ('hyb_gga_xc_sb98_2b', 424000),
        'Schmider-Becke 98 parameterization 2c': ('hyb_gga_xc_sb98_2c', 425000),
        'van Leeuwen & Baerends modified': ('gga_x_lbm', 182),
        'von Weiszaecker correction to Thomas-Fermi': ('gga_k_vw', 500),
        'Second-order gradient expansion (l = 1/9)': ('gga_k_ge2', 501),
        'TF-lambda-vW form by Golden (l = 13/45)': ('gga_k_golden', 502),
        'TF-lambda-vW form by Yonei and Tomishima (l = 1/5)': ('gga_k_yt65', 503),
        'TF-lambda-vW form by Baltin (l = 5/9)': ('gga_k_baltin', 504),
        'TF-lambda-vW form by Lieb (l = 0.185909191)': ('gga_k_lieb', 505),
        'gamma-TFvW form by Acharya et al [g = 1 - 1.412/N^(1/3)]': (
            'gga_k_absr1',
            506,
        ),
        'gamma-TFvW form by Acharya et al [g = 1 - 1.332/N^(1/3)]': (
            'gga_k_absr2',
            507,
        ),
        'gamma-TFvW form by Gázquez and Robles': ('gga_k_gr', 508),
        'gamma-TFvW form by Ludeña': ('gga_k_ludena', 509),
        'gamma-TFvW form by Ghosh and Parr': ('gga_k_gp85', 510),
        'Pearson': ('gga_k_pearson', 511),
        'Ou-Yang and Levy v.1': ('gga_k_ol1', 512),
        'Ou-Yang and Levy v.2': ('gga_k_ol2', 513),
        'Fuentealba & Reyes (B88 version)': ('gga_k_fr_b88', 514),
        'Fuentealba & Reyes (PW86 version)': ('gga_k_fr_pw86', 515),
        'The original hybrid proposed by Becke': ('hyb_gga_xc_b3pw91', 401000),
        'The (in)famous B3LYP': ('hyb_gga_xc_b3lyp', 402000),
        'Perdew 86 hybrid similar to B3PW91': ('hyb_gga_xc_b3p86', 403000),
        'hybrid using the optx functional': ('hyb_gga_xc_o3lyp', 404000),
        'mixture of mPW91 and PW91 optimized for kinetics': (
            'hyb_gga_xc_mpw1k',
            405000,
        ),
        'aka PBE0 or PBE1PBE': ('hyb_gga_xc_pbeh', 406000),
        'maybe the best hybrid': ('hyb_gga_xc_x3lyp', 411000),
        'Becke 1-parameter mixture of WC and PBE': ('hyb_gga_xc_b1wc', 412000),
        'mixture with the mPW functional': ('hyb_gga_xc_mpw3pw', 415000),
        'Becke 1-parameter mixture of B88 and LYP': ('hyb_gga_xc_b1lyp', 416000),
        'Becke 1-parameter mixture of B88 and PW91': ('hyb_gga_xc_b1pw91', 417000),
        'Becke 1-parameter mixture of mPW91 and PW91': (
            'hyb_gga_xc_mpw1pw',
            418000,
        ),
        'mixture of mPW and LYP': ('hyb_gga_xc_mpw3lyp', 419000),
        'Local tau approximation of Ernzerhof & Scuseria': ('mgga_x_lta', 201),
        'Perdew, Tao, Staroverov & Scuseria exchange': ('mgga_x_tpss', 202),
        'Zhao, Truhlar exchange': ('mgga_x_m06l', 203),
        'GVT4 from Van Voorhis and Scuseria (exchange part)': ('mgga_x_gvt4', 204),
        'tau-HCTH from Boese and Handy': ('mgga_x_tau_hcth', 205),
        'Becke-Roussel 89': ('mgga_x_br89', 206),
        'Becke & Johnson correction to Becke-Roussel 89': ('mgga_x_bj06', 207),
        'Tran & Blaha correction to Becke & Johnson': ('mgga_x_tb09', 208),
        'Rasanen, Pittalis, and Proetto correction to Becke & Johnson': (
            'mgga_x_rpp09',
            209,
        ),
        'Pittalis, Rasanen, Helbig, Gross Exchange Functional': (
            'mgga_x_2d_prhg07',
            210,
        ),
        'PRGH07 with PRP10 correction': ('mgga_x_2d_prhg07_prp10', 211),
        'Perdew, Tao, Staroverov & Scuseria correlation': ('mgga_c_tpss', 231000),
        'VSxc from Van Voorhis and Scuseria (correlation part)': (
            'mgga_c_vsxc',
            232000,
        ),
        'Orestes, Marcasso & Capelle': ('lca_omc', 301),
        'Lee, Colwell & Handy': ('lca_lch', 302),
        'OEP: Exact exchange': ('oep_x', 901),
    }
    _units_mapping = dict(
        ev=ureg.eV, hartree=ureg.hartree, angstrom=ureg.angstrom, bohr=ureg.bohr
    )
    _info = None
    _initial_system = None

    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER

    def init_parser(self):
        self._info = None
        self._initial_system = None
        self._maindir = os.path.dirname(self.filepath)
        self.log_parser.mainfile = os.path.join(self._maindir, 'exec/parser.log')
        self.inp_parser.mainfile = os.path.join(self._maindir, 'inp')

    @property
    def info(self):
        if self._info is not None:
            return self._info

        self._info = dict()
        self._info.update(self.inp_parser.info)
        # the variables read from log are more reliable as these are converted
        # to proper octopus values, i.e. we do not worry about compatibility issue
        self._info.update(self.log_parser.info)
        # add units for energy, length
        # TODO not sure what defaults are for units and expected values for units(input)
        # and their meaning
        energyunit = 'hartree'
        lengthunit = 'bohr'
        units = self._info.get('Units', self._info.get('UnitsInput', 0))
        if isinstance(units, str):
            units = units.lower()
            energyunit = 'eV' if 'ev' in units else energyunit
            lengthunit = 'angstrom' if 'angs' in units else lengthunit
        elif isinstance(units, int):
            energyunit = 'eV' if units == 1 else energyunit
            lengthunit = 'angstrom' if units == 1 else lengthunit
        elif units is None:
            cell = self.out_parser.get('grid', {}).get('cell')
            if cell is not None:
                lengthunit = cell.units
                energyunit = 'eV' if lengthunit == 'angstrom' else energyunit
        self._info.update(dict(energyunit=energyunit, lengthunit=lengthunit))
        return self._info

    def get_header(self, header: list[list[str]]) -> dict[str, str]:
        return {key: val for key, val in header.get('options', [])}

    def get_xc_functionals(self, theory_level: dict[str, Any]) -> list[str]:
        xc_functionals = []
        xc_int_sum = 0
        xc_functional = self.info.get('XCFunctional', 0)
        if isinstance(xc_functional, str):
            xc_functionals = [xc.strip() for xc in xc_functional.upper().split('+')]
        else:
            # we try to resolve from output string
            for key in ['exchange', 'correlation']:
                val = theory_level.get(key, None)
                name, number = self._xc_functionals.get(val, (None, 0))
                if name is not None:
                    xc_functionals.append(name.upper())
                    xc_int_sum += number
            # get it from log
            diff = xc_functional - xc_int_sum
            if diff > 0:
                names_numbers = list(zip(*self._xc_functionals.values()))
                while diff > 0:
                    index = np.argmin(abs(np.array(names_numbers[1]) - diff))
                    xc_functionals.append(names_numbers[0][index].upper())
                    diff -= names_numbers[1][index]
                if diff < 0:
                    xc_functionals = []
                    self.logger.error(
                        'Error resolving xc functional',
                        data=dict(XCfunctional=xc_functional),
                    )
        return xc_functionals

    @property
    def initial_system(self) -> dict[str, Any]:  # noqa: PLR0912
        if self._initial_system is not None:
            return self._initial_system

        self._initial_system = {}
        atoms = None  # Initialize atoms to None for all paths
        symbols, coordinates = self.log_parser.get_coordinates()
        if len(coordinates) == 0:
            # get if from inp
            symbols, coordinates = self.inp_parser.get_coordinates()
        if len(coordinates) == 0:
            # try to read from file
            file_types = {
                'PDBCoordinates': ['proteindatabank'],
                'XYZCoordinates': ['extxyz', 'xyz'],
                'XSFCoordinates': ['xsf'],
            }
            filenames = []
            for ftype, fformats in file_types.items():
                filename = self.info.get(ftype)
                if filename is None:
                    continue
                filenames.append(filename)
                filepath = os.path.join(self._maindir, filename)
                if not os.path.isfile(filepath):
                    continue
                for fformat in fformats:
                    try:
                        atoms = read(filepath, format=fformat)
                    except Exception:
                        continue
                    if atoms is not None:
                        break
                if atoms is not None:
                    symbols = atoms.get_chemical_symbols()
                    coordinates = atoms.get_positions()
                    break
            if atoms is None and filenames:
                self.logger.error(
                    'Error reading coordinates files', data=dict(filenames=filenames)
                )

        if len(coordinates) == 0:
            self.logger.error('Error parsing atom positions and labels.')
            return

        cell = self.data.get('grid', {}).get('cell')
        if cell is not None:
            npbc = self.data.get('grid', {}).get('npbc', 3)
            self._initial_system['pbc'] = [True for _ in range(npbc)]
            self._initial_system['lattice_vectors'] = cell
        elif atoms is not None:
            ase_cell = atoms.get_cell()
            if ase_cell is not None and np.asarray(ase_cell).size > 0:
                self._initial_system['lattice_vectors'] = (
                    np.asarray(ase_cell) * self._units_mapping['angstrom']
                )
            if atoms.pbc is not None:
                self._initial_system['pbc'] = [
                    bool(val) for val in np.asarray(atoms.pbc)
                ]

        if self.info.get('ReducedCoordinates', None) is not None and cell is not None:
            coordinates = np.dot(coordinates, cell.magnitude)
            units = cell.units
        elif self.info.get('Coordinates', None) is not None:
            units = self.info.get('lengthunit')
        else:
            # read from ase atoms (in angstroms)
            units = 'angstrom'

        self._initial_system['positions'] = coordinates * self._units_mapping.get(
            str(units).lower()
        )
        self._initial_system['labels'] = symbols
        return self._initial_system

    def get_systems(
        self, source: dict[str, Any] | list[dict[str, Any]] | None = None
    ) -> list[dict[str, Any]]:
        # Always include initial system, even for static calculations
        initial = self.initial_system
        if not initial:
            return []

        systems = [initial]

        # Extract minimization data from source (dict or list)
        minimization = None
        if isinstance(source, dict):
            minimization = source.get('minimization')
        elif isinstance(source, list):
            minimization = source

        # Add geometry optimization steps if present
        unreadable_paths = []
        for mini in minimization or []:
            number = mini.get('numbr')
            if number is None:
                continue
            path = os.path.join(self._maindir, f'geom/go.{number:04d}.xyz')
            if not os.path.isfile(path):
                unreadable_paths.append(path)
                continue
            try:
                atoms = read(path, format='xyz')
            except Exception:
                unreadable_paths.append(path)
                continue
            systems.append(
                dict(
                    positions=atoms.get_positions() * ureg.angstrom,
                    labels=atoms.get_chemical_symbols(),
                    pbc=systems[0].get('pbc'),
                )
            )
        if unreadable_paths:
            self.logger.warning(
                'Could not read geometry optimization steps',
                data=dict(paths=unreadable_paths),
            )
        return systems

    def get_outputs(self, sources: list[dict[str]]) -> list[dict[str, Any]]:
        outputs = []
        # for source in [s for ss in sources for s in ss]:
        for source in sources:
            energy = source.get('energy', source.get('energy_total'))
            if energy is None:
                continue
            outputs.append(
                dict(
                    energy=energy
                    * self._units_mapping.get(
                        self.info.get('energyunit', 'hartree').lower()
                    )
                )
            )
        return outputs


class OctopusMetainfoParser(MetainfoParser):
    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER


class OctopusEigenvalueParser(TextParser):
    unit = ureg.hartree

    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER

    def get_eigenvalues(self, source: list[np.ndarray]) -> list[dict[str, Any]]:
        eigen_section = self.data.get('eigenvalues') if hasattr(self, 'data') else None
        fermi = (
            eigen_section.get('fermi_energy')
            if isinstance(eigen_section, dict)
            else None
        )
        reference_energy = fermi.to(self.unit) if fermi is not None else None
        kpts, eigs, occs = list(zip(*[e for e in source if e is not None]))
        eigs = np.transpose(eigs, axes=(2, 0, 1))
        occs = np.transpose(occs, axes=(2, 0, 1))
        return [
            dict(
                eigenvalues=eig * self.unit,
                occupations=occs[n],
                kpoints=kpts,
                highest_occupied=reference_energy,
            )
            for n, eig in enumerate(eigs)
        ]

    def get_band_structures(self, source: list[np.ndarray]) -> list[dict[str, Any]]:
        return [
            dict(
                value=eigenvalue_data['eigenvalues'],
                highest_occupied=eigenvalue_data.get('highest_occupied'),
            )
            for eigenvalue_data in self.get_eigenvalues(source)
            if eigenvalue_data.get('kpoints') is not None
        ]

    def get_band_gaps(self, source: list[np.ndarray]) -> list[dict[str, Any]]:
        """Calculate band gaps from eigenvalues using common utility."""
        band_gaps = []
        for spin_channel, eigenvalue_data in enumerate(self.get_eigenvalues(source)):
            energies = eigenvalue_data.get('eigenvalues')
            occupations = eigenvalue_data.get('occupations')

            # Use common utility for band gap calculation (handles units automatically)
            gap_result = calculate_band_gap_from_occupations(
                energies, occupations, spin_channel=spin_channel
            )
            if gap_result is not None:
                band_gaps.append(gap_result)

        return band_gaps


class OctopusInfoParser(OctopusEigenvalueParser):
    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER


class OctopusArchiveWriter(ArchiveWriter):
    mainfile_parser = OctopusMainfileParser(text_parser=OutParser())
    archive_parser = OctopusMetainfoParser()
    info_parser = OctopusInfoParser(text_parser=InfoParser())
    eigenvalues_parser = OctopusEigenvalueParser(text_parser=EigenvalueParser())

    def write_to_archive(self) -> None:
        self.mainfile_parser.filepath = self.mainfile
        # initialize auxiliary file parsers
        self.mainfile_parser.init_parser()
        self.archive.data = Simulation(program=Program(name='Octopus'))

        self.archive_parser.data_object = self.archive.data
        self.archive_parser.annotation_key = octopus.OUT_KEY

        self.mainfile_parser.convert(self.archive_parser)

        maindir = os.path.dirname(self.mainfile)

        # read quantities from info file
        self.info_parser.filepath = os.path.join(maindir, 'static/info')
        # pass the energy unit
        energy_unit = self.mainfile_parser._units_mapping.get(
            self.mainfile_parser.info.get('energyunit')
        )
        self.info_parser.unit = energy_unit
        self.archive_parser.annotation_key = octopus.INFO_KEY
        self.info_parser.convert(self.archive_parser, update_mode='merge@-1')

        # read eigenvalues from eigenvalues file
        self.eigenvalues_parser.filepath = os.path.join(maindir, 'static/eigenvalues')
        # pass the energy unit
        self.eigenvalues_parser.unit = energy_unit
        self.archive_parser.annotation_key = octopus.EIGENVALUES_KEY
        self.eigenvalues_parser.convert(self.archive_parser, update_mode='merge@-1')

        self.eigenvalues_parser.close()
        self.archive_parser.close()
        self.mainfile_parser.close()
        self.info_parser.close()


class OctopusParser(MatchingParser):
    archive_writer = OctopusArchiveWriter()

    """
    Main parser interface to NOMAD.
    """

    def parse(
        self,
        mainfile: str,
        archive: EntryArchive,
        logger: BoundLogger,
        child_archives: dict[str, EntryArchive] = {},
    ) -> None:
        self.archive_writer.write(mainfile, archive, logger, child_archives)

        # run the old parser
        # TODO remove
        from electronicparsers.octopus.parser import OctopusParser

        OctopusParser().parse(mainfile, archive, logger)