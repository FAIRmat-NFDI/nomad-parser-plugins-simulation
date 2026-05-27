from typing import Any

import ase
import numpy as np
import pint
from ase.io.ulm import Reader
from nomad_file_parser import DataTextParser, FileParser, TarParser, XMLParser
from nomad.units import ureg


class GPWTarParser(TarParser):
    def __init__(self) -> None:
        super().__init__()
        self._version_map = {6: '1.1.0', 5: '0.11.0', 3: '0.10.0'}
        self._info_map = {
            'energy_total': 'Epot',
            'energy_XC': 'Exc',
            'electronic_kinetic_energy': 'Ekin',
            'energy_correction_entropy': 'S',
            'atom_forces_free': 'CartesianForce',
            'atom_positions': 'CartesianPositions',
            'occupation': 'OccupationNumbers',
            'kpoints': 'IBZKPoints',
        }

    def reset(self) -> None:
        super().reset()
        self._info = None

    @property
    def info(self) -> dict[str, Any]:
        if self._info is None:
            self._info = dict()
            xml_file = self.get('info.xml', None)
            if xml_file is None:
                return
            self.xml_file = xml_file
            xml_parser = XMLParser(xml_file, logger=self.logger)

            def convert(val):
                if isinstance(val, list):
                    return [convert(v) for v in val]
                try:
                    if val in ['True', 'False']:
                        return val == 'True'
                    else:
                        val = float(val)
                        if val % 1.0 == 0.0:
                            val = int(val)
                        return val
                except Exception:
                    return val

            # parameters
            self._info['parameter'] = {
                'lengthunit': 'angstrom',
                'energyunit': 'eV',
                'timeunit': 'femtosecond',
            }
            self._info['parameter'].update(
                {
                    p['name'].lower(): convert(p['value'])
                    for p in xml_parser.get('parameter/', [])
                }
            )

            # array shapes, types, dimensions
            self._info['array'] = dict()
            dimension = dict()
            for arr in xml_parser.root.findall('./array'):
                name = arr.attrib.get('name', None)
                dtype = arr.attrib.get('type', None)
                if name is None or dtype is None:
                    continue
                shape = []
                for dim in arr.findall('./dimension'):
                    length = int(dim.attrib.get('length', 0))
                    shape.append(length)
                    dimension[dim.attrib.get('name')] = length
                self._info['array'][name.lower()] = dict(dtype=dtype, shape=shape)
            self._info['array_dimension'] = dimension

            self._info['bytesswap'] = (
                xml_parser.root.attrib['endianness'] == 'little'
            ) != np.little_endian

            xml_parser.close()

        return self._info

    def get_parameter(self, key: str, unit=None) -> Any:
        key = self._info_map.get(key, key)
        return self.info['parameter'].get(key.lower(), None)

    def get_array(self, key: str, unit=None) -> np.ndarray:
        key = self._info_map.get(key, key)
        file_object = self.get(key)
        if file_object is None:
            return

        key = key.lower()
        shape = self.info['array'].get(key, {}).get('shape', None)
        dtype = self.info['array'].get(key, {}).get('dtype', None)
        dtype = np.dtype({'int': 'int32'}.get(dtype, dtype))
        size = np.prod(shape) * dtype.itemsize

        file_object.seek(0)
        parser = DataTextParser(
            mainfile_contents=file_object.read(size), logger=self.logger, dtype=dtype
        )
        if parser.data is None:
            return

        array = parser.data
        if self._info['bytesswap']:
            array = array.byteswap()
        if dtype == np.int32:
            array = np.asarray(array, int)
        array.shape = shape
        return array

    def get_array_dimension(self, key: str) -> list[int]:
        if key == 'ngpts':
            val = [self.get_array_dimension(f'ngpts{n}') for n in ['x', 'y', 'z']]
        else:
            val = self.info['array_dimension'].get(key)
        return val

    def get_program_version(self) -> str:
        return self._version_map.get(self.get_parameter('version'), '0.9.0')

    def get_smearing_width(self) -> float:
        return self.get_parameter('fermiwidth')


class GPW2FileParser(FileParser):
    def reset(self):
        super().reset()
        self._info = None

    @property
    def ulm(self):
        if self._file_handler is None:
            try:
                self._file_handler = Reader(self.mainfile)
            except Exception:
                pass
        return self._file_handler

    @property
    def info(self):
        if self._info is None:
            self._info = dict()
            self._info['parameter'] = {
                'mode': 'fd',
                'xc': 'LDA',
                'occupations': None,
                'poissonsolver': None,
                'h': None,
                'gpts': None,
                'kpts': [(0.0, 0.0, 0.0)],
                'nbands': None,
                'charge': 0,
                'setups': {},
                'basis': {},
                'spinpol': None,
                'fixdensity': False,
                'filter': None,
                'mixer': None,
                'eigensolver': None,
                'background_charge': None,
                'external': None,
                'random': False,
                'hund': False,
                'maxiter': 333,
                'idiotproof': True,
                'symmetry': {
                    'point_group': True,
                    'time_reversal': True,
                    'symmorphic': True,
                    'tolerance': 1e-7,
                },
                'convergence': {
                    'energy': 0.0005,
                    'density': 1.0e-4,
                    'eigenstates': 4.0e-8,
                    'bands': 'occupied',
                    'forces': np.inf,
                },
                'dtype': None,
                'width': None,
                'verbose': 0,
                'lengthunit': 'angstrom',
                'energyunit': 'eV',
                'timeunit': 'femtosecond',
            }
            if self.ulm is not None:
                self._info['parameter'].update(self.ulm.parameters.asdict())

        self._info.update(
            {
                'planewavecutoff': self._info['parameter']
                .get('mode', {})
                .get('ecut', None),
                'basisset': self._info['parameter'].get('basis'),
                'energyerror': self._info['parameter']
                .get('convergence', {})
                .get('energy', None),
                'xcfunctional': self._info['parameter'].get('xc'),
            }
        )

        for key, f in {
            'energy_total': lambda: self.ulm.hamiltonian.e_total_extrapolated,
            'energy_free': lambda: self.ulm.hamiltonian.e_total_free,
            'energy_XC': lambda: self.ulm.hamiltonian.e_xc,
            'electronic_kinetic_energy': lambda: self.ulm.hamiltonian.e_kinetic,
            'energy_correction_entropy': lambda: self.ulm.hamiltonian.e_entropy,
            'fermilevel': lambda: self.ulm.occupations.fermilevel,
            'split': lambda: self.ulm.occupations.split,
            'converged': lambda: self.ulm.scf.converged,
        }.items():
            try:
                self._info[key] = f()
            except Exception:
                pass

        return self._info

    def get_parameter(self, key: str) -> Any:
        return self.info.get(key, self.info['parameter'].get(key))

    def get_array(self, key: str) -> np.ndarray:
        if self.ulm is None:
            return
        values = {
            'unitcell': self.ulm.atoms.cell,
            'atomicnumbers': self.ulm.atoms.numbers,
            'atom_positions': self.ulm.atoms.positions,
            'boundaryconditions': self.ulm.atoms.pbc,
            'momenta': self.ulm.atoms.momenta,
            'atom_forces_free_raw': self.ulm.results.forces,
            'magneticmoments': self.ulm.results.magmoms,
            'eigenvalues': self.ulm.wave_functions.eigenvalues,
            'occupation': self.ulm.wave_functions.occupations,
            # TODO no koints data in ulm?
            'kpoints': lambda: self.ulm.IBZKPoints,
            'density': lambda: self.ulm.density.density,
            'potential_effective': lambda: self.ulm.hamiltonian.potential,
            'band_paths': self.ulm.wave_functions.band_paths.asdict,
        }
        try:
            if key in values:
                val = values.get(key)
            else:
                val = self.ulm.asdict().get(key, None)
        except Exception:
            val = None
        return val

    def get_array_dimension(self, key: str) -> list[int]:
        if key == 'ngpts':
            val = self.ulm.density.density.shape
        else:
            val = self.ulm.asdict().get(key, None)
        return val

    def get_program_version(self) -> str:
        return self.ulm.gpaw_version

    def get_smearing_width(self) -> float:
        if self.get_parameter('occupations') is None:
            return 0.0 if tuple(self.get_parameter('kpts')) == (1, 1, 1) else 0.1
        else:
            return self.get_parameter('occupations').get('width')

    def parse(self, key=None) -> None:
        pass


class GPWFileParser(FileParser):
    _units_map = {
        'ev': ureg.eV,
        'hartree': ureg.hartree,
        'angstrom': ureg.angstrom,
        'bohr': ureg.bohr,
        'femtosecond': ureg.fs,
    }
    _xc_map = {
        'LDA': ['LDA_X', 'LDA_C_PW'],
        'PW91': ['GGA_X_PW91', 'GGA_C_PW91'],
        'PBE': ['GGA_X_PBE', 'GGA_C_PBE'],
        'PBEsol': ['GGA_X_PBE_SOL', 'GGA_C_PBE_SOL'],
        'revPBE': ['GGA_X_PBE_R', 'GGA_C_PBE'],
        'RPBE': ['GGA_X_RPBE', 'GGA_C_PBE'],
        'BLYP': ['GGA_X_B88', 'GGA_C_LYP'],
        'HCTH407': ['GGA_XC_HCTH_407'],
        'WC': ['GGA_X_WC', 'GGA_C_PBE'],
        'AM05': ['GGA_X_AM05', 'GGA_C_AM05'],
        'M06-L': ['MGGA_X_M06_L', 'MGGA_C_M06_L'],
        'TPSS': ['MGGA_X_TPSS', 'MGGA_C_TPSS'],
        'revTPSS': ['MGGA_X_REVTPSS', 'MGGA_C_REVTPSS'],
        'mBEEF': ['MGGA_X_MBEEF', 'GGA_C_PBE_SOL'],
    }
    parser = GPWTarParser()

    def apply_unit(self, val: np.ndarray | float, unit: str) -> pint.Quantity:
        if val is None:
            return

        p_unit = self.parser.info['parameter'].get(unit, '').lower()
        unit = self._units_map.get(p_unit, p_unit) if p_unit else unit
        return val * unit

    def get_mode(self) -> str:
        mode = self.parser.get_parameter('mode')
        if isinstance(mode, dict):
            mode = mode.get('name')
        return mode

    def parse(self, key: str = None):
        self.parser = GPWTarParser()
        self.parser.mainfile = self.mainfile
        if self.parser.mainfile_obj is None:
            self.parser = GPW2FileParser()
            self.parser.mainfile = self.mainfile

        self._results = {'program_version': self.parser.get_program_version()}
        self._results['unitcell'] = self.apply_unit(
            self.parser.get_array('unitcell'), 'lengthunit'
        )
        self._results['atom_positions'] = self.apply_unit(
            self.parser.get_array('atom_positions'), 'lengthunit'
        )
        self._results['labels'] = [
            ase.data.chemical_symbols[z] for z in self.parser.get_array('atomicnumbers')
        ]

        pbc = np.ones(3, bool) if self.get_mode() == 'pw' else np.zeros(3, bool)
        if self.parser.get_array('boundary_conditions') is not None:
            bc = np.array(self.parser.get_array('boundary_conditions'), bool)
            bc.shape = [bc.size]
            pbc[: bc.size] = bc
        self._results['boundary_conditions'] = pbc
        xc_functional = self.parser.get_parameter('xcfunctional')
        self._results['xcfunctional'] = [
            xc for xc in self._xc_map.get(xc_functional, [xc_functional])
        ]
