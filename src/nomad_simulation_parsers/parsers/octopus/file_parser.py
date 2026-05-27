import re
from typing import Any

import numpy as np
import pint
from ase import units as ase_units
from nomad_file_parser.text_parser import Quantity, TextParser
from nomad.units import ureg

RE_FLOAT = r'[\d\.\-\+Ee]+'


class EigenvalueParser(TextParser):
    def init_quantities(self) -> None:
        def str_to_eigenvalues(
            val_in: str,
        ) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
            val = [v.split() for v in val_in.strip().split('\n')]
            kpoint = np.array(
                [0.0, 0.0, 0.0]
                if val[0][0] == '#st'
                else [float(v.rstrip(',')) for v in val[0]]
            )
            if len(val) == 1:
                return

            eigenvalues = np.array([[v[0], v[2], v[3]] for v in val[1:]], dtype=float)
            eigenvalues = np.transpose(eigenvalues)
            nspin = 2 if eigenvalues[0][1] == 1.0 else 1
            nbands = int(max(eigenvalues[0]))
            eigenvalues = eigenvalues[1:]
            eigenvalues.shape = (2, nbands, nspin)
            return kpoint, eigenvalues[0], eigenvalues[1]

        def str_to_fermi_energy(val_in: str) -> pint.Quantity:
            val = val_in.split()
            unit = ureg.eV if val[1].startswith('e') else ureg.hartree
            return float(val[0]) * unit

        self._quantities = [
            Quantity(
                'eigenvalues',
                r'(Eigenvalues \[[\s\S]+?)(?:\n\n|\Z)',
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'eigenvalues',
                            r'(?:(#st)\s*Spin\s*Eigenvalue\s*Occupation|#k =\s*\d+, '
                            r'k = \(([\-\d\.,\s]+)\))([\d\s\.\-updn]+)',
                            str_operation=str_to_eigenvalues,
                            repeats=True,
                            convert=False,
                        ),
                        Quantity(
                            'unit',
                            r'Eigenvalues \[(.*)\]',
                            convert=False,
                            str_operation=lambda x: 'eV' if x == 'eV' else 'hartree',
                        ),
                        Quantity(
                            'fermi_energy',
                            rf'Fermi energy =\s*({RE_FLOAT} .*)',
                            convert=False,
                            str_operation=str_to_fermi_energy,
                        ),
                    ]
                ),
            )
        ]


class InfoParser(EigenvalueParser):
    def init_quantities(self) -> None:
        super().init_quantities()

        def str_to_energies(val_in: str) -> dict[str, pint.Quantity]:
            val = [v.split('=') for v in val_in.strip().split('\n')]
            unit = ureg.eV if '[eV]:' in val[0] else ureg.hartree
            val_n = 2
            return {v[0].strip(): float(v[1]) * unit for v in val if len(v) == val_n}

        def str_to_forces(val_in: str) -> pint.Quantity:
            val = [v.split() for v in val_in.strip().split('\n')]
            unit = (
                ureg.eV / ureg.angstrom
                if '[eV/A]' in val[0]
                else ureg.hartree / ureg.bohr
            )
            nval = 5
            return np.array([v[2:5] for v in val if len(v) == nval], dtype=float) * unit

        self._quantities.extend(
            [
                Quantity(
                    'brillouin_zone_sampling',
                    r'Brillouin zone sampling([\s\S]+?)\*+\n\n',
                    sub_parser=TextParser(
                        quantities=[
                            Quantity(
                                'kgrid',
                                r'Dimensions of the k\-point grid\s*=(.*)',
                                dtype=int,
                            ),
                            Quantity(
                                'n_kpoints',
                                r'Total number of k\-points\s*=(.*)',
                                dtype=int,
                            ),
                            Quantity(
                                'n_kpoints_reduced',
                                r'Number of symmetry-reduced k\-points\s*=(.*)',
                                dtype=int,
                            ),
                            Quantity(
                                'kpoints',
                                rf'\d+\s*({RE_FLOAT})\s*({RE_FLOAT})\s*({RE_FLOAT})\s*({RE_FLOAT})',
                                dtype=float,
                                repeats=True,
                            ),
                        ]
                    ),
                ),
                Quantity(
                    'energies',
                    r'Energy (\[[\s\S]+?)\n\n',
                    str_operation=str_to_energies,
                ),
                Quantity(
                    'total_magnetic_moment',
                    rf'Total Magnetic Moment:\s*mz\s*=\s*({RE_FLOAT})',
                ),
                Quantity(
                    'local_magnetic_moments',
                    r'Ion\s*mz\s*([\w\s\.\-]+?)\n\n',
                    convert=False,
                    str_operation=lambda x: np.array(
                        [v.split()[2] for v in x.strip().split('\n')], dtype=float
                    ),
                ),
                Quantity(
                    'dipole',
                    rf'Dipole:.*\[Debye\]\s*<x> =\s*\S+\s*({RE_FLOAT})\s*'
                    rf'<y> =\s*\S+\s*({RE_FLOAT})\s*<z> =\s*\S+\s*({RE_FLOAT})',
                    dtype=float,
                ),
                Quantity(
                    'forces',
                    r'Forces on the ions (\[.*\]\s*)Ion\s*x\s*y\s*z\s*([\s\S]*)',
                    str_operation=str_to_forces,
                    convert=False,
                ),
            ]
        )


class ControlParser(TextParser):
    _constants = {
        'pi': np.pi,
        'angstrom': 1.0 / ase_units.Bohr,
        'ev': 1.0 / ase_units.Hartree,
        'yes': True,
        'no': False,
        't': True,
        'f': False,
        'true': True,
        'false': False,
        'i': 1j,
    }
    _re_sqrt = re.compile(r'sqrt([\w\. ]+\Z)')

    def __init__(self):
        super().__init__()
        self._info = None
        self._keys_mapping = dict()

    def init_quantities(self) -> None:
        def str_to_line(val_in: str) -> list[str]:
            val = val_in.replace('"', '').replace("'", '').split('=', 1)
            return [v.strip().split('#')[0] for v in val]

        self._quantities = [
            Quantity(
                'line', r'(\w.*\s*=\s*.*)#?', str_operation=str_to_line, repeats=True
            )
        ]

    def reset(self):
        super().reset()
        self._info = None
        self._keys_mapping = dict()

    def evaluate_value(self, value: Any) -> Any:
        """
        Evaluate octopus parameterized variables, e.g. 2*angstrom
        """
        # TODO implement all operations and parameters supported by Octopus

        if isinstance(value, list):
            return [self.evaluate_value(v) for v in value]

        if not isinstance(value, str):
            return value

        try:
            return float(value)
        except Exception:
            pass

        value = value.strip()
        if value == '':
            output = 1

        elif (val := self._constants.get(value.lower())) is not None:
            output = val

        elif (val := self._info.get(self._keys_mapping.get(value.lower()))) is not None:
            output = val

        elif (open_p := value.rfind('(')) > -1:
            n_groups = value.count('(')
            if n_groups != value.count(')'):
                output = value
            else:
                for _ in range(n_groups):
                    part = value[open_p + 1 :]
                    part = part[: part.find(')')]
                    value = value.replace(f'({part})', str(self.evaluate_value(part)))
                    open_p = value.rfind('(')
                output = self.evaluate_value(value)
        else:
            output = self._evaluate(value)

        return output

    def _evaluate(self, value: str) -> float:
        if sqrt := self._re_sqrt.match(value):
            val = self.evaluate_value(sqrt.group(1))
            output = np.sqrt(val)

        elif len(vals := value.split('**')) > 1:
            vals = [self.evaluate_value(v) for v in vals]
            val = vals[0]
            for v in reversed(vals[1:]):
                val = val**v
            output = val

        elif len(vals := value.split('*')) > 1:
            vals = [self.evaluate_value(v) for v in vals]
            output = np.prod(vals)

        elif len(vals := value.split('/')) > 1:
            vals = [self.evaluate_value(v) for v in vals]
            val = vals[0]
            for v in vals[1:]:
                val /= v
            output = val

        elif '+' in value:
            vals = value.split('+')
            vals = [self.evaluate_value(v) for v in vals]
            val = 0.0
            for v in vals:
                val += v
            output = val

        elif '-' in value:
            vals = value.split('-')
            vals = [self.evaluate_value(v) for v in vals]
            val = 0.0
            for v in vals:
                val -= v
            output = val

        return output

    @property
    def info(self):
        if self._info is None:
            self._info = {v[0].strip(): v[1] for v in self.get('line', [])}
            self._info.update({v[0].strip(): v[1:] for v in self.get('block', [])})
            self._keys_mapping = {k.lower(): k for k in self._info.keys()}
            for key, val in self._info.items():
                try:
                    self._info[key] = self.evaluate_value(val)
                except Exception:
                    self._info[key] = val
                self._keys_mapping[key.lower()] = key
        return self._info


class InpParser(ControlParser):
    def init_quantities(self) -> None:
        def str_to_block(val_in: str) -> list[str]:
            val = [v.split('#')[0] for v in val_in.strip().split('\n')]
            val = [v.replace('"', '').replace("'", '').split('|') for v in val if v]
            val[0] = val[0][0]
            return val

        super().init_quantities()
        self._quantities.extend(
            [
                Quantity(
                    'block', r'%([\s\S]+?)%', repeats=True, str_operation=str_to_block
                )
            ]
        )

    def get_coordinates(self) -> tuple[list[str, np.ndarray]]:
        val = self.info.get('Coordinates', self.info.get('ReducedCoordinates', []))

        symbols = []
        coordinates = []
        for v in val:
            symbols.append(v[0].strip())
            coordinates.append(v[1:])

        coordinates = np.array(coordinates, dtype=float)
        coordinates.shape = (len(symbols), 3)
        return symbols, coordinates


class LogParser(ControlParser):
    def init_quantities(self) -> None:
        def str_to_block(val_in: str) -> tuple[str, list[str]]:
            val = val_in.strip().split('\n')
            name = val[0].strip().replace('"', '').replace("'", '')
            val = [
                v.split('#')[0].split('=')[1].replace('"', '').replace("'", '').strip()
                for v in val[1:]
            ]
            return name, val

        super().init_quantities()
        self._quantities.extend(
            [
                Quantity(
                    'block',
                    r'Opened block([\s\S]+?)Closed block',
                    repeats=True,
                    str_operation=str_to_block,
                )
            ]
        )

    def get_coordinates(self) -> tuple[str, np.ndarray]:
        symbols = []
        coordinates = []

        val = self.info.get('Coordinates', self.info.get('ReducedCoordinates', [[]]))[0]

        for v in val:
            if v[0].isdecimal() or not isinstance(v, str):
                coordinates.append(v)
            else:
                symbols.append(v.strip())

        coordinates = np.array(coordinates, dtype=float)
        coordinates.shape = (len(symbols), 3)

        return symbols, coordinates


class OutParser(TextParser):
    def init_quantities(self) -> None:
        def str_to_option(val_in: str) -> list[str]:
            val = val_in.strip().split(':')
            return [val[0].strip(), ''.join(val[1:]).strip()]

        def str_to_cell(val_in: str) -> pint.Quantity:
            unit = ureg.angstrom if val_in.lower().startswith('a') else ureg.bohr
            val = [v.split() for v in val_in.strip().split('\n')[1:]]
            return np.array(val, dtype=float) * unit

        def str_to_spacing(val_in: str) -> pint.Quantity:
            unit = ureg.angstrom if val_in.startswith('A') else ureg.bohr
            return np.array(val_in.split()[1:4], dtype=float) * unit

        def str_to_energy(val_in: str) -> pint.Quantity:
            val = val_in.split()
            unit = ureg.eV if val[1].startswith('e') else ureg.hartree
            return float(val[0]) * unit

        def str_to_td_iteration(val_in: str) -> dict[str, Any]:
            val = val_in.strip().split()
            return dict(
                iter=int(val[0]),
                time=float(val[1]),
                energy=float(val[2]),
                scfsteps=int(val[3]),
                elapsed_time=float(val[4]),
            )

        iteration_quantities = [
            Quantity('energy_total', rf'etot\s*=\s*({RE_FLOAT})'),
            # TODO scf_iteration eigenvalues are sometimes truncated and unusable
            Quantity(
                'fermi_level',
                rf'Fermi energy\s*=\s*({RE_FLOAT} .*)',
                str_operation=str_to_energy,
                convert=False,
            ),
            Quantity(
                'time',
                r'Elapsed time for SCF step\s*\d+:\s*([\d\.]+)',
                unit='s',
                dtype=float,
            ),
        ]

        self._quantities = [
            Quantity(
                'header',
                r'Running octopus([\s\S]+?)\*{10}',
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'options',
                            r'\n *([\w ]+?\s*:\s*.*)',
                            str_operation=str_to_option,
                            repeats=True,
                        ),
                    ]
                ),
            ),
            Quantity(
                'grid',
                r'\*\s*Grid\s*\*+\s*([\s\S]+?)\*{10}',
                sub_parser=TextParser(
                    quantities=[
                        Quantity('boxshape', r'Type\s*=\s*(.*)'),
                        Quantity(
                            'npbc',
                            r'Octopus will treat the system as periodic in (\S+) dim',
                            dtype=int,
                        ),
                        Quantity(
                            'cell',
                            r'Lattice Vectors \[(.*)\]([-\d\s\.]+)',
                            str_operation=str_to_cell,
                            convert=False,
                        ),
                        Quantity(
                            'spacing',
                            r'Spacing \[(.*)\] = \(\s*(\S+), (\S+), (\S+)\s*\)',
                            str_operation=str_to_spacing,
                            convert=False,
                        ),
                    ]
                ),
            ),
            Quantity(
                'theory_level',
                r'\*\s*Theory Level\s*\*+\s*([\s\S]+?)\*{10}',
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'theory_level', r'\[TheoryLevel = (.+)\]', flatten=False
                        ),
                        Quantity('exchange', r'Exchange\s+(.*) \(', flatten=False),
                        Quantity(
                            'correlation', r'Correlation\s+(.*) \(', flatten=False
                        ),
                    ]
                ),
            ),
            Quantity(
                'self_consistent',
                r'Info: Starting SCF iteration\.\s*([\s\S]+?)Info: SCF',
                repeats=True,
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'iteration',
                            r'SCF CYCLE ITER #\s*(\d+\s*\*+[\s\S]+?)\*{10}',
                            repeats=True,
                            sub_parser=TextParser(quantities=iteration_quantities),
                        )
                    ]
                ),
            ),
            Quantity(
                'time_dependent',
                r'Time\-Dependent Simulation \*+([\s\S]+?)Info: '
                r'Finished writing information',
                repeats=True,
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'iteration',
                            rf'\n *(\d+\s*{RE_FLOAT}\s*{RE_FLOAT}\s*\d+\s*{RE_FLOAT})'
                            rf' *\n',
                            str_operation=str_to_td_iteration,
                            repeats=True,
                            convert=False,
                        )
                    ]
                ),
            ),
            Quantity(
                'x_octopus_info_scf_converged_iterations',
                r'SCF converged in\s*(\d+) iterations',
                dtype=int,
            ),
            Quantity(
                'minimization',
                r'(MINIMIZATION ITER #:\s*\d+\s*\++\s*Energy[\s\S]+?\+{10})',
                repeats=True,
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'energy_total',
                            rf'Energy\s*=\s*({RE_FLOAT} .*)',
                            str_operation=str_to_energy,
                            convert=False,
                        ),
                        Quantity('number', r'ITER #:\s*(\d+)', dtype=int),
                    ]
                ),
            ),
            # calculation results are not printed in outfile but in info
        ]

        self._header = None

    @property
    def header(self):
        return {k: v for k, v in self.get('header', {}).get('options', [])}
