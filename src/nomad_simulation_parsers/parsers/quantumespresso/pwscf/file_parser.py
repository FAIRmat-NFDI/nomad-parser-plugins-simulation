import re
from typing import Any

import numpy as np
import pint
from nomad.units import ureg
from nomad_file_parser.text_parser import DataTextParser, Quantity, TextParser

from ..common import (
    RE_FLOAT,
    calculation_quantities,
    general_quantities,
    header_quantities,
    scf_iteration_quantities,
    str_to_atom_data,
    tail_quantities,
)


def str_to_energy_contributions(val_in: str) -> dict[str, pint.Quantity]:
    val = [v.split('=') for v in val_in.strip().split('\n')]
    return {v[0].strip(): float(v[1].split()[0]) * ureg.rydberg for v in val}


def str_to_forces(val_in: str) -> pint.Quantity:
    val = val_in.strip().split('\n')
    val = [v.split('=')[1].split() for v in val if 'force =' in v]
    return np.array(val, dtype=float) * ureg.rydberg / ureg.bohr


def str_to_stress(val_in: str) -> tuple[pint.Quantity, pint.Quantity]:
    val = [v.split() for v in val_in.strip().split('\n')]
    pressure = float(val[0][0]) * ureg.kilobar
    stress = np.array([v[3:6] for v in val[1:4]], dtype=float) * ureg.kilobar
    return pressure, stress


def str_to_labels_positions(val_in: str) -> dict[str, Any]:
    units = re.search(r'ATOMIC_POSITIONS \((.+)\)', val_in)
    out = dict()
    if units:
        out['units'] = units.group(1)
    val = [v.split()[:4] for v in val_in.strip().split('\n')[1:]]
    val = np.transpose([v for v in val if v[1][-1].isdecimal()])
    out['labels'] = val[0]
    out['positions'] = np.array(val[1:], dtype=float).T
    return out


def str_to_algorithm(val_in: str) -> str:
    val = val_in.strip().lower()
    return {'cg style': 'conjugate_gradient'}.get(val, val)


def str_to_band_energies(val_in: str) -> np.ndarray:
    out = [line.strip().split() for line in val_in.strip().splitlines()]
    return np.array([item for line in out for item in line], dtype=np.float64)


class PWSCFFileParser(TextParser):
    def init_quantities(self) -> None:
        diagonalization_quantities = [
            Quantity(
                'diagonalization_algorithm',
                r'(Davidson|CG style)\s*diagonalization',
                str_operation=str_to_algorithm,
                convert=False,
            ),
            Quantity('diagonalization_ethr', r'ethr =\s*([\d\.\-E]+)', dtype=float),
            Quantity(
                'diagonalization_iteration_avg',
                r'avg # of iterations\s*=\s*([\d\.]+)',
                dtype=float,
            ),
            Quantity(
                'diagonalization_c_bands_n_unconverged_eigenvalues',
                r'c_bands:\s*(\d+)\s*eigenvalues not converged',
                dtype=int,
                repeats=True,
            ),
        ]

        output_quantities = (
            general_quantities
            + calculation_quantities
            + [
                Quantity('spin_pol', r'SPIN (UP|DOWN)'),
                Quantity(
                    'k_points',
                    rf'k\s*=\s*({RE_FLOAT})\s*({RE_FLOAT})\s*({RE_FLOAT})',
                    repeats=True,
                    dtype=float,
                ),
                Quantity(
                    'number_of_planewaves',
                    r'\(\s*(\d+)\s*PWs\)',
                    dftype=int,
                    repeats=True,
                ),
                Quantity(
                    'band_energies',
                    r'band(?:s|\s+energies)\s*\(\s*[eE][vV]\s*\)\s*:\s*([\d\.\-\s]+)',
                    str_operation=str_to_band_energies,
                    repeats=True,
                    convert=False,
                ),
                Quantity(
                    'occupation_numbers',
                    r'occupation numbers\s*([\d\.\-\s]+?)\n *\n',
                    repeats=True,
                    dtype=float,
                ),
                Quantity(
                    'homo_lumo',
                    r'highest occupied(?:, lowest unoccupied)* level '
                    r'\(ev\):\s*([\-\d\. ]+)',
                    dtype=float,
                ),
                Quantity(
                    'fermi_energy',
                    r'(?:the Fermi energy is|the spin up/dw Fermi energies are)'
                    r'\s*([\-\d\. ]+)',
                    dtype=float,
                ),
                Quantity(
                    'energy_contributions',
                    r'The total energy is.+?the sum of the following terms:'
                    r'\s*([\s\S]+?Ry\n\s*\n)',
                    str_operation=str_to_energy_contributions,
                    convert=False,
                ),
                Quantity(
                    'magnetization_total',
                    rf'total magnetization\s*=\s*({RE_FLOAT})\s*Bohr mag/cell',
                    dtype=float,
                    unit='bohr_magneton',
                    repeats=True,
                ),
                Quantity(
                    'magnetization_absolute',
                    rf'absolute magnetization\s*=\s*({RE_FLOAT})\s*Bohr mag/cell',
                    dtype=float,
                    unit='bohr_magneton',
                    repeats=True,
                ),
                Quantity(
                    'convergence_iterations',
                    r'convergence has been achieved in\s*([\d]+) iterations',
                    dtype=int,
                ),
                Quantity(
                    'forces',
                    r'Forces acting on atoms \((?:cartesian axes, )?'
                    r'Ry\/au\):\s*([\s\S]+?)(?:The|\n\s*\n)',
                    str_operation=str_to_forces,
                    convert=False,
                ),
                Quantity(
                    'total_force',
                    rf'Total force\s*=\s*({RE_FLOAT})\s*Total SCF correction\s*=\s*'
                    rf'({RE_FLOAT})',
                    dtype=float,
                    unit='rydberg/bohr',
                ),
                Quantity(
                    'forces_dispersion',
                    r'Dispersion forces acting on atoms:\s*([\s\S]+?)(?:The|\n\s*\n)',
                    str_operation=str_to_forces,
                    convert=False,
                ),
                Quantity(
                    'total_force_dispersion',
                    rf'Total Dispersion Force =\s*({RE_FLOAT})',
                    dtype=float,
                    unit='rydberg/bohr',
                ),
                Quantity(
                    'stress',
                    r'total\s*stress\s*\(Ry/bohr\*\*3\)\s*\(kbar\)\s*P=\s*([\s\S]+?)\n\s*\n',
                    str_operation=str_to_stress,
                    convert=False,
                ),
                Quantity(
                    'units', r'crystal axes: \(cart\. coord\. in units of ([\w ]+)\)\s*'
                ),
                Quantity(
                    'simulation_cell',
                    r'a\(1\) = \(([\-\d\. ]+)\)\s*a\(2\) = \(([\-\d\. ]+)\)\s*a\(3\) = '
                    r'\(([\-\d\. ]+)\)\s*',
                    dtype=float,
                    shape=(3, 3),
                ),
                Quantity(
                    'reciprocal_cell_units',
                    r'reciprocal axes: \(cart\. coord\. in units ([\w \/]+)\)',
                    flatten=False,
                ),
                Quantity(
                    'reciprocal_cell',
                    r'b\(1\) = \(([\-\d\. ]+)\)\s*b\(2\) = \(([\-\d\. ]+)\)\s*b\(3\) = '
                    r'\(([\-\d\. ]+)\)\s*',
                    dtype=float,
                    shape=(3, 3),
                ),
                Quantity(
                    'labels_positions',
                    r'(ATOMIC_POSITIONS \(.+\)[\s\S]+?)\n\s*\n',
                    str_operation=str_to_labels_positions,
                    convert=False,
                ),
                Quantity(
                    'starting_magnetization',
                    r'Starting magnetic structure\s*atomic species\s*magnetization'
                    r'([\s\S]+?)\n\s*\n',
                    str_operation=str_to_atom_data,
                    convert=False,
                ),
                Quantity(
                    'exx_refine', r'(EXX: now go back to refine exchange calculation)'
                ),
            ]
        )

        scf_quantities = [
            Quantity(
                'iteration',
                r'( # *\d+[\s\S]+?(?:\n *iter|End))',
                repeats=True,
                sub_parser=TextParser(
                    quantities=scf_iteration_quantities
                    + calculation_quantities
                    + diagonalization_quantities
                ),
            )
        ] + output_quantities

        # TODO add electric field calculation

        bandstructure_quantities = diagonalization_quantities + output_quantities

        sampling_quantities = [
            Quantity(
                'self_consistent',
                r'(consistent Calculation[\s\S]+?)(?:Self-|init_run|\Z)',
                repeats=True,
                sub_parser=TextParser(quantities=scf_quantities),
            ),
            Quantity('dynamics', r'(Entering Dynamics)', repeats=False),
        ]

        bfgs_quantities = [
            Quantity(
                'final_energy',
                rf'Final energy\s*=\s*({RE_FLOAT})\s*Ry',
                dtype=float,
                unit='rydberg',
            ),
            Quantity(
                'convergence',
                r'bfgs converged in\s*(\d+)\s*scf cycles and\s*(\d+)\s*bfgs steps',
                dtype=int,
            ),
        ] + sampling_quantities

        md_quantities = [
            Quantity(
                'diffusion_coefficients',
                r'atom\s*\d+\s*D =\s*({re_float})\s*cm\^2/s',
                repeats=True,
                unit='cm**2/s',
            ),
            Quantity(
                'diffusion_coefficient_mean',
                r'< D > =\s*({re_float})\s*cm\^2/s',
                unit='cm**2/s',
            ),
        ] + sampling_quantities

        langevin_quantities = sampling_quantities

        vcs_quantities = sampling_quantities

        damped_quantities = sampling_quantities

        self._quantities = [
            Quantity(
                'header',
                r'([Pp]rogram PWSCF[\s\S]+?)(?:Self\-|Band)',
                repeats=False,
                sub_parser=TextParser(
                    quantities=header_quantities + general_quantities
                ),
            ),
            Quantity(
                'self_consistent',
                r'(consistent Calculation[\s\S]+?(?:Self-|init_run|\Z))',
                repeats=False,
                sub_parser=TextParser(quantities=scf_quantities),
            ),
            Quantity(
                'bandstructure',
                r'(Structure Calculation[\s\S]+?)(?:init_run|\Z)',
                repeats=False,
                sub_parser=TextParser(quantities=bandstructure_quantities),
            ),
            Quantity(
                'bfgs_geometry_optimization',
                r'(S Geometry Optimization[\s\S]+?)(?:init_run|\Z)',
                repeats=False,
                sub_parser=TextParser(quantities=bfgs_quantities),
            ),
            Quantity(
                'molecular_dynamics',
                r'(r Dynamics Calculation[\s\S]+?)(?:init_run|\Z)',
                repeats=False,
                sub_parser=TextParser(quantities=md_quantities),
            ),
            Quantity(
                'damped_dynamics',
                r'(d Dynamics Calculation[\s\S]+?)(?:init_run|\Z)',
                repeats=False,
                sub_parser=TextParser(quantities=damped_quantities),
            ),
            Quantity(
                'langevin_dynamics',
                r'(d Langevin Dynamics Calculation[\s\S]+?)(?:init_run|\Z)',
                repeats=False,
                sub_parser=TextParser(quantities=langevin_quantities),
            ),
            Quantity(
                'vcs_wentzcovitch_damped_minimization',
                r'(h Damped Cell[\- ]*Dynamics Minimization:[\s\S]+?)(?:init_run|\Z)',
                repeats=False,
                sub_parser=TextParser(quantities=vcs_quantities),
            ),
        ] + tail_quantities


class PWSCFDOSTextParser(DataTextParser):
    MIN_DOS_COLUMNS = 2
    DOS_ARRAY_NDIM = 2

    def parse(self, key=None):
        super().parse(key)
        data = self._results.pop('data', None)
        if data is not None:
            if data.ndim == 1 and data.size >= self.MIN_DOS_COLUMNS:
                data = data.reshape(1, -1)
            if (
                data.ndim == self.DOS_ARRAY_NDIM
                and data.shape[1] >= self.MIN_DOS_COLUMNS
            ):
                self._results['energies'] = data[:, 0] * ureg.eV
                self._results['value'] = np.abs(data[:, 1]) / ureg.eV
