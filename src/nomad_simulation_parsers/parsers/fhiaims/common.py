from typing import Any

import numpy as np
from ase import Atoms
from nomad.parsing.file_parser.text_parser import Quantity, TextParser


class GeometryParser(TextParser):
    """
    Parser for the FHI-aims geometry file.
    """

    def init_quantities(self) -> None:
        self._quantities = [
            Quantity(
                'lattice_vector',
                r'lattice_vector +(\S+) +(\S+) +(\S+)',
                repeats=True,
                dtype=np.float64,
            ),
            Quantity(
                'atom',
                r'(atom.+)',
                repeats=True,
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'fractional', r'(atom_frac)', str_operation=lambda x: True
                        ),
                        Quantity(
                            'position',
                            r' +(\S+) +(\S+) +(\S+) +',
                        ),
                        Quantity('label', r'([A-Z][a-z]*)'),
                    ]
                ),
            ),
            Quantity('magmom', r'initial_moment (.+)', dtype=np.float64),
        ]

    def get_atoms(self) -> Atoms:
        atoms = Atoms(
            cell=self.lattice_vector,
            pbc=self.lattice_vector is not None,
            positions=[
                atom.position
                if not atom.fractional
                else np.dot(atom.position, self.lattice_vector)
                for atom in self.atom
            ],
            symbols=[atom.label for atom in self.atom],
            magmoms=self.magmom,
        )
        # TODO compatibility with PhonopyAtoms
        add_attribs = dict(
            masses=atoms.get_masses(),
            magnetic_moments=None,
            scaled_positions=atoms.get_scaled_positions(),
        )
        for key, val in add_attribs.items():
            setattr(atoms, key, val)
        return atoms


class ControlParser(TextParser):
    def __init__(self):
        super().__init__()

    def init_quantities(self) -> None:
        def str_to_nac(val_in: str) -> dict[str, Any]:
            val = val_in.strip().split()
            nac = dict(file=val[0], method=val[1].lower())
            n_val = 2
            if len(val) > n_val:
                nac['delta'] = [float(v) for v in val[3:6]]
            return nac

        def str_to_supercell(val_in: str) -> np.ndarray:
            val = [int(v) for v in val_in.strip().split()]
            n_val = 3
            if len(val) == n_val:
                return np.diag(val)
            else:
                return np.reshape(val, (3, 3))

        def str_to_int_array(val_in: str) -> np.ndarray:
            """Parse space-separated integers into numpy array."""
            return np.array([int(v) for v in val_in.strip().split()], dtype=np.int32)

        def str_to_float_array(val_in: str) -> np.ndarray:
            """Parse space-separated floats into numpy array."""
            return np.array(
                [float(v) for v in val_in.strip().split()], dtype=np.float64
            )

        self._quantities = [
            Quantity(
                'displacement', r'\n *phonon\s*displacement\s*([\d\.]+)', dtype=float
            ),
            Quantity(
                'symmetry_thresh',
                r'\n *phonon symmetry_thresh\s*([\d\.]+)',
                dtype=float,
            ),
            Quantity('frequency_unit', r'\n *phonon frequency_unit\s*(\S+)'),
            Quantity(
                'supercell',
                r'\n *phonon\s*supercell\s*(.+)',
                str_operation=str_to_supercell,
            ),
            Quantity('nac', r'\n *phonon nac\s*(.+)', str_operation=str_to_nac),
            Quantity(
                'k_grid',
                r'\n *k_grid\s*([\d ]+)',
                str_operation=str_to_int_array,
                convert=False,
            ),
            Quantity(
                'k_offset',
                r'\n *k_offset\s*([-+\d\. ]+)',
                str_operation=str_to_float_array,
                convert=False,
            ),
        ]
