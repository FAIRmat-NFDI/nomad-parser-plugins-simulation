import numpy as np
from nomad.parsing.file_parser.text_parser import Quantity, TextParser


def str_to_eigenvalues(val_in: str) -> dict[str, np.ndarray]:
    val = val_in[: val_in.rfind('\n \n')].strip()
    val = np.array([v.split() for v in val.split('\n')], dtype=float)
    val = np.transpose(val)
    occs = val[-1]
    eigs = val[-2]

    tol = 0.1
    nspin = 1 if np.any(occs > 1 + tol) else 2
    data = dict()
    data['occupancies'] = np.reshape(occs, (nspin, len(occs) // nspin))
    data['eigenvalues'] = np.reshape(eigs, (nspin, len(eigs) // nspin))
    return data


class EigvalReader(TextParser):
    def init_quantities(self):
        self._quantities = [
            Quantity('k_points', r'\s*\d+\s*([\d\.Ee\- ]+):\s*k\-point', repeats=True),
            Quantity(
                'eigenvalues_occupancies',
                r'\(state\, eigenvalue and occupancy below\)\s*'
                r'([\d\.Ee\-\s]+?(?:\n *\n))',
                str_operation=str_to_eigenvalues,
                repeats=True,
            ),
            Quantity('n_k_points', r'(\d+) +\: +nkpt', dtype=int),
            Quantity('n_states', r'(\d+) +\: +nstsv', dtype=int),
        ]
