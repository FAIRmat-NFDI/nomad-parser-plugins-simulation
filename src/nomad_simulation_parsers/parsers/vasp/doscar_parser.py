# ruff: noqa: PLR2004
import os
import re
from typing import Any

import numpy as np
from nomad_file_parser import FileParser
from nomad_file_parser.mapping_parser import MappingParser


class DOSCARFileParser(FileParser):
    """Parser for VASP and Lobster DOSCAR files."""

    re_atomic_number = re.compile(r'.+?Z=\s*(\d+).+')

    def parse(self, key=None):  # noqa: PLR0912, PLR0915, PLR2004
        if self._results is None:
            self._results = {}

        n_header = 5
        n_dos = 0
        dos = np.array([], dtype=np.float64)
        pdos = []
        pdos_atom = np.array([], dtype=np.float64)
        atomic_numbers = []
        with self.open_mainfile_obj() as f:
            for n, line in enumerate(f):
                if n < n_header:
                    continue

                line_split = line.strip().split()
                if (
                    n > n_header + n_dos
                    and len(line_split) > 2
                    and line_split[2].isdigit()
                    and int(line_split[2]) == n_dos
                    and 'Z=' not in line
                ):
                    continue

                try:
                    value = np.array(line_split, dtype=np.float64)
                except Exception:
                    atomic_number = self.re_atomic_number.match(line)
                    if atomic_number:
                        atomic_numbers.append(int(atomic_number.group(1)))
                    continue

                if n == n_header:
                    # VASP's DOSCAR header is Emax, Emin, NEDOS, Efermi
                    # (https://www.vasp.at/wiki/index.php/DOSCAR).
                    # Older Lobster files may omit the latter fields.
                    self._results['e_fermi'] = float(
                        value[3] if len(value) > 3 else value[1]
                    )
                    n_dos = int(value[2])
                    continue

                if n <= n_dos + n_header:
                    dos = value if not len(dos) else np.vstack((dos, value))
                else:
                    value = value[1:]
                    pdos_atom = (
                        value if not len(pdos_atom) else np.vstack((pdos_atom, value))
                    )
                    if len(pdos_atom) == n_dos:
                        pdos.append(pdos_atom)
                        pdos_atom = np.array([], dtype=np.float64)

        if len(dos) == 0:
            return

        # DOSCAR format: energy, DOS (up/down), integrated DOS (up/down).
        dos = np.transpose(dos)
        # Standard total DOS has one energy column and two columns per spin
        # (DOS and integrated DOS). Keep one channel for reduced test/minimal
        # files that only contain energy and DOS.
        n_spin = max(1, (np.shape(dos)[0] - 1) // 2)

        if atomic_numbers:
            self._results['atomic_numbers'] = atomic_numbers
            self._results['pbc'] = [True, True, True]
        self._results['energies'] = dos[0]
        self._results['total_dos'] = dos[1:]
        self._results['projected_dos'] = []

        if len(pdos) == 0:
            return

        for atom_index, pdos_n in enumerate(pdos):
            pdos_atom = np.transpose(pdos_n)
            n_lm = np.shape(pdos_atom)[0] // n_spin
            pdos_atom = np.reshape(pdos_atom, (n_lm, n_spin, n_dos))
            pdos_dict = dict(
                dos=np.transpose(pdos_atom, (1, 0, 2)), atom_index=atom_index
            )
            if n_lm == 1:
                fields = ['s']
            elif n_lm == 3:
                fields = ['s', 'p', 'd']
            elif n_lm == 9:
                fields = ['s', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2', 'dxz', 'dx2']
            elif n_lm == 16:
                fields = [
                    's',
                    'py',
                    'pz',
                    'px',
                    'dxy',
                    'dyz',
                    'dz2',
                    'dxz',
                    'dx2',
                    'f-3',
                    'f-2',
                    'f-1',
                    'f0',
                    'f1',
                    'f2',
                    'f3',
                ]
            else:
                fields = [None] * n_lm
            pdos_dict['dos_fields'] = fields
            self._results['projected_dos'].append(pdos_dict)


class DOSCARParser(MappingParser):
    """Mapping parser for the DOSCAR associated with a VASP calculation."""

    def load_file(self) -> FileParser | None:
        maindir = os.path.dirname(self.filepath)
        outcar_suffix = os.path.basename(self.filepath).removeprefix('OUTCAR')
        doscar_path = os.path.join(maindir, f'DOSCAR{outcar_suffix}')
        if not os.path.isfile(doscar_path):
            doscar_path = os.path.join(maindir, 'DOSCAR')
        if not os.path.isfile(doscar_path):
            return None

        doscar_parser = DOSCARFileParser()
        doscar_parser.mainfile = doscar_path
        return doscar_parser

    def to_dict(self, **kwargs) -> dict:
        if not self.data_object:
            return {}

        self.data_object.parse()
        return self.data_object._results or {}

    def get_dos(
        self, total: np.ndarray, projected: list[dict[str, Any]]
    ) -> list[dict[str, Any]]:
        if len(total) % 2 != 0:
            self.logger.error(
                'Incorrect shape for total DOS.', extra=dict(file=self.filepath)
            )
            return []

        dos = []
        n_spin = len(total) // 2
        has_negative_dos = False
        for spin in range(n_spin):
            if np.any(total[spin] < 0):
                has_negative_dos = True
            dct = dict(
                dos=np.maximum(total[spin], 0),
                integrated=total[n_spin + spin],
                spin=spin if n_spin == 2 else None,
                projected=[],
            )
            for atom, pdos_dict in enumerate(projected):
                for lm, pdos_lm in enumerate(pdos_dict['dos'][spin]):
                    if np.any(pdos_lm < 0):
                        has_negative_dos = True
                    dct['projected'].append(
                        dict(
                            dos=np.maximum(pdos_lm, 0),
                            lm=lm,
                            atom_index=atom,
                            spin=spin if n_spin == 2 else None,
                        )
                    )
            dos.append(dct)
        if has_negative_dos:
            self.logger.warning('Found negative values in DOS, setting them to zero.')
        return dos

    def from_dict(self, dct: dict):
        pass
