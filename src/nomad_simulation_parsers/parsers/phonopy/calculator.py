#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD.
# See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
import re
from fractions import Fraction
from typing import Any

import numpy as np
from nomad_simulations.schema_packages.numerical_settings import KSpaceFunctionalities
from phonopy import Phonopy
from phonopy.phonon.band_structure import BandStructure
from phonopy.physical_units import get_physical_units
from phonopy.structure.atoms import PhonopyAtoms

_PHONOPY_UNITS = get_physical_units()
EvTokJmol = _PHONOPY_UNITS.EvTokJmol
VaspToTHz = _PHONOPY_UNITS.DefaultToTHz


def generate_kpath_parameters(
    points: dict[str, np.ndarray], paths: list[list[str]], npoints: int
) -> list[dict[str, Any]]:
    k_points: list[list[np.ndarray]] = []
    for p in paths:
        k_points.append([points[k] for k in p])
        for index in range(len(p)):
            if p[index] in ('G', 'GAMMA'):
                p[index] = 'Γ'
    parameters: list[dict[str, Any]] = []
    n_k = 2
    for h, seg in enumerate(k_points):
        for i, path in enumerate(seg):
            parameter: dict[str, Any] = {}
            parameter['npoints'] = npoints
            parameter['startname'] = paths[h][i]
            if i == 0 and len(seg) > n_k:
                parameter['kstart'] = path
                parameter['kend'] = seg[i + 1]
                parameter['endname'] = paths[h][i + 1]
                parameters.append(parameter)
            elif i == (len(seg) - 2):
                parameter['kstart'] = path
                parameter['kend'] = seg[i + 1]
                parameter['endname'] = paths[h][i + 1]
                parameters.append(parameter)
                break
            else:
                parameter['kstart'] = path
                parameter['kend'] = seg[i + 1]
                parameter['endname'] = paths[h][i + 1]
                parameters.append(parameter)
    return parameters


def read_kpath(filename: str) -> list[dict[str, Any]]:
    with open(filename) as f:
        string = f.read()

        labels_extracted = re.search(r'BAND_LABELS\s*=\s*(.+)', string)
        try:
            labels = labels_extracted.group(1).strip().split()
        except Exception:
            return []

        points_extracted = re.search(r'BAND\s*=\s*(.+)', string)
        try:
            points = points_extracted.group(1)
            points = [float(Fraction(p)) for p in points.split()]
            points = np.reshape(points, (len(labels), 3))
            points = {labels[i]: points[i] for i in range(len(labels))}
        except Exception:
            return []

        npoints_extracted = re.search(r'BAND_POINTS\s*\=\s*(\d+)', string)
        npoints = 100 if npoints_extracted is None else int(npoints_extracted.group(1))

    return generate_kpath_parameters(points, [labels], npoints)


def generate_kpath_seekpath(
    atoms: PhonopyAtoms, symprec: float, logger=None
) -> list[dict[str, Any]]:
    structure = (
        atoms.cell.tolist(),
        atoms.scaled_positions.tolist(),
        atoms.numbers.tolist(),
    )
    resolved = KSpaceFunctionalities.resolve_special_points_and_path(
        structure, logger, eps=symprec
    )
    if resolved is None:
        return []
    return generate_kpath_parameters(resolved['points'], resolved['path'], 100)


class PhononProperties:
    def __init__(self, phonopy_obj, logger, **kwargs) -> None:
        self.logger = logger
        self.phonopy_obj: Phonopy = phonopy_obj
        self.t_max = kwargs.get('t_max', 1000)
        self.t_min = kwargs.get('t_min', 0)
        self.t_step = kwargs.get('t_step', 100)
        self.band_conf = kwargs.get('band_conf')

        self.n_atoms = len(phonopy_obj.unitcell)

        k_mesh = kwargs.get('k_mesh', 30)
        mesh_density = (2 * k_mesh**3) / self.n_atoms
        mesh_number = int(np.round(mesh_density ** (1.0 / 3.0)))
        self.mesh = [mesh_number, mesh_number, mesh_number]

        self.n_atoms_supercell = len(phonopy_obj.supercell)

    def get_bandstructure(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        phonopy_obj = self.phonopy_obj

        frequency_unit_factor = VaspToTHz
        is_eigenvectors = False

        sym_tol = phonopy_obj.symmetry.tolerance
        if self.band_conf is not None:
            parameters = read_kpath(self.band_conf)
        else:
            parameters = generate_kpath_seekpath(
                phonopy_obj.unitcell, sym_tol, self.logger
            )
        if not parameters:
            return None, None, None

        # Distances calculated in phonopy.band_structure.BandStructure object
        # are based on absolute positions of q-points in reciprocal space
        # as calculated by using the cell which is handed over during instantiation.
        # Fooling that object by handing over a "unit cell" diag(1,1,1) instead clashes
        # with calculation of non-analytical terms.
        # Hence generate appropriate distances and special k-points list based on
        # fractional coordinates in reciprocal space (to keep backwards compatibility
        # with previous FHI-aims phonon implementation).
        bands = []
        bands_distances = []
        distance = 0.0
        bands_special_points = [distance]
        bands_labels = []
        label = parameters[0]['startname']
        for b in parameters:
            kstart = np.array(b['kstart'])
            kend = np.array(b['kend'])
            npoints = b['npoints']
            dk = (kend - kstart) / (npoints - 1)
            bands.append([(kstart + dk * n) for n in range(npoints)])
            dk_length = np.linalg.norm(dk)

            for n in range(npoints):
                bands_distances.append(distance + dk_length * n)

            distance += dk_length * (npoints - 1)
            bands_special_points.append(distance)
            label = [b['startname'], b['endname']]
            bands_labels.append(label)

        bs_obj = BandStructure(
            bands,
            phonopy_obj.dynamical_matrix,
            with_eigenvectors=is_eigenvectors,
            factor=frequency_unit_factor,
        )

        freqs = bs_obj.frequencies

        return np.array(freqs), np.array(bands), np.array(bands_labels)

    def get_dos(self) -> tuple[np.ndarray, np.ndarray]:
        phonopy_obj = self.phonopy_obj
        mesh = self.mesh

        phonopy_obj.run_mesh(mesh, is_gamma_center=True)
        frequencies = phonopy_obj.mesh.frequencies
        self.frequencies = np.array(frequencies)
        min_freq = min(np.ravel(frequencies))
        max_freq = max(np.ravel(frequencies)) + max(np.ravel(frequencies)) * 0.05

        phonopy_obj.run_total_dos(
            freq_min=min_freq,
            freq_max=max_freq,
            use_tetrahedron_method=True,
        )
        dos_result = phonopy_obj.total_dos
        f, dos = dos_result.frequency_points, dos_result.dos

        return f, dos

    def get_thermodynamical_properties(
        self,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        phonopy_obj = self.phonopy_obj

        phonopy_obj.run_mesh(self.mesh, is_gamma_center=True)
        phonopy_obj.run_thermal_properties(
            t_step=self.t_step, t_max=self.t_max, t_min=self.t_min
        )
        thermal_properties = phonopy_obj.thermal_properties
        T = thermal_properties.temperatures
        if hasattr(thermal_properties, 'free_energy'):
            fe = thermal_properties.free_energy
            entropy = thermal_properties.entropy
            cv = thermal_properties.heat_capacity
        else:
            T, fe, entropy, cv = thermal_properties.thermal_properties
        kJmolToEv = 1.0 / EvTokJmol
        fe = fe * kJmolToEv
        JmolToEv = kJmolToEv / 1000
        cv = JmolToEv * cv
        return T, fe, entropy, cv
