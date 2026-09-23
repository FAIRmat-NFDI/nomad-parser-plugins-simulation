from types import SimpleNamespace

import numpy as np
import pytest

import nomad_simulation_parsers.parsers.phonopy.calculator as calculator
from nomad_simulation_parsers.parsers.phonopy.calculator import (
    EvTokJmol,
    PhononProperties,
    generate_kpath_parameters,
    read_kpath,
)

@pytest.mark.unit
def test_generates_kpath_parameters_for_each_segment():
    points = {'G': np.array([0.0, 0.0, 0.0]), 'X': np.array([0.5, 0.0, 0.0])}

    parameters = generate_kpath_parameters(points, [['G', 'X']], 10)

    assert len(parameters) == 1
    assert parameters[0]['npoints'] == 10
    assert parameters[0]['startname'] == 'Γ'
    assert parameters[0]['endname'] == 'X'
    np.testing.assert_allclose(parameters[0]['kstart'], [0.0, 0.0, 0.0])
    np.testing.assert_allclose(parameters[0]['kend'], [0.5, 0.0, 0.0])


@pytest.mark.unit
def test_reads_fractional_band_path_file(tmp_path):
    band_conf = tmp_path / 'band.conf'
    band_conf.write_text(
        'BAND = 0 0 0 1/2 0 0\n'
        'BAND_LABELS = G X\n'
        'BAND_POINTS = 25\n'
    )

    parameters = read_kpath(str(band_conf))

    assert len(parameters) == 1
    assert parameters[0]['npoints'] == 25
    np.testing.assert_allclose(parameters[0]['kend'], [0.5, 0.0, 0.0])


class FakePhonopyObject:
    unitcell = [1, 2, 3, 4]
    supercell = list(range(32))


@pytest.mark.unit
def test_initializes_mesh_from_atom_count():
    properties = PhononProperties(FakePhonopyObject(), logger=None, k_mesh=2)

    assert properties.n_atoms == 4
    assert properties.n_atoms_supercell == 32
    assert properties.mesh == [2, 2, 2]


class FakeCell:
    cell = np.eye(3)

    def __len__(self):
        return 4


class FakeSymmetry:
    tolerance = 1e-5


class FakeBandStructure:
    def __init__(self, bands, dynamical_matrix, with_eigenvectors, factor):
        self.bands = bands
        self.dynamical_matrix = dynamical_matrix
        self.with_eigenvectors = with_eigenvectors
        self.factor = factor
        self.frequencies = np.array([[[1.0], [2.0]]])


class FakePropertiesObject:
    unitcell = FakeCell()
    supercell = [1, 2, 3, 4]
    symmetry = FakeSymmetry()
    dynamical_matrix = object()


@pytest.mark.unit
def test_get_bandstructure_uses_generated_kpath(monkeypatch):
    properties = PhononProperties(FakePropertiesObject(), logger=None)
    parameters = [
        {
            'npoints': 2,
            'startname': 'G',
            'endname': 'X',
            'kstart': np.array([0.0, 0.0, 0.0]),
            'kend': np.array([0.5, 0.0, 0.0]),
        }
    ]
    monkeypatch.setattr(calculator, 'generate_kpath_ase', lambda *args: parameters)
    monkeypatch.setattr(calculator, 'BandStructure', FakeBandStructure)

    frequencies, bands, labels = properties.get_bandstructure()

    assert frequencies.shape == (1, 2, 1)
    assert bands.shape == (1, 2, 3)
    assert labels.tolist() == [['G', 'X']]


class FakeDosResult:
    frequency_points = np.array([-1.0, 0.0, 1.0])
    dos = np.array([0.5, 1.0, 0.5])


class FakeThermalProperties:
    temperatures = np.array([0.0, 300.0])
    free_energy = np.array([2.0, 4.0])
    entropy = np.array([0.0, 1.0])
    heat_capacity = np.array([1.0, 3.0])


class FakeResultsObject(FakePropertiesObject):
    def __init__(self):
        self.supercell = list(range(8))
        self.mesh_calls = []
        self.dos_calls = []
        self.thermal_calls = []
        self.total_dos = FakeDosResult()
        self.thermal_properties = FakeThermalProperties()

    def run_mesh(self, mesh, is_gamma_center):
        self.mesh_calls.append((mesh, is_gamma_center))
        self.mesh = SimpleNamespace(frequencies=np.array([-1.0, 0.0, 1.0]))

    def run_total_dos(self, freq_min, freq_max, use_tetrahedron_method):
        self.dos_calls.append((freq_min, freq_max, use_tetrahedron_method))

    def run_thermal_properties(self, t_step, t_max, t_min):
        self.thermal_calls.append((t_step, t_max, t_min))


@pytest.mark.unit
def test_get_dos_runs_mesh_and_tetrahedron_calculation():
    phonopy_obj = FakeResultsObject()
    properties = PhononProperties(phonopy_obj, logger=None, k_mesh=2)

    frequencies, dos = properties.get_dos()

    np.testing.assert_allclose(frequencies, [-1.0, 0.0, 1.0])
    np.testing.assert_allclose(dos, [0.5, 1.0, 0.5])
    assert phonopy_obj.mesh_calls == [([2, 2, 2], True)]
    assert phonopy_obj.dos_calls == [(-1.0, 1.05, True)]


@pytest.mark.unit
def test_get_thermodynamical_properties_converts_and_scales_values():
    phonopy_obj = FakeResultsObject()
    properties = PhononProperties(
        phonopy_obj, logger=None, k_mesh=2, t_min=100, t_max=500, t_step=200
    )

    temperatures, free_energies, entropy, heat_capacities = (
        properties.get_thermodynamical_properties()
    )

    assert phonopy_obj.mesh_calls == [([2, 2, 2], True)]
    assert phonopy_obj.thermal_calls == [(200, 500, 100)]
    np.testing.assert_allclose(temperatures, [0.0, 300.0])
    np.testing.assert_allclose(entropy, [0.0, 1.0])
    np.testing.assert_allclose(
        free_energies,
        np.array([2.0, 4.0]) / EvTokJmol,
    )
    np.testing.assert_allclose(
        heat_capacities,
        np.array([1.0, 3.0]) / EvTokJmol / 1000,
    )
