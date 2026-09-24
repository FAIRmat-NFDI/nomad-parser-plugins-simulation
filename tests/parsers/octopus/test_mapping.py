import numpy as np
import pytest
from nomad.units import ureg

from nomad_simulation_parsers.parsers.octopus.parser import (
    OctopusEigenvalueParser,
    OctopusInfoParser,
    OctopusMainfileParser,
)
from tests.parsers.common import approx, assert_approx


@pytest.mark.unit
class TestOctopusMainfileMapping:
    def test_maps_named_xc_functional(self):
        parser = OctopusMainfileParser()
        parser._info = {'XCFunctional': 'PBE+LDA_C_PZ'}

        assert parser.get_xc_functionals({}) == ['PBE', 'LDA_C_PZ']

    def test_maps_outputs_to_energy_values(self):
        parser = OctopusMainfileParser()
        parser._info = {'energyunit': 'hartree'}
        outputs = parser.get_outputs([{'energy': -1.5}, {'energy_total': -1.25}])

        assert [value['energy'] for value in outputs] == [
            -1.5 * ureg.hartree,
            -1.25 * ureg.hartree,
        ]

    def test_maps_band_gaps_from_occupations(self):
        parser = OctopusEigenvalueParser()
        parser.unit = ureg.hartree
        gaps = parser.get_band_gaps(
            [
                (
                    np.array([0.0, 0.0, 0.0]),
                    np.array([[-2.0], [1.0]]),
                    np.array([[1.0], [0.0]]),
                )
            ]
        )

        assert gaps[0]['value'].to('hartree').magnitude == approx(3.0)


@pytest.mark.unit
class TestOctopusInfoMapping:
    def test_maps_eigenvalues_occupations_and_band_structures(self):
        parser = OctopusInfoParser()
        parser.unit = ureg.eV
        source = [
            (
                np.array([0.0, 0.0, 0.0]),
                np.array([[-2.0], [1.0]]),
                np.array([[1.0], [0.0]]),
            )
        ]

        eigenvalues = parser.get_eigenvalues(source)
        band_structures = parser.get_band_structures(source)

        assert len(eigenvalues) == 1
        assert_approx(
            eigenvalues[0]['eigenvalues'].to('eV').magnitude,
            [[-2.0, 1.0]],
        )
        assert_approx(eigenvalues[0]['occupations'], [[1.0, 0.0]])
        assert_approx(eigenvalues[0]['kpoints'], [[0.0, 0.0, 0.0]])
        assert_approx(
            band_structures[0]['value'].to('eV').magnitude,
            [[-2.0, 1.0]],
        )
