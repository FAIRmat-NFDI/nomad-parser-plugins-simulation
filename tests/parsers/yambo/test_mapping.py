from datetime import datetime

import numpy as np
import pytest
from nomad.units import ureg

from nomad_simulation_parsers.parsers.yambo.parser import (
    YamboMainfileParser,
    YamboNetCDFParser,
)
from tests.parsers.common import assert_approx


@pytest.mark.unit
class TestYamboMainfileMapping:
    @pytest.fixture(scope='class')
    def parser(self):
        return YamboMainfileParser()

    def test_maps_eigenvalues_and_reference_energy(self, parser):
        outputs = parser.get_outputs(
            {
                'eigenenergies': {
                    'kpoints': np.array([[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]]),
                    'energies': np.array([1.0, 2.0, 1.5, 2.5]),
                },
                'valence_conduction': [1.2, 2.3],
            },
            [],
            {},
        )

        assert len(outputs) == 1
        assert outputs[0]['eigenvalues'][0]['energies'].shape == (2, 2)
        assert outputs[0]['eigenvalues'][0][
            'highest_occupied'
        ].to('eV').magnitude == pytest.approx(1.2)

    def test_maps_band_gap_from_valence_and_conduction(self, parser):
        result = parser.get_band_gaps(None, 1.0 * ureg.eV, 3.5 * ureg.eV)

        assert len(result) == 1
        assert result[0]['value'].to('eV').magnitude == pytest.approx(2.5)

    def test_maps_program_name_and_start_time(self, parser):
        assert parser.get_wallstart('29/11/2021 20:57') == datetime.strptime(
            '29/11/2021 20:57', '%d/%m/%Y %H:%M'
        ).timestamp()


@pytest.mark.unit
class TestYamboNetCDFMapping:
    @pytest.fixture(scope='class')
    def parser(self):
        parser = YamboNetCDFParser.__new__(YamboNetCDFParser)
        parser._data = {
            'MAX_ATOMS': np.array([2]),
            'N_ATOMS': np.array([1, 2]),
            'ATOM_POS': np.arange(12, dtype=float),
            'atomic_numbers': np.array([6, 1]),
            'LATTICE_VECTORS': np.eye(3),
            'K-POINTS': np.array([[0.0, 0.5], [0.0, 0.0], [0.0, 0.0]]),
        }
        return parser

    def test_maps_structure_and_kspace_helpers(self, parser):
        assert_approx(parser.get_positions(), [[0, 1, 2], [6, 7, 8], [9, 10, 11]])
        assert parser.get_labels() == ['C', 'H', 'H']
        assert_approx(parser.get_lattice_vectors(), np.eye(3))
        assert_approx(parser.get_kpoints(), [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]])

    def test_maps_eigenvalues(self, parser):
        parser._data['EIGENVALUES'] = np.array([[1.0, 2.0], [3.0, 4.0]])

        eigenvalues = parser.get_eigenvalues()

        assert len(eigenvalues) == 2
        assert_approx(eigenvalues[0]['energies'].to('eV').magnitude, [1.0, 2.0])
        assert_approx(eigenvalues[1]['energies'].to('eV').magnitude, [3.0, 4.0])
