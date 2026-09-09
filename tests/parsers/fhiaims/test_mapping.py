import numpy as np
import pytest
from nomad.units import ureg

from tests.parsers.common import approx
from nomad_simulation_parsers.parsers.fhiaims.parser import FHIAimsOutMappingParser


@pytest.mark.unit
class TestFHIAimsMapping:
    def test_maps_scf_criteria(self):
        parser = FHIAimsOutMappingParser()
        source = {
            'convergence_energy': 1e-6 * ureg.eV,
            'convergence_density': 1e-5,
            'convergence_eigenvalues': 1e-3 * ureg.eV,
            'max_scf_iterations': 100,
        }

        criteria = parser.get_all_criteria(source)

        assert [item['name'] for item in criteria] == [
            'total_energy_change',
            'charge_density_change',
            'sum_eigenvalues_change',
        ]
        assert criteria[0]['threshold_change'].to('eV').magnitude == approx(1e-6)
        assert criteria[1]['n_max_iterations'] == 100

    def test_maps_default_and_explicit_k_offset(self):
        parser = FHIAimsOutMappingParser()
        np.testing.assert_allclose(parser.get_k_offset_with_default(None), [0, 0, 0])
        np.testing.assert_allclose(
            parser.get_k_offset_with_default(np.array([0.5, 0.25, 0])),
            [0.5, 0.25, 0],
        )

    def test_maps_spin_resolved_dos(self, tmp_path):
        parser = FHIAimsOutMappingParser()
        (tmp_path / 'KS_DOS_total.dat').write_text('0 1 2\n1 3 4\n')
        parser.filepath = str(tmp_path / 'aims.out')

        dos = parser.get_dos([['KS_DOS_total.dat']], [], [])

        assert len(dos) == 2
        np.testing.assert_allclose(dos[0]['energies'], [0, 1])
        np.testing.assert_allclose(dos[0]['values'], [1, 3])
        assert dos[1]['spin'] == 1

    def test_maps_gw_workflow_names(self):
        parser = FHIAimsOutMappingParser()
        assert parser.get_gw_flag('gw') == 'G0W0'
        assert parser.get_gw_flag('ev_scgw') == 'ev-scGW'

    def test_maps_unique_kpoints_from_spin_resolved_eigenvalues(self):
        parser = FHIAimsOutMappingParser()
        source = {
            'array_size_parameters': {
                'parameter': [{'Number of spin channels': 2}]
            },
            'geometry_optimization': [
                {
                    'eigenvalues': [
                        {
                            'kpoints': [
                                [0, 0, 0],
                                [0, 0, 0],
                                [0, 0, 0.5],
                                [0, 0, 0.5],
                            ],
                            'occupation_eigenvalue': [
                                [1, 0],
                                [1, 0],
                                [1, 0],
                                [1, 0],
                            ],
                        }
                    ]
                }
            ],
        }

        np.testing.assert_allclose(
            parser.get_kpoints(source), [[0, 0, 0], [0, 0, 0.5]]
        )

    def test_maps_band_gap_from_eigenvalues_and_occupations(self):
        parser = FHIAimsOutMappingParser()
        source = [
            {
                'kpoints': [[0, 0, 0]],
                'occupation_eigenvalue': [
                    [2.0, -1.0],
                    [2.0, -0.5],
                    [0.0, 0.25],
                ],
            }
        ]

        band_gaps = parser.get_band_gaps(source, {})

        assert len(band_gaps) == 1
        assert band_gaps[0]['value'].to('hartree').magnitude == approx(0.75)

    def test_maps_kpoints_in_section_order(self):
        parser = FHIAimsOutMappingParser()
        eigenvalues = {
            'kpoints': [[0, 0, 0]],
            'occupation_eigenvalue': [[1, 0]],
        }
        source = {
            'full_scf': [{'eigenvalues': [eigenvalues]}],
            'geometry_optimization': [
                {
                    'eigenvalues': [
                        {
                            'kpoints': [[0, 0, 0.5]],
                            'occupation_eigenvalue': [[1, 0]],
                        }
                    ]
                }
            ],
        }

        np.testing.assert_allclose(parser.get_kpoints(source), [[0, 0, 0]])
