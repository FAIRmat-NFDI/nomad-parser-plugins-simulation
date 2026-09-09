import gzip

import numpy as np
import pytest
from nomad.units import ureg

from nomad_simulation_parsers.parsers.abinit.file_parser import AbinitOutParser
from nomad_simulation_parsers.parsers.abinit.parser import DosParser

MINIMAL_OUTPUT = """
.Version 9.10.4 of ABINIT
.Starting date : Mon 01 Jan 2024.
- ( at 12h30 )
-outvars: echo values of preprocessed input variables -----
 natom 2
 nsppol 1
 ixc 11
 ionmov 2
==========
== DATASET 1 ==
========== iter Etot
 ETOT 1 -1.000000 1.0E-2
 ETOT 2 -1.100000 1.0E-3
==========
== END DATASET
"""

MINIMAL_RESULTS_OUTPUT = """
.Version 9.10.4 of ABINIT
- Total cpu time (s,m,h): 1.4
- Total wall clock time (s,m,h): 0.5
Calculation completed
== DATASET 1 ==
----------iterations are completed or convergence reached----------
 >>>>>>>>> Etotal = -8.866223895975928
 Kinetic energy = 3.06072461520598
 Fermi (or HOMO) energy (hartree) = -0.1
 Cartesian components of stress tensor (hartree/bohr^3)
 sigma(1 1) = 1.0
 sigma(1 2) = 2.0
 sigma(2 2) = 3.0
 sigma(1 3) = 4.0
 sigma(3 3) = 5.0
 sigma(2 3) = 6.0
================================================================================
== END DATASET
"""


def read_source(path):
    parser = AbinitOutParser()
    parser.mainfile = str(path)
    return parser.to_dict()


@pytest.mark.unit
class TestAbinitReader:
    def test_extracts_header_inputs_and_scf_rows(self, tmp_path):
        mainfile = tmp_path / 'minimal.out'
        mainfile.write_text(MINIMAL_OUTPUT)

        source = read_source(mainfile)

        assert source['program_version'] == '9.10.4'
        assert source['x_abinit_start_date'] == 'Mon 01 Jan 2024'
        assert source['x_abinit_start_time'] == '12h30'
        assert source['input_variables']['key_value'] == [
            ['natom', 2],
            ['nsppol', 1],
            ['ixc', 11],
            ['ionmov', 2],
        ]
        assert source['dataset'][0]['x_abinit_dataset_number'] == 1
        np.testing.assert_allclose(
            source['dataset'][0]['self_consistent']['energy_total_scf_iteration'],
            [[-1.0, 1.0e-2], [-1.1, 1.0e-3]],
        )


    def test_compression_does_not_change_extracted_values(self, tmp_path):
        plain = tmp_path / 'minimal.out'
        compressed = tmp_path / 'minimal.out.gz'
        plain.write_text(MINIMAL_OUTPUT)
        with gzip.open(compressed, 'wt') as handle:
            handle.write(MINIMAL_OUTPUT)

        plain_source = read_source(plain)
        compressed_source = read_source(compressed)

        assert compressed_source['program_version'] == plain_source['program_version']
        np.testing.assert_allclose(
            compressed_source['dataset'][0]['self_consistent'][
                'energy_total_scf_iteration'
            ],
            plain_source['dataset'][0]['self_consistent'][
                'energy_total_scf_iteration'
            ],
        )


    def test_truncated_output_is_returned_as_partial_source(self, tmp_path):
        mainfile = tmp_path / 'truncated.out'
        mainfile.write_text('.Version 9.10.4 of ABINIT\n== DATASET 1 ==\n')

        source = read_source(mainfile)

        assert source['program_version'] == '9.10.4'
        assert source.get('dataset') is None


    def test_extracts_eigenvalues_and_occupations_from_results_block(self, tmp_path):
        mainfile = tmp_path / 'bands.out'
        mainfile.write_text(
            """
.Version 9.10.4 of ABINIT
== DATASET 1 ==
 ----iterations are completed or convergence reached----
 kpt# 1, nband= 2, wtk= 1.00000, kpt= 0.0 0.0 0.0 (reduced coord)
 -1.0 1.0
 occupation numbers for kpt# 1
 2.0 0.0
================================================================================
== END DATASET
"""
        )

        results = read_source(mainfile)['dataset'][0]['results']

        np.testing.assert_allclose(
            results['eigenvalues'], [[1.0, 2.0, 1.0, 0.0, 0.0, 0.0, -1.0, 1.0]]
        )
        np.testing.assert_allclose(results['occupation_numbers'], [[2.0, 0.0]])


    def test_extracts_numeric_dos_table(self, tmp_path):
        dos_file = tmp_path / 'calculation_o_DS2_DOS'
        dos_file.write_text(
            '# energy DOS integrated-DOS\n-1.0 10.0 0.0\n0.0 20.0 1.0\n'
        )
        parser = DosParser()
        parser.filepath = str(dos_file)

        np.testing.assert_allclose(
            parser.data['data'], [[-1.0, 10.0, 0.0], [0.0, 20.0, 1.0]]
        )


    def test_extracts_runtime_and_results(self, tmp_path):
        mainfile = tmp_path / 'minimal-results.out'
        mainfile.write_text(MINIMAL_RESULTS_OUTPUT)
        source = read_source(mainfile)

        assert source['x_abinit_total_cpu_time'] == pytest.approx(1.4)
        assert source['x_abinit_total_wallclock_time'] == pytest.approx(0.5)
        assert source['run_clean_end'] == 'Calculation completed'

        results = source['dataset'][0]['results']
        assert results['energy_total'].units == ureg.hartree
        assert results['fermi_energy'].units == ureg.hartree
        assert results['stress_tensor'].units == ureg.hartree / ureg.bohr**3
        assert results['energy_kinetic_electronic'].units == ureg.hartree
        assert results['energy_total'].magnitude == pytest.approx(-8.866223895975928)
