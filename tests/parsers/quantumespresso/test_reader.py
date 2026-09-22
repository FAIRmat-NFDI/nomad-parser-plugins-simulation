from pathlib import Path

import pytest

from nomad_simulation_parsers.parsers.quantumespresso.epw.file_parser import (
    EPWFileParser,
)
from nomad_simulation_parsers.parsers.quantumespresso.gipaw.file_parser import (
    GIPAWFileParser,
)
from nomad_simulation_parsers.parsers.quantumespresso.phonon.file_parser import (
    PhononFileParser,
)
from nomad_simulation_parsers.parsers.quantumespresso.pwscf.file_parser import (
    PWSCFFileParser,
)
from nomad_simulation_parsers.parsers.quantumespresso.xspectra.file_parser import (
    XSpectraFileParser,
)

from nomad_simulation_parsers.parsers.quantumespresso.file_parser import (
    QuantumEspressoFileParser,
)
from tests.parsers.common import approx, assert_approx


DATA_DIR = Path(__file__).resolve().parents[2] / 'data' / 'quantumespresso'


def read(parser_class, path):
    parser = parser_class()
    parser.mainfile = str(DATA_DIR / path)
    return parser.to_dict()


@pytest.mark.unit
class TestQuantumEspressoFileReader:
    def test_reads_program_blocks(self, tmp_path):
        mainfile = tmp_path / 'pw.out'
        mainfile.write_text(
            'Program PWSCF v.7.3 starts\n'
            'some output\n'
            'JOB DONE.\n'
        )

        parser = QuantumEspressoFileParser()
        parser.mainfile = str(mainfile)
        source = parser.to_dict()

        assert len(source['program']) == 1
        assert source['program'][0]['header'].startswith('Program PWSCF v.7.3')


@pytest.mark.unit
class TestEPWReader:
    def test_reads_epw_header_and_meshes(self):
        source = read(EPWFileParser, 'epw/epw.out')

        assert source['header']['program_name_version'] == ['EPW', 'v.5.4.1']
        assert source['header']['number_of_atoms'] == 12
        assert_approx(source['q_mesh'], [7, 7, 2])
        assert_approx(source['k_mesh'], [14, 14, 4])
        assert source['n_q_mesh'] == 98
        assert source['n_k_mesh'] == 1568
        assert source['e_fermi_coarse_grid'].to('eV').magnitude == approx(12.256609)


@pytest.mark.unit
class TestGIPAWReader:
    def test_reads_gipaw_chemical_shifts_and_tensors(self):
        source = read(
            GIPAWFileParser,
            'gipaw/scf_out_nmr_out_741/quartz-nmr.out',
        )

        assert source['header']['program_name_version'] == ['GIPAW', 'v.7.4.1']
        assert len(source['ms_list']) == 9
        assert source['ms_list'][0][:2] == ['Si', 1]
        assert_approx(
            source['ms_list'][0][2:],
            [
                432.5142,
                0.5825,
                6.1764,
                0.5819,
                431.8332,
                3.5599,
                -11.3083,
                -6.5295,
                429.5852,
            ],
        )
        assert_approx(
            source['chi_bare_pGv'],
            [
                [-6.48376e-11, 0, 0],
                [0, -6.48840e-11, -4.95e-14],
                [0, 5.6e-14, -6.51570e-11],
            ],
        )


@pytest.mark.unit
class TestPhononReader:
    def test_reads_phonon_calculations(self):
        source = read(PhononFileParser, 'phonon/ph.out')

        assert source['header']['program_name_version'] == ['PHONON', 'v.6.8']
        assert source['header']['nproc'] == 24
        assert len(source['calculation']) == 40
        calculation = source['calculation'][0]
        assert calculation['number_of_atoms'] == 2
        assert calculation['number_of_species'] == 2
        assert calculation['xc_functional'] == 'PBE'
        assert_approx(
            calculation['alat'].to('angstrom').magnitude,
            2.923545335,
        )


@pytest.mark.unit
class TestPWSCFReader:
    def test_reads_pwscf_scf_and_geometry_optimization(self):
        source = read(PWSCFFileParser, 'pwscf/TiO2_opt/pw.out')

        assert source['header']['program_name_version'] == ['PWSCF', 'v.7.3']
        assert source['header']['number_of_atoms'] == 6
        assert source['header']['number_of_species'] == 2
        assert len(source['self_consistent']['iteration']) == 13
        assert source['bfgs_geometry_optimization']['convergence'].tolist() == [5, 4]
        energy = source['self_consistent']['iteration'][0]['energies']['energy_total']
        assert energy.to('eV').magnitude == approx(-5513.223)


@pytest.mark.unit
class TestXSpectraReader:
    def test_reads_xspectra_header_and_xanes_steps(self):
        source = read(
            XSpectraFileParser,
            'xspectra/ms-10734/Spectra-1-1-1/0/dipole1/xanes.out',
        )

        assert source['header']['program_name_version'] == ['XSpectra', 'v.6.7MaX']
        assert source['header']['number_of_atoms'] == 72
        assert source['xanes']['algorithm'] == ['Lanczos', 'recursion']
        assert source['xanes']['step_1']['k_calculation'][0]['converged'] is True
        assert source['xanes']['step_2']['xnepoint'] == 400
        assert (
            source['xanes']['step_2']['energy_zero'].to('eV').magnitude
            == approx(15.157)
        )
