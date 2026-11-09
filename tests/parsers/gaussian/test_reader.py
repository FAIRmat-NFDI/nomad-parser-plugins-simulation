import pytest

from nomad_simulation_parsers.parsers.gaussian.text_parser import GaussianOutReader
from tests.parsers.common import approx, assert_approx

HEADER = """ Gaussian 09: EM64L-G09RevD.01 24-Apr-2013 16-May-2014
%NProc=8
"""
ORIENTATION = """ Standard orientation:
 ---------------------------------------------------------------------
 Center     Atomic      Atomic             Coordinates (Angstroms)
 Number     Number       Type             X           Y           Z
 ---------------------------------------------------------------------
      1         13           0        0.000000    0.000000    0.000000
      2          8           0        0.000000    0.000000    1.250000
 ---------------------------------------------------------------------
"""


def read_source(tmp_path, fragment: str) -> dict:
    mainfile = tmp_path / 'calculation.out'
    mainfile.write_text(fragment)
    reader = GaussianOutReader()
    reader.mainfile = str(mainfile)
    return reader.to_dict()


def gaussian_run(body: str, route: str = '#B3LYP/aug-cc-pVTZ') -> str:
    return (
        HEADER
        + f'----------\n {route}\n----------\n'
        + ' Charge = -1 Multiplicity = 1\n'
        + body
        + ' Normal termination of Gaussian 09 at Thu 1-Jan-2026.\n'
    )


@pytest.mark.unit
class TestGaussianReader:
    def test_extracts_program_route_and_spin(self, tmp_path):
        source = read_source(tmp_path, gaussian_run(''))

        assert source['program'] == [
            '09',
            'EM64L-G09RevD.01',
            '24-Apr-2013',
            '16-May-2014',
        ]
        assert source['x_gaussian_number_of_processors'] == '8'
        assert len(source['run']) == 1
        run = source['run'][0]
        assert run['x_gaussian_settings_corrected'] == '#B3LYP/aug-cc-pVTZ'
        assert run['charge'] == -1
        assert run['spin_target'] == 1

    def test_extracts_orientations_and_scf_energy(self, tmp_path):
        source = read_source(
            tmp_path,
            gaussian_run(
                ORIENTATION
                + ' SCF Done:  E(RB3LYP) =  -242.389700911934     A.U.\n'
                + ' Leave Link  502\n'
            ),
        )
        system = source['run'][0]['system'][0]

        assert system['standard_orientation'].shape == (2, 6)
        assert_approx(
            system['standard_orientation'][:, -3:],
            [[0.0, 0.0, 0.0], [0.0, 0.0, 1.25]],
        )
        final_scf = system['scf_iteration_final']
        assert final_scf['energy_total'].to('hartree').magnitude == approx(
            -242.389700911934
        )

    def test_extracts_forces_and_frequencies(self, tmp_path):
        body = (
            ORIENTATION
            + """ Forces (Hartrees/Bohr)
 Number     Number              X              Y              Z
 -------------------------------------------------------------------
      1         13        0.000000000    0.000000000   -0.026613983
      2          8        0.000000000    0.000000000    0.026613983
 --
 Harmonic frequencies (CM**-1)
 Frequencies --   148.4167
 Red. masses --    26.9815
 Temperature   298.150 Kelvin.  Pressure   1.00000 Atm.
"""
        )
        source = read_source(tmp_path, gaussian_run(body))
        system = source['run'][0]['system'][0]

        assert system['forces'].shape == (2, 3)
        assert system['forces'][0, 2] == approx(-0.026613983)
        assert_approx(system['frequencies'], [148.4167])
        assert_approx(system['reduced_masses'], [26.9815])
        assert system['temperature_pressure'][0] == approx(298.15)

    def test_extracts_post_hf_energies(self, tmp_path):
        body = (
            ORIENTATION
            + """ E2 = -3.100000D-01 EUMP2 = -258.83218315785D+00
 EUMP3 = -258.78043052D+00
 UMP4(DQ) = -258.80457235D+00
 UMP4(SDQ) = -258.82670817D+00
 CCSD(T) = -258.84363978D+00
"""
        )
        source = read_source(
            tmp_path, gaussian_run(body, route='#CCSD(T)/6-31G* scf=xqc')
        )
        system = source['run'][0]['system'][0]

        assert len(system['energy_total_mp']) == 4
        assert system['energy_total_mp'][-1].to('hartree').magnitude == approx(
            -258.82670817
        )
        assert system['energy_total_cc'][0].to('hartree').magnitude == approx(
            -258.84363978
        )
