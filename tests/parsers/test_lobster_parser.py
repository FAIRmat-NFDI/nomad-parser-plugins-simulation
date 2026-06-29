import ase
import numpy as np
import pytest
from nomad.datamodel import EntryArchive
from nomad.utils import get_logger
from packaging.version import Version

from nomad_simulation_parsers.parsers.lobster.parser import LobsterParser

logger = get_logger(__name__)


@pytest.fixture
def parser():
    return LobsterParser()


def approx(value):
    return pytest.approx(value, abs=0, rel=1e-6)


def test_Fe(parser):  # noqa: PLR0915
    """
    Tests spin-polarized Fe calculation with LOBSTER 4.0.0
    """

    archive = EntryArchive()
    parser.parse('tests/data/lobster/Fe/lobsterout', archive, logger)

    data = archive.data
    assert data.program.version == '4.0.0'
    assert data.wall_start.magnitude == 1619680785.0

    assert len(data.model_system) == 1
    model_system = data.model_system[0]
    assert model_system.positions.to('angstrom').magnitude[0][1] == approx(1.41588779)
    assert len(model_system.particle_states) == 2
    assert model_system.particle_states[0].chemical_symbol == 'Fe'
    assert model_system.lattice_vectors.to('angstrom').magnitude[0][0] == approx(
        2.83177559
    )
    assert model_system.periodic_boundary_conditions == [True, True, True]

    assert len(data.model_method) == 1
    assert len(data.model_method[0].numerical_settings) == 1
    assert len(data.model_method[0].numerical_settings[0].basis_set_components) == 1
    assert (
        data.model_method[0].numerical_settings[0].basis_set_components[0].basis_set
        == 'pbeVaspFit2015'
    )

    assert len(data.outputs) == 1
    output = data.outputs[0]
    # ICOHPLIST.lobster
    cohp = output.x_lobster_section_cohp
    assert cohp.x_lobster_number_of_cohp_pairs == 20

    assert cohp.x_lobster_number_of_cohp_pairs == 20
    assert len(cohp.x_lobster_cohp_atom1_labels) == 20
    assert cohp.x_lobster_cohp_atom1_labels[19] == 'Fe2'
    assert len(cohp.x_lobster_cohp_atom2_labels) == 20
    assert cohp.x_lobster_cohp_atom1_labels[3] == 'Fe1'
    assert len(cohp.x_lobster_cohp_distances) == 20
    assert cohp.x_lobster_cohp_distances[0].to('angstrom').magnitude == approx(2.83178)
    assert cohp.x_lobster_cohp_distances[13].to('angstrom').magnitude == approx(2.45239)
    assert cohp.x_lobster_cohp_distances[19].to('angstrom').magnitude == approx(2.83178)
    assert np.shape(cohp.x_lobster_cohp_translations) == (20, 3)
    assert (cohp.x_lobster_cohp_translations[0] == [0, 0, -1]).all()
    assert (cohp.x_lobster_cohp_translations[13] == [0, 0, 0]).all()
    assert (cohp.x_lobster_cohp_translations[19] == [0, 0, 1]).all()
    assert np.shape(cohp.x_lobster_integrated_cohp_at_fermi_level) == (2, 20)
    assert cohp.x_lobster_integrated_cohp_at_fermi_level[0, 0].to(
        'eV'
    ).magnitude == approx(-0.08672)
    assert cohp.x_lobster_integrated_cohp_at_fermi_level[0, 19].to(
        'eV'
    ).magnitude == approx(-0.08672)
    assert cohp.x_lobster_integrated_cohp_at_fermi_level[1, 19].to(
        'eV'
    ).magnitude == approx(-0.16529)
    assert cohp.x_lobster_integrated_cohp_at_fermi_level[1, 7].to(
        'eV'
    ).magnitude == approx(-0.48790)

    # # COHPCAR.lobster
    assert len(cohp.x_lobster_cohp_energies) == 201
    assert cohp.x_lobster_cohp_energies[0].to('eV').magnitude == approx(-10.06030)
    assert cohp.x_lobster_cohp_energies[200].to('eV').magnitude == approx(3.00503)
    assert np.shape(cohp.x_lobster_average_cohp_values) == (2, 201)
    assert cohp.x_lobster_average_cohp_values[0][196] == approx(0.02406)
    assert cohp.x_lobster_average_cohp_values[1][200] == approx(0.014390)
    assert np.shape(cohp.x_lobster_average_integrated_cohp_values) == (2, 201)
    assert cohp.x_lobster_average_integrated_cohp_values[0][200].to(
        'eV'
    ).magnitude == approx(-0.06616)
    assert cohp.x_lobster_average_integrated_cohp_values[1][200].to(
        'eV'
    ).magnitude == approx(-0.11366)
    assert np.shape(cohp.x_lobster_cohp_values) == (20, 2, 201)
    assert cohp.x_lobster_cohp_values[10][1][200] == approx(0.02291)
    assert cohp.x_lobster_cohp_values[19][0][200] == approx(0.01816)
    assert np.shape(cohp.x_lobster_integrated_cohp_values) == (20, 2, 201)
    assert cohp.x_lobster_integrated_cohp_values[10][0][200].to(
        'eV'
    ).magnitude == approx(-0.11401)
    assert cohp.x_lobster_integrated_cohp_values[19][1][200].to(
        'eV'
    ).magnitude == approx(-0.06876)

    # ICOOPLIST.lobster
    coop = output.x_lobster_section_coop
    assert coop.x_lobster_number_of_coop_pairs == 20
    assert len(coop.x_lobster_coop_atom1_labels) == 20
    assert coop.x_lobster_coop_atom1_labels[19] == 'Fe2'
    assert len(coop.x_lobster_coop_atom2_labels) == 20
    assert coop.x_lobster_coop_atom1_labels[3] == 'Fe1'
    assert len(coop.x_lobster_coop_distances) == 20
    assert coop.x_lobster_coop_distances[0].to('angstrom').magnitude == approx(2.83178)
    assert coop.x_lobster_coop_distances[13].to('angstrom').magnitude == approx(2.45239)
    assert coop.x_lobster_coop_distances[19].to('angstrom').magnitude == approx(2.83178)
    assert np.shape(coop.x_lobster_coop_translations) == (20, 3)
    assert (coop.x_lobster_coop_translations[0] == [0, 0, -1]).all()
    assert (coop.x_lobster_coop_translations[13] == [0, 0, 0]).all()
    assert (coop.x_lobster_coop_translations[19] == [0, 0, 1]).all()
    assert np.shape(coop.x_lobster_integrated_coop_at_fermi_level) == (2, 20)
    assert coop.x_lobster_integrated_coop_at_fermi_level[0, 0].magnitude == approx(
        -0.06882
    )
    assert coop.x_lobster_integrated_coop_at_fermi_level[0, 19].magnitude == approx(
        -0.06882
    )
    assert coop.x_lobster_integrated_coop_at_fermi_level[1, 19].magnitude == approx(
        -0.11268
    )
    assert coop.x_lobster_integrated_coop_at_fermi_level[1, 7].magnitude == approx(
        -0.05179
    )

    # COOPCAR.lobster
    assert len(coop.x_lobster_coop_energies) == 201
    assert coop.x_lobster_coop_energies[0].to('eV').magnitude == approx(-10.06030)
    assert coop.x_lobster_coop_energies[200].to('eV').magnitude == approx(3.00503)
    assert np.shape(coop.x_lobster_average_coop_values) == (2, 201)
    assert coop.x_lobster_average_coop_values[0][196] == approx(-0.04773)
    assert coop.x_lobster_average_coop_values[1][200] == approx(-0.00788)
    assert np.shape(coop.x_lobster_average_integrated_coop_values) == (2, 201)
    assert coop.x_lobster_average_integrated_coop_values[0][200].magnitude == approx(
        -0.12265
    )
    assert coop.x_lobster_average_integrated_coop_values[1][200].magnitude == approx(
        -0.10557
    )
    assert np.shape(coop.x_lobster_coop_values) == (20, 2, 201)
    assert coop.x_lobster_coop_values[3][1][200] == approx(-0.01346)
    assert coop.x_lobster_coop_values[0][0][200] == approx(-0.04542)
    assert np.shape(coop.x_lobster_integrated_coop_values) == (20, 2, 201)
    assert coop.x_lobster_integrated_coop_values[10][0][199].magnitude == approx(
        -0.11299
    )
    assert coop.x_lobster_integrated_coop_values[19][1][200].magnitude == approx(
        -0.13041
    )

    # CHARGE.lobster
    charges = output.x_lobster_section_charges
    assert len(charges) == 2
    mulliken = charges[0]
    assert mulliken.type == 'mulliken'
    assert len(mulliken.contributions) == 2
    assert mulliken.contributions[0].value.to(
        'elementary_charge'
    ).magnitude == pytest.approx(0.0, abs=1e-6)
    assert mulliken.contributions[1].value.to(
        'elementary_charge'
    ).magnitude == pytest.approx(0.0, abs=1e-6)

    loewdin = charges[1]
    assert loewdin.type == 'loewdin'
    assert len(loewdin.contributions) == 2
    assert loewdin.contributions[0].value.to(
        'elementary_charge'
    ).magnitude == pytest.approx(0.0, abs=1e-6)
    assert loewdin.contributions[1].value.to(
        'elementary_charge'
    ).magnitude == pytest.approx(0.0, abs=1e-6)

    # DOSCAR.lobster total and integrated DOS
    assert len(output.electronic_dos) == 2
    dos_up = output.electronic_dos[0]
    dos_down = output.electronic_dos[1]
    # assert dos_up.n_energies == 201
    assert len(dos_up.energies.points) == 201
    assert dos_up.energies.points[0].to('eV').magnitude == approx(-10.06030)
    assert dos_up.energies.points[16].to('eV').magnitude == approx(-9.01508)
    assert dos_up.energies.points[200].to('eV').magnitude == approx(3.00503)
    assert len(dos_up.value) == len(dos_down.value) == 201
    assert dos_up.value[6].magnitude == pytest.approx(0.0, abs=1e-30)
    assert dos_up.value[200].to('1/eV').magnitude == approx(0.26779)
    assert dos_down.value[195].to('1/eV').magnitude == approx(0.37457)
    # TODO implement value_integrated in schema
    # assert np.shape(dos_up.total[0].value_integrated) == (201,)
    # assert dos_up.total[0].value_integrated[10] == approx(0.0 + 18)
    # assert dos_up.total[0].value_integrated[188] == approx(11.07792 + 18)
    # assert dos_down.total[0].value_integrated[200] == approx(10.75031 + 18)

    # DOSCAR.lobster atom and lm-projected dos
    assert len(dos_up.projected_dos) == 12 and len(dos_down.projected_dos) == 12
    # assert (
    #     dos_up.atom_projected[0].atom_index == 0
    #     and dos_up.atom_projected[6].atom_index == 1
    # )
    # assert dos_up.atom_projected[0].m_kind == 'real_orbital'
    # assert (dos_up.atom_projected[4].lm == [2, 1]).all()
    assert np.shape(dos_up.projected_dos[11].value) == (201,)
    assert dos_up.projected_dos[5].value[190].to('1/eV').magnitude == approx(0.00909)
    assert dos_down.projected_dos[5].value[190].to('1/eV').magnitude == approx(0.29205)


def test_NaCl(parser):  # noqa: PLR0915
    """
    Test non-spin-polarized NaCl calculation with LOBSTER 3.2.0
    """

    archive = EntryArchive()
    parser.parse('tests/data/lobster/NaCl/lobsterout', archive, logger)

    data = archive.data
    assert data.program.name == 'LOBSTER'
    assert data.program.version == '3.2.0'
    assert data.wall_start.magnitude == 1619705848.0

    method = data.model_method
    assert len(method) == 1
    assert method[0].x_lobster_code == 'VASP'
    assert len(method[0].numerical_settings) == 1
    assert len(method[0].numerical_settings[0].basis_set_components) == 2
    assert (
        method[0].numerical_settings[0].basis_set_components[0].basis_set
        == 'pbeVaspFit2015'
    )

    assert len(data.outputs) == 1
    output = data.outputs[0]
    assert len(output.x_lobster_abs_total_spilling) == 1
    assert output.x_lobster_abs_total_spilling[0] == approx(9.29)
    assert len(output.x_lobster_abs_charge_spilling) == 1
    assert output.x_lobster_abs_charge_spilling[0] == approx(0.58)

    # ICOHPLIST.lobster
    cohp = output.x_lobster_section_cohp
    assert cohp.x_lobster_number_of_cohp_pairs == 72
    assert len(cohp.x_lobster_cohp_atom1_labels) == 72
    assert cohp.x_lobster_cohp_atom1_labels[71] == 'Cl7'
    assert len(cohp.x_lobster_cohp_atom2_labels) == 72
    assert cohp.x_lobster_cohp_atom2_labels[43] == 'Cl6'
    assert len(cohp.x_lobster_cohp_distances) == 72
    assert cohp.x_lobster_cohp_distances[0].to('angstrom').magnitude == approx(3.99586)
    assert cohp.x_lobster_cohp_distances[47].to('angstrom').magnitude == approx(2.82550)
    assert cohp.x_lobster_cohp_distances[71].to('angstrom').magnitude == approx(3.99586)
    assert np.shape(cohp.x_lobster_cohp_translations) == (72, 3)
    assert (cohp.x_lobster_cohp_translations[0] == [-1, 0, 0]).all()
    assert (cohp.x_lobster_cohp_translations[54] == [0, -1, 0]).all()
    assert (cohp.x_lobster_cohp_translations[71] == [0, 1, 0]).all()
    assert np.shape(cohp.x_lobster_integrated_cohp_at_fermi_level) == (1, 72)
    assert cohp.x_lobster_integrated_cohp_at_fermi_level[0, 0].to(
        'eV'
    ).magnitude == approx(-0.02652)
    assert cohp.x_lobster_integrated_cohp_at_fermi_level[0, 71].to(
        'eV'
    ).magnitude == approx(-0.02925)

    # COHPCAR.lobster
    assert len(cohp.x_lobster_cohp_energies) == 201
    assert cohp.x_lobster_cohp_energies[0].to('eV').magnitude == approx(-12.02261)
    assert cohp.x_lobster_cohp_energies[200].to('eV').magnitude == approx(2.55025)
    assert np.shape(cohp.x_lobster_average_cohp_values) == (1, 201)
    assert cohp.x_lobster_average_cohp_values[0][0] == pytest.approx(0.0)
    assert cohp.x_lobster_average_cohp_values[0][151] == approx(-0.03162)
    assert np.shape(cohp.x_lobster_average_integrated_cohp_values) == (1, 201)
    assert cohp.x_lobster_average_integrated_cohp_values[0][0].to(
        'eV'
    ).magnitude == approx(-0.15834)
    assert cohp.x_lobster_average_integrated_cohp_values[0][200].to(
        'eV'
    ).magnitude == approx(-0.24310)
    assert np.shape(cohp.x_lobster_cohp_values) == (72, 1, 201)
    assert cohp.x_lobster_cohp_values[1][0][200] == pytest.approx(0.0)
    assert cohp.x_lobster_cohp_values[71][0][140] == approx(-0.00403)
    assert np.shape(cohp.x_lobster_integrated_cohp_values) == (72, 1, 201)
    assert cohp.x_lobster_integrated_cohp_values[2][0][200].to(
        'eV'
    ).magnitude == approx(-0.02652)
    assert cohp.x_lobster_integrated_cohp_values[67][0][199].to(
        'eV'
    ).magnitude == approx(-0.04137)

    # ICOOPLIST.lobster
    coop = output.x_lobster_section_coop
    assert coop.x_lobster_number_of_coop_pairs == 72
    assert len(coop.x_lobster_coop_atom1_labels) == 72
    assert coop.x_lobster_coop_atom1_labels[71] == 'Cl7'
    assert len(coop.x_lobster_coop_atom2_labels) == 72
    assert coop.x_lobster_coop_atom2_labels[0] == 'Na2'
    assert len(coop.x_lobster_coop_distances) == 72
    assert coop.x_lobster_coop_distances[0].to('angstrom').magnitude == approx(3.99586)
    assert coop.x_lobster_coop_distances[12].to('angstrom').magnitude == approx(2.82550)
    assert coop.x_lobster_coop_distances[71].to('angstrom').magnitude == approx(3.99586)
    assert np.shape(coop.x_lobster_coop_translations) == (72, 3)
    assert (coop.x_lobster_coop_translations[0] == [-1, 0, 0]).all()
    assert (coop.x_lobster_coop_translations[13] == [0, 1, 0]).all()
    assert (coop.x_lobster_coop_translations[71] == [0, 1, 0]).all()
    assert np.shape(coop.x_lobster_integrated_coop_at_fermi_level) == (1, 72)
    assert coop.x_lobster_integrated_coop_at_fermi_level[0, 0].magnitude == approx(
        -0.00519
    )
    assert coop.x_lobster_integrated_coop_at_fermi_level[0, 71].magnitude == approx(
        -0.00580
    )

    # COOPCAR.lobster
    assert len(coop.x_lobster_coop_energies) == 201
    assert coop.x_lobster_coop_energies[0].to('eV').magnitude == approx(-12.02261)
    assert coop.x_lobster_coop_energies[200].to('eV').magnitude == approx(2.55025)
    assert np.shape(coop.x_lobster_average_coop_values) == (1, 201)
    assert coop.x_lobster_average_coop_values[0][0] == pytest.approx(0.0)
    assert coop.x_lobster_average_coop_values[0][145] == approx(0.03178)
    assert np.shape(coop.x_lobster_average_integrated_coop_values) == (1, 201)
    assert coop.x_lobster_average_integrated_coop_values[0][0].magnitude == approx(
        0.00368
    )
    assert coop.x_lobster_average_integrated_coop_values[0][200].magnitude == approx(
        0.00682
    )
    assert np.shape(coop.x_lobster_coop_values) == (72, 1, 201)
    assert coop.x_lobster_coop_values[1][0][200] == pytest.approx(0.0)
    assert coop.x_lobster_coop_values[71][0][143] == approx(0.01862)
    assert np.shape(coop.x_lobster_integrated_coop_values) == (72, 1, 201)
    assert coop.x_lobster_integrated_coop_values[2][0][200].magnitude == approx(
        -0.00519
    )
    assert coop.x_lobster_integrated_coop_values[71][0][199].magnitude == approx(
        -0.00580
    )

    # CHARGE.lobster
    charges = output.x_lobster_section_charges
    assert len(charges) == 2
    mulliken = charges[0]
    assert mulliken.type == 'mulliken'
    # here the approx is not really working (changing the 0.78 to for example
    # 10 makes the test still pass)
    assert mulliken.contributions[0].value.to('elementary_charge').magnitude == approx(
        0.78
    )
    assert mulliken.contributions[7].value.to('elementary_charge').magnitude == approx(
        -0.78
    )

    loewdin = charges[1]
    assert loewdin.type == 'loewdin'
    assert loewdin.contributions[0].value.to('elementary_charge').magnitude == approx(
        0.67
    )
    assert loewdin.contributions[7].value.to('elementary_charge').magnitude == approx(
        -0.67
    )

    # DOSCAR.lobster total and integrated DOS
    assert len(output.electronic_dos) == 1
    dos = output.electronic_dos[0]
    # assert dos.n_energies == 201
    assert len(dos.energies.points) == 201
    assert dos.energies.points[0].to('eV').magnitude == approx(-12.02261)
    assert dos.energies.points[25].to('eV').magnitude == approx(-10.20101)
    assert dos.energies.points[200].to('eV').magnitude == approx(2.55025)
    assert np.shape(dos.value) == (201,)
    assert dos.value[6].to('1/eV').magnitude == pytest.approx(0.0, abs=1e-30)
    assert dos.value[162].to('1/eV').magnitude == approx(20.24722)
    assert dos.value[200].to('1/eV').magnitude == pytest.approx(0.0, abs=1e-30)
    # assert np.shape(dos.value_integrated) == (201,)
    # assert dos.value_integrated[10] == approx(7.99998 + 80)
    # assert dos.value_integrated[160] == approx(27.09225 + 80)
    # assert dos.value_integrated[200] == approx(31.99992 + 80)

    # DOSCAR.lobster atom and lm-projected dos
    assert len(dos.projected_dos) == 20
    # dos.projected_dos[0].atom_index == 0
    # dos.projected_dos[19].atom_index == 7
    # assert dos.projected_dos[5].m_kind == 'real_orbital'
    # assert (dos.projected_dos[17].lm == [1, 2]).all()
    assert np.shape(dos.projected_dos[13].value) == (201,)
    assert np.shape(dos.projected_dos[8].value) == (201,)
    assert dos.projected_dos[0].value[190].to('1/eV').magnitude == pytest.approx(
        0.0, abs=1e-30
    )
    assert dos.projected_dos[19].value[141].to('1/eV').magnitude == approx(0.32251)
    assert dos.projected_dos[16].value[152].to('1/eV').magnitude == approx(0.00337)


def test_HfV(parser):  # noqa: PLR0915
    """
    Test non-spin-polarized HfV2 calculation with LOBSTER 2.0.0,
    it has different ICOHPLIST.lobster and ICOOPLIST.lobster scheme.
    Also test backup structure parsing when no CONTCAR is present.
    """

    archive = EntryArchive()
    parser.parse('tests/data/lobster/HfV2/lobsterout', archive, logger)

    data = archive.data
    assert data.program.name == 'LOBSTER'
    assert data.program.version == '2.0.0'

    # backup partial system parsing
    assert len(data.model_system) == 1
    system = data.model_system[0]
    assert len(system.particle_states) == 12
    assert [s.atomic_number for s in system.particle_states] == [
        72,
        72,
        72,
        72,
        23,
        23,
        23,
        23,
        23,
        23,
        23,
        23,
    ]
    assert system.periodic_boundary_conditions == [True, True, True]

    # method
    assert len(data.model_method) == 1
    method = data.model_method[0]
    assert len(method.numerical_settings[0].basis_set_components) == 2
    assert method.numerical_settings[0].basis_set_components[0].basis_set == 'Koga'

    assert len(data.outputs) == 1
    output = data.outputs[0]
    assert len(output.x_lobster_abs_total_spilling) == 1
    assert output.x_lobster_abs_total_spilling[0] == approx(4.39)
    assert len(output.x_lobster_abs_charge_spilling) == 1
    assert output.x_lobster_abs_charge_spilling[0] == approx(2.21)

    # ICOHPLIST.lobster
    cohp = output.x_lobster_section_cohp
    assert cohp.x_lobster_number_of_cohp_pairs == 56
    assert len(cohp.x_lobster_cohp_atom1_labels) == 56
    assert cohp.x_lobster_cohp_atom1_labels[41] == 'V6'
    assert len(cohp.x_lobster_cohp_atom2_labels) == 56
    assert cohp.x_lobster_cohp_atom2_labels[16] == 'V9'
    assert len(cohp.x_lobster_cohp_distances) == 56
    assert cohp.x_lobster_cohp_distances[0].to('angstrom').magnitude == approx(3.17294)
    assert cohp.x_lobster_cohp_distances[47].to('angstrom').magnitude == approx(2.60684)
    assert cohp.x_lobster_cohp_distances[55].to('angstrom').magnitude == approx(2.55809)
    assert cohp.x_lobster_cohp_translations is None
    assert len(cohp.x_lobster_cohp_number_of_bonds) == 56
    assert cohp.x_lobster_cohp_number_of_bonds[0] == 2
    assert cohp.x_lobster_cohp_number_of_bonds[53] == 1
    assert np.shape(cohp.x_lobster_integrated_cohp_at_fermi_level) == (1, 56)
    assert cohp.x_lobster_integrated_cohp_at_fermi_level[0, 0].to(
        'eV'
    ).magnitude == approx(-1.72125)
    assert cohp.x_lobster_integrated_cohp_at_fermi_level[0, 55].to(
        'eV'
    ).magnitude == approx(-1.62412)

    # ICOOPLIST.lobster
    coop = output.x_lobster_section_coop
    assert coop.x_lobster_number_of_coop_pairs == 56
    assert len(coop.x_lobster_coop_atom1_labels) == 56
    assert coop.x_lobster_coop_atom1_labels[41] == 'V6'
    assert len(coop.x_lobster_coop_atom2_labels) == 56
    assert coop.x_lobster_coop_atom2_labels[11] == 'Hf4'
    assert len(coop.x_lobster_coop_distances) == 56
    assert coop.x_lobster_coop_distances[0].to('angstrom').magnitude == approx(3.17294)
    assert coop.x_lobster_coop_distances[47].to('angstrom').magnitude == approx(2.60684)
    assert coop.x_lobster_coop_distances[55].to('angstrom').magnitude == approx(2.55809)
    assert coop.x_lobster_coop_translations is None
    assert len(coop.x_lobster_coop_number_of_bonds) == 56
    assert coop.x_lobster_coop_number_of_bonds[0] == 2
    assert coop.x_lobster_coop_number_of_bonds[53] == 1
    assert np.shape(coop.x_lobster_integrated_coop_at_fermi_level) == (1, 56)
    assert coop.x_lobster_integrated_coop_at_fermi_level[0, 0].magnitude == approx(
        -0.46493
    )
    assert coop.x_lobster_integrated_coop_at_fermi_level[0, 55].magnitude == approx(
        -0.50035
    )


@pytest.mark.skipif(
    Version(ase.__version__) > Version('3.22'), reason='Incompatible with ase v26'
)
def test_QE_Ni(parser):
    """
    Check that basic info is parsed properly when LOBSTER is run on top
    of Quantum Espresso calculations.
    """

    archive = EntryArchive()
    parser.parse('tests/data/lobster/Ni/lobsterout', archive, logger)

    data = archive.data

    # QE system parsing
    assert len(data.model_system) == 1
    system = data.model_system[0]
    assert len(system.particle_states) == 1
    assert [s.chemical_symbol for s in system.particle_states] == ['Ni']
    assert system.periodic_boundary_conditions == [True, True, True]
    assert len(system.positions) == 1
    assert (system.positions[0].magnitude == [0, 0, 0]).all()

    assert len(data.model_method) == 1
    method = data.model_method[0]
    assert method.x_lobster_code == 'Quantum Espresso'
    assert len(method.numerical_settings[0].basis_set_components) == 1
    assert method.numerical_settings[0].basis_set_components[0].basis_set == 'Bunge'
    # assert method.x_lobster_basis_functions == {
    #     'Ni': [
    #         '4s',
    #         '3p_y',
    #         '3p_z',
    #         '3p_x',
    #         '3d_xy',
    #         '3d_yz',
    #         '3d_z^2',
    #         '3d_xz',
    #         '3d_x^2-y^2',
    #     ]
    # }

    assert len(data.outputs) == 1
    output = data.outputs[0]
    assert len(output.x_lobster_abs_total_spilling) == 2
    assert output.x_lobster_abs_total_spilling[0] == approx(36.14)
    assert output.x_lobster_abs_total_spilling[1] == approx(36.11)
    assert len(output.x_lobster_abs_charge_spilling) == 2
    assert output.x_lobster_abs_charge_spilling[0] == approx(4.02)
    assert output.x_lobster_abs_charge_spilling[1] == approx(3.37)


def test_Si(parser):  # noqa: PLR0915
    """
    Test spin-polarized orbitalwise Si calculation with LOBSTER 4.1.0,
    it has different ICOHPLIST.lobster and ICOOPLIST.lobster scheme.
    """

    archive = EntryArchive()
    parser.parse('tests/data/lobster/Si/lobsterout.gz', archive, logger)

    data = archive.data
    assert data.program.name == 'LOBSTER'
    assert data.program.version == '4.1.0'

    assert len(data.outputs) == 1
    output = data.outputs[0]
    assert len(output.x_lobster_abs_total_spilling) == 2
    assert output.x_lobster_abs_total_spilling[0] == approx(17.91)
    assert len(output.x_lobster_abs_charge_spilling) == 2
    assert output.x_lobster_abs_charge_spilling[0] == approx(1.42)

    # backup partial system parsing
    assert len(data.model_system) == 1
    system = data.model_system[0]
    assert len(system.particle_states) == 2
    assert [s.chemical_symbol for s in system.particle_states] == ['Si', 'Si']
    assert system.periodic_boundary_conditions == [True, True, True]

    # method
    assert len(data.model_method) == 1
    method = data.model_method[0]
    assert len(method.numerical_settings[0].basis_set_components) == 1
    assert (
        method.numerical_settings[0].basis_set_components[0].basis_set
        == 'pbeVaspFit2015'
    )
    # assert method.x_lobster_basis_functions == {
    #     'Si': [
    #         '3s',
    #         '3p_y',
    #         '3p_z',
    #         '3p_x',
    #     ]
    # }

    # ICOHPLIST.lobster
    cohp = output.x_lobster_section_cohp
    assert cohp.x_lobster_number_of_cohp_pairs == 64
    assert len(cohp.x_lobster_cohp_atom1_labels) == 64
    assert len(cohp.x_lobster_cohp_atom2_labels) == 64
    assert len(cohp.x_lobster_cohp_distances) == 64
    assert cohp.x_lobster_cohp_distances[0].to('angstrom').magnitude == approx(5.468728)
    assert cohp.x_lobster_cohp_distances[47].to('angstrom').magnitude == approx(3.86697)
    assert cohp.x_lobster_cohp_distances[23].to('angstrom').magnitude == approx(4.53443)
    assert np.array_equal(cohp.x_lobster_cohp_translations[26], [-1, 1, -1])
    assert np.shape(cohp.x_lobster_integrated_cohp_at_fermi_level) == (2, 64)
    assert cohp.x_lobster_integrated_cohp_at_fermi_level[0, 0].to(
        'eV'
    ).magnitude == approx(-0.00058)
    assert cohp.x_lobster_integrated_cohp_at_fermi_level[0, 31].to(
        'eV'
    ).magnitude == approx(-2.24755)

    # ICOOPLIST.lobster
    coop = output.x_lobster_section_coop
    assert coop.x_lobster_number_of_coop_pairs == 64
    assert len(coop.x_lobster_coop_atom1_labels) == 64
    assert len(coop.x_lobster_coop_atom2_labels) == 64
    assert len(coop.x_lobster_coop_distances) == 64

    # check if ICOBILIST.lobster is correctly read
    cobi = output.x_lobster_section_cobi
    assert cobi.x_lobster_number_of_cobi_pairs == 64
    assert len(cobi.x_lobster_cobi_atom1_labels) == 64
    assert len(cobi.x_lobster_cobi_atom2_labels) == 64
    assert len(cobi.x_lobster_cobi_distances) == 64

    # check if orbital-wise data is correctly read
    assert len(coop.x_lobster_coop_orbital_per_label) == 64
    assert len(cohp.x_lobster_cohp_orbital_per_label) == 64
    assert len(cobi.x_lobster_cobi_orbital_per_label) == 64
    assert (
        coop.x_lobster_coop_orbital_per_label[-1]
        .x_lobster_orbital_pairs[0]
        .x_lobster_atom1_orbital
        == 'Si2_3s'
    )
    assert (
        coop.x_lobster_coop_orbital_per_label[-1]
        .x_lobster_orbital_pairs[0]
        .x_lobster_atom2_orbital
        == 'Si2_3s'
    )
    assert (
        coop.x_lobster_coop_orbital_per_label[10]
        .x_lobster_orbital_pairs[1]
        .x_lobster_atom1_orbital
        == 'Si1_3p_y'
    )
    assert (
        coop.x_lobster_coop_orbital_per_label[10]
        .x_lobster_orbital_pairs[1]
        .x_lobster_atom2_orbital
        == 'Si1_3s'
    )
    assert cohp.x_lobster_cohp_orbital_per_label[24].x_lobster_orbital_pairs[
        1
    ].x_lobster_integrated_cohp_orbital_values[0][5] == approx(-0.2004)
    assert cobi.x_lobster_cobi_orbital_per_label[20].x_lobster_orbital_pairs[
        1
    ].x_lobster_integrated_cobi_orbital_values[0][5] == approx(0.00052)
    assert (
        coop.x_lobster_integrated_coop_values[24, 1, 5].magnitude
        == coop.x_lobster_integrated_coop_at_fermi_level[0][24].magnitude
    )

    # test if data is parsed correctly by matching data from icoxplist with coxpcar
    for spin in [0, 1]:
        for ix, icohp in enumerate(cohp.x_lobster_integrated_cohp_at_fermi_level[spin]):
            assert np.isclose(
                icohp.magnitude,
                cohp.x_lobster_integrated_cohp_values[ix, spin, 5].magnitude,
            )
