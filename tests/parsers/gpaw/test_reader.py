import numpy as np
import pytest

from nomad_simulation_parsers.parsers.gpaw.gpw_parser import (
    GPW2FileParser,
    GPWFileParser,
    GPWTarParser,
)
from tests.parsers.common import approx, assert_approx

from .conftest import DATA_DIR


@pytest.mark.unit
class TestGPWTarReader:
    @pytest.fixture
    def parser(self):
        parser = GPWTarParser()
        parser.mainfile = str(DATA_DIR / 'Fe2.gpw')
        return parser

    def test_reads_metadata(self, parser):
        assert parser.info['parameter']['energyunit'] == 'Hartree'
        assert parser.get_program_version() == '1.1.0'
        assert parser.get_parameter('version') == 6
        assert parser.get_parameter('energy_total') == approx(-0.2683163515925333)
        assert parser.get_smearing_width() == approx(0.003674932247495664)

    def test_reads_array_data_and_dimensions(self, parser):
        np.testing.assert_equal(parser.get_array('atomicnumbers'), [26, 26])
        np.testing.assert_equal(parser.get_array('atom_positions').shape, (2, 3))
        np.testing.assert_equal(parser.get_array('unitcell').shape, (3, 3))
        np.testing.assert_equal(parser.get_array('boundaryconditions'), [1, 1, 1])
        np.testing.assert_equal(parser.get_array('eigenvalues').shape, (2, 10, 18))
        np.testing.assert_equal(parser.get_array('occupation').shape, (2, 10, 18))
        assert parser.get_array_dimension('ngpts') == [12, 12, 12]


@pytest.mark.unit
class TestGPW2Reader:
    @pytest.fixture
    def parser(self):
        parser = GPW2FileParser()
        parser.mainfile = str(DATA_DIR / 'Si_pw.gpw2')
        return parser

    def test_reads_metadata(self, parser):
        assert parser.info['parameter']['mode']['name'] == 'pw'
        assert parser.get_program_version()
        assert parser.get_parameter('xc') == 'LDA'
        assert parser.get_parameter('energyerror') == approx(0.0005)
        assert parser.get_smearing_width() == approx(0.1)

    def test_reads_array_data_and_dimensions(self, parser):
        np.testing.assert_equal(parser.get_array('atomicnumbers').shape, (2,))
        np.testing.assert_equal(parser.get_array('atom_positions').shape, (2, 3))
        np.testing.assert_equal(
            np.asarray(parser.get_array('boundary_conditions')).shape, (3,)
        )
        np.testing.assert_equal(parser.get_array('eigenvalues').shape, (1, 10, 8))
        np.testing.assert_equal(parser.get_array('occupation').shape, (1, 10, 8))
        assert parser.get_array_dimension('ngpts') == (1, 14, 14, 14)

    def test_reads_lcao_density_and_effective_potential(self):
        parser = GPW2FileParser()
        parser.mainfile = str(DATA_DIR / 'Si_lcao.gpw2')

        density = parser.get_array('density')
        potential = parser.get_array('potential_effective')
        assert density.shape == (1, 16, 16, 16)
        assert potential.shape == (1, 16, 16, 16)
        assert density[0, 13, 2, 8] == approx(0.08703530984156761)
        assert potential[0, 5, 14, 1] == approx(-7.709631839033944)


@pytest.mark.unit
class TestGPWFileReader:
    def test_parses_gpw_and_applies_source_units(self):
        parser = GPWFileParser()
        parser.mainfile = str(DATA_DIR / 'Fe2.gpw')
        parser.parse()

        assert parser.get_mode() == 'pw'
        assert parser.results['labels'] == ['Fe', 'Fe']
        assert parser.results['unitcell'].to('angstrom').magnitude.shape == (3, 3)
        assert_approx(
            parser.apply_unit(2.0, 'energyunit').to('hartree').magnitude,
            2.0,
        )
