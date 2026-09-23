import numpy as np
import pytest
from nomad.datamodel import EntryArchive
from nomad.utils import get_logger
from phonopy.structure.atoms import PhonopyAtoms

from nomad_simulation_parsers.parsers.phonopy.calculator import (
    generate_kpath_seekpath,
)
from nomad_simulation_parsers.parsers.phonopy.parser import PhonopyParser

LOGGER = get_logger(__name__)


def test_parse_file():
    parser = PhonopyParser()
    archive = EntryArchive()
    parser.parse('tests/data/phonopy/vasp/phonopy.yaml', archive, LOGGER)


@pytest.mark.parametrize(
    'cell, expected_n_segments',
    [
        pytest.param(np.diag([4.0] * 3), 6, id='cubic_canonical'),
        # Non-canonical orthorhombic cell (b < a < c) from a crashing upload:
        # ASE's canonical-cell matcher rejected this lattice-vector ordering,
        # while SeeKpath classifies it via spglib regardless of the ordering.
        pytest.param(
            np.diag([9.91421472, 4.39275307, 11.33757198]),
            12,
            id='orthorhombic_non_canonical',
        ),
    ],
)
def test_generate_kpath_seekpath(cell: np.ndarray, expected_n_segments: int):
    atoms = PhonopyAtoms(symbols=['Si'], cell=cell, scaled_positions=[[0.0, 0.0, 0.0]])

    parameters = generate_kpath_seekpath(atoms, 1e-5, LOGGER)

    assert len(parameters) == expected_n_segments
    for segment in parameters:
        assert {'npoints', 'startname', 'kstart', 'endname', 'kend'} <= segment.keys()
        assert segment['npoints'] == 100
    assert parameters[0]['startname'] == 'Γ'
