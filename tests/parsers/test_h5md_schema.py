from pathlib import Path

import numpy as np
import pytest

h5py = pytest.importorskip('h5py')
pytest.importorskip('yaml')

from nomad_simulation_parsers.parsers.h5md.schema import (  # noqa: E402
    load_schema,
    validate_hdf5_file,
)


def _write_minimal_h5md_file(path: Path) -> None:
    with h5py.File(path, 'w') as h5_file:
        h5md = h5_file.create_group('h5md')
        h5md.attrs['version'] = np.array([1, 0], dtype=np.int32)

        author = h5md.create_group('author')
        author.attrs['name'] = 'Ada Lovelace'
        author.attrs['email'] = 'ada@example.org'

        creator = h5md.create_group('creator')
        creator.attrs['name'] = 'h5py'
        creator.attrs['version'] = '3.0.0'

        particles = h5_file.create_group('particles')
        particles_all = particles.create_group('all')
        position = particles_all.create_group('position')
        position.attrs['unit'] = 'angstrom'
        position.create_dataset('step', data=np.array([0, 1], dtype=np.int64))
        position.create_dataset('time', data=np.array([0.0, 0.5]))
        position.create_dataset('value', data=np.ones((2, 3, 3)))
        particles_all.create_dataset(
            'species_label',
            data=np.array([b'H', b'O', b'H']),
        )


def test_h5md_schema_validates_minimal_core_file(tmp_path):
    mainfile = tmp_path / 'minimal.h5'
    _write_minimal_h5md_file(mainfile)

    result = validate_hdf5_file(mainfile, load_schema())

    assert result.is_valid, [
        f'{issue.path}: {issue.message}' for issue in result.issues
    ]


def test_h5md_schema_reports_missing_required_header(tmp_path):
    mainfile = tmp_path / 'missing_creator.h5'
    _write_minimal_h5md_file(mainfile)
    with h5py.File(mainfile, 'a') as h5_file:
        del h5_file['h5md/creator']

    result = validate_hdf5_file(mainfile, load_schema())

    assert not result.is_valid
    assert any(issue.path == '/h5md/creator' for issue in result.issues)


def test_h5md_schema_validates_reference_fixture():
    mainfile = (
        Path(__file__).resolve().parents[1]
        / 'data'
        / 'h5md'
        / 'test_traj_openmm_reduced-SOL_5frames_07-10-25.h5'
    )

    result = validate_hdf5_file(mainfile, load_schema())

    assert result.is_valid, [
        f'{issue.path}: {issue.message}' for issue in result.issues
    ]
