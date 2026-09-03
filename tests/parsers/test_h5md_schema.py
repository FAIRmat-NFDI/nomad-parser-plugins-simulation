from pathlib import Path

import h5py
import numpy as np

from nomad_simulation_parsers.parsers.h5md.schema import (
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
        box = particles_all.create_group('box')
        box.attrs['dimension'] = np.int32(3)
        box.attrs['boundary'] = np.array([False, False, False])

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


def test_h5md_schema_checks_attribute_dimensions_against_datasets(tmp_path):
    mainfile = tmp_path / 'invalid_boundary_dimension.h5'
    _write_minimal_h5md_file(mainfile)
    with h5py.File(mainfile, 'a') as h5_file:
        box = h5_file['particles/all/box']
        box.attrs['boundary'] = np.array([True, True])

    result = validate_hdf5_file(mainfile, load_schema())

    assert not result.is_valid
    assert any(
        issue.path == '/particles/all/position/value'
        and "Expected dimension 'spatial_dimension'=2, got 3." == issue.message
        for issue in result.issues
    )


def test_h5md_schema_allows_optional_time_dataset(tmp_path):
    mainfile = tmp_path / 'without_time.h5'
    _write_minimal_h5md_file(mainfile)
    with h5py.File(mainfile, 'a') as h5_file:
        del h5_file['particles/all/position/time']

    result = validate_hdf5_file(mainfile, load_schema())

    assert result.is_valid, [
        f'{issue.path}: {issue.message}' for issue in result.issues
    ]


def test_h5md_schema_allows_different_element_frame_counts(tmp_path):
    mainfile = tmp_path / 'different_frame_counts.h5'
    _write_minimal_h5md_file(mainfile)
    with h5py.File(mainfile, 'a') as h5_file:
        velocity = h5_file['particles/all'].create_group('velocity')
        velocity.create_dataset('step', data=np.array([0], dtype=np.int64))
        velocity.create_dataset('value', data=np.ones((1, 3, 3)))

    result = validate_hdf5_file(mainfile, load_schema())

    assert result.is_valid, [
        f'{issue.path}: {issue.message}' for issue in result.issues
    ]


def test_h5md_schema_allows_string_boundary_and_static_edges(tmp_path):
    mainfile = tmp_path / 'string_boundary_static_edges.h5'
    _write_minimal_h5md_file(mainfile)
    with h5py.File(mainfile, 'a') as h5_file:
        box = h5_file['particles/all/box']
        box.attrs['boundary'] = np.array([b'periodic', b'periodic', b'none'])
        box.create_dataset('edges', data=np.eye(3))

    result = validate_hdf5_file(mainfile, load_schema())

    assert result.is_valid, [
        f'{issue.path}: {issue.message}' for issue in result.issues
    ]


def test_h5md_schema_requires_box_for_particles_group(tmp_path):
    mainfile = tmp_path / 'missing_box.h5'
    _write_minimal_h5md_file(mainfile)
    with h5py.File(mainfile, 'a') as h5_file:
        del h5_file['particles/all/box']

    result = validate_hdf5_file(mainfile, load_schema())

    assert not result.is_valid
    assert any(issue.path == '/particles/all/box' for issue in result.issues)


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
