"""Frame-order regressions for LAMMPS model systems."""

import pytest
from nomad.client import normalize_all
from nomad.datamodel import EntryArchive, EntryMetadata
from nomad.utils import get_logger

from nomad_simulation_parsers.parsers.lammps.parser import LammpsParser


def _dump_frame(step: int, x_position: float, velocity: float | None = None) -> str:
    columns = 'id type x y z'
    atom = f'1 1 {x_position} 2 3'
    if velocity is not None:
        columns += ' vx vy vz'
        atom += f' {velocity} 0 0'
    return (
        'ITEM: TIMESTEP\n'
        f'{step}\n'
        'ITEM: NUMBER OF ATOMS\n'
        '1\n'
        'ITEM: BOX BOUNDS pp pp pp\n'
        '0 10\n0 10\n0 10\n'
        f'ITEM: ATOMS {columns}\n'
        f'{atom}\n'
    )


@pytest.mark.integration
@pytest.mark.parametrize('steps', [(10, 10), (20, 10)])
def test_model_systems_keep_file_order_when_steps_repeat_or_decrease(tmp_path, steps):
    mainfile = tmp_path / 'log.lammps'
    mainfile.write_text(
        'LAMMPS (test)\nunits real\n'
        'dump 1 all custom 1 trajectory.lammpstrj id type x y z\n'
        'run 20\n'
    )
    (tmp_path / 'trajectory.lammpstrj').write_text(
        _dump_frame(steps[0], 1.0) + _dump_frame(steps[1], 4.0)
    )
    archive = EntryArchive()

    LammpsParser().parse(str(mainfile), archive, get_logger(__name__))

    systems = archive.data.model_system
    assert len(systems) == 2
    assert [system.positions[0][0].to('angstrom').magnitude for system in systems] == [
        pytest.approx(1.0),
        pytest.approx(4.0),
    ]
    assert not systems[0].is_representative
    assert systems[1].is_representative

    normalized = EntryArchive.m_from_dict(archive.m_to_dict())
    normalized.metadata = EntryMetadata()
    normalize_all(normalized, logger=get_logger(__name__))
    normalized_systems = normalized.data.model_system
    assert len(normalized_systems) == 2
    assert normalized_systems[1].positions[0][0].to('angstrom').magnitude == (
        pytest.approx(4.0)
    )


@pytest.mark.integration
@pytest.mark.parametrize(
    ('unit_style', 'velocity_m_per_s'), [('real', 20000), ('metal', 20)]
)
def test_model_systems_store_source_velocities_with_physical_units(
    tmp_path, unit_style, velocity_m_per_s
):
    mainfile = tmp_path / 'log.lammps'
    mainfile.write_text(
        f'LAMMPS (test)\nunits {unit_style}\n'
        'dump 1 all custom 1 trajectory.lammpstrj id type x y z vx vy vz\n'
        'run 1\n'
    )
    (tmp_path / 'trajectory.lammpstrj').write_text(_dump_frame(1, 1.0, 0.2))
    archive = EntryArchive()

    LammpsParser().parse(str(mainfile), archive, get_logger(__name__))

    system = archive.data.model_system[0]
    assert system.velocities.shape == (1, 3)
    assert system.velocities[0][0].to('meter/second').magnitude == pytest.approx(
        velocity_m_per_s
    )
    serialized = EntryArchive.m_from_dict(archive.m_to_dict())
    assert serialized.data.model_system[0].velocities[0][0].to(
        'meter/second'
    ).magnitude == pytest.approx(velocity_m_per_s)


@pytest.mark.integration
def test_fallback_bonds_use_trajectory_row_indices_not_lammps_atom_ids(tmp_path):
    mainfile = tmp_path / 'log.lammps'
    mainfile.write_text('LAMMPS (test)\nunits real\nread_data data.test\nrun 1\n')
    (tmp_path / 'data.test').write_text(
        'LAMMPS data file\n\n'
        '2 atoms\n1 bonds\n1 atom types\n1 bond types\n\n'
        'Masses\n\n1 12.011\n\n'
        'Atoms # full\n\n'
        '10 1 1 0.0 1 2 3\n'
        '20 1 1 0.0 4 5 6\n\n'
        'Bonds\n\n1 1 10 20\n'
    )
    (tmp_path / 'trajectory.lammpstrj').write_text(
        'ITEM: TIMESTEP\n1\n'
        'ITEM: NUMBER OF ATOMS\n2\n'
        'ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n'
        'ITEM: ATOMS id type x y z\n'
        '20 1 4 5 6\n'
        '10 1 1 2 3\n'
    )
    archive = EntryArchive()

    LammpsParser().parse(str(mainfile), archive, get_logger(__name__))

    system = archive.data.model_system[0]
    assert system.n_particles == 2
    assert system.positions[:, 0].to('angstrom').magnitude.tolist() == [1, 4]
    assert system.bond_list.tolist() == [[0, 1]]
    serialized = EntryArchive.m_from_dict(archive.m_to_dict())
    assert serialized.data.model_system[0].bond_list.tolist() == [[0, 1]]


@pytest.mark.integration
def test_fallback_bonds_are_omitted_without_trajectory_atom_ids(tmp_path):
    mainfile = tmp_path / 'log.lammps'
    mainfile.write_text('LAMMPS (test)\nunits real\nread_data data.test\nrun 1\n')
    (tmp_path / 'data.test').write_text(
        'LAMMPS data file\n\n'
        '2 atoms\n1 bonds\n1 atom types\n1 bond types\n\n'
        'Masses\n\n1 12.011\n\n'
        'Atoms # full\n\n'
        '1 1 1 0.0 1 2 3\n'
        '2 1 1 0.0 4 5 6\n\n'
        'Bonds\n\n1 1 1 2\n'
    )
    (tmp_path / 'trajectory.lammpstrj').write_text(
        'ITEM: TIMESTEP\n1\n'
        'ITEM: NUMBER OF ATOMS\n2\n'
        'ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n'
        'ITEM: ATOMS type x y z\n'
        '1 1 2 3\n'
        '1 4 5 6\n'
    )
    archive = EntryArchive()

    LammpsParser().parse(str(mainfile), archive, get_logger(__name__))

    system = archive.data.model_system[0]
    assert system.n_particles == 2
    assert system.bond_list is None
