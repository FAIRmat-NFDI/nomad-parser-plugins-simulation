import numpy as np
import pytest
from nomad.datamodel import EntryArchive
from nomad.utils import get_logger

from nomad_simulation_parsers.parsers.lammps.parser import LammpsParser


@pytest.mark.integration
def test_stores_input_structure_without_a_trajectory(tmp_path):
    (tmp_path / 'input.data').write_text(
        'LAMMPS data file\n\n'
        '2 atoms\n1 atom types\n\n'
        '0 10 xlo xhi\n0 12 ylo yhi\n0 14 zlo zhi\n\n'
        'Masses\n\n1 12.01\n\n'
        'Atoms\n\n1 1 0.5 1.0 1.5\n2 1 2.0 3.0 4.0\n'
    )
    mainfile = tmp_path / 'log.lammps'
    mainfile.write_text(
        'LAMMPS (test)\nunits real\natom_style atomic\nboundary p p f\n'
        'read_data input.data\nminimize 1e-4 1e-6 100 1000\n'
    )

    archive = EntryArchive()
    LammpsParser().parse(str(mainfile), archive, get_logger(__name__))

    assert archive.workflow2 is None
    assert len(archive.data.model_system) == 1
    system = archive.data.model_system[0]
    assert system.name == 'Input structure'
    assert not system.is_representative
    assert system.n_particles == 2
    assert system.periodic_boundary_conditions == [True, True, False]
    assert np.diag(system.lattice_vectors.to('angstrom').magnitude) == pytest.approx(
        [10, 12, 14]
    )
    assert system.positions[0].to('angstrom').magnitude == pytest.approx(
        [0.5, 1.0, 1.5]
    )
    assert len(system.particle_states) == 2


@pytest.mark.integration
def test_missing_input_data_does_not_invent_a_system(tmp_path):
    mainfile = tmp_path / 'log.lammps'
    mainfile.write_text(
        'LAMMPS (test)\nunits real\natom_style atomic\n'
        'read_data missing.data\nminimize 1e-4 1e-6 100 1000\n'
    )

    archive = EntryArchive()
    LammpsParser().parse(str(mainfile), archive, get_logger(__name__))

    assert archive.workflow2 is None
    assert archive.data.model_system == []


@pytest.mark.integration
def test_uses_redirected_log_data_when_dump_is_missing(tmp_path):
    (tmp_path / 'input.data').write_text(
        'LAMMPS data file\n\n1 atoms\n1 atom types\n\n'
        '0 10 xlo xhi\n0 10 ylo yhi\n0 10 zlo zhi\n\n'
        'Masses\n\n1 12.01\n\nAtoms\n\n1 1 1.0 2.0 3.0\n'
    )
    mainfile = tmp_path / 'log.lammps'
    mainfile.write_text(
        'LAMMPS (test)\nunits metal\natom_style atomic\nlog redirected.log\n'
    )
    (tmp_path / 'redirected.log').write_text(
        'LAMMPS (test)\nread_data input.data\n'
        'dump d all xyz 10 missing.xyz\nminimize 1e-4 1e-6 100 1000\n'
    )

    archive = EntryArchive()
    LammpsParser().parse(str(mainfile), archive, get_logger(__name__))

    assert archive.workflow2 is None
    assert len(archive.data.model_system) == 1
    system = archive.data.model_system[0]
    assert system.n_particles == 1
    assert system.positions[0].to('angstrom').magnitude == pytest.approx(
        [1.0, 2.0, 3.0]
    )
    assert not system.is_representative


@pytest.mark.integration
def test_missing_data_and_dump_leave_no_system(tmp_path):
    mainfile = tmp_path / 'log.lammps'
    mainfile.write_text(
        'LAMMPS (test)\nunits metal\natom_style atomic\n'
        'read_data missing.data\n'
        'dump d all custom 10 missing.lammpstrj id type x y z\n'
        'fix md all nve\nrun 100\n'
    )

    archive = EntryArchive()
    LammpsParser().parse(str(mainfile), archive, get_logger(__name__))

    assert archive.workflow2.method.n_steps == 100
    assert archive.data.model_system == []


@pytest.mark.integration
def test_input_structure_uses_declared_distance_units(tmp_path):
    (tmp_path / 'input.data').write_text(
        'LAMMPS data file\n\n1 atoms\n1 atom types\n\n'
        '0 10 xlo xhi\n0 10 ylo yhi\n0 10 zlo zhi\n\n'
        'Masses\n\n1 1.0\n\nAtoms\n\n1 1 0.5 1.0 1.5\n'
    )
    mainfile = tmp_path / 'log.lammps'
    mainfile.write_text(
        'LAMMPS (test)\nunits si\natom_style atomic\n'
        'read_data input.data\nminimize 1e-4 1e-6 100 1000\n'
    )

    archive = EntryArchive()
    LammpsParser().parse(str(mainfile), archive, get_logger(__name__))

    system = archive.data.model_system[0]
    assert system.positions[0].to('meter').magnitude == pytest.approx([0.5, 1.0, 1.5])
    assert np.diag(system.lattice_vectors.to('meter').magnitude) == pytest.approx(
        [10, 10, 10]
    )


@pytest.mark.integration
def test_nearby_data_without_read_data_is_not_assumed_to_be_the_input(tmp_path):
    (tmp_path / 'other.data').write_text(
        'LAMMPS data file\n\n1 atoms\n1 atom types\n\n'
        '0 10 xlo xhi\n0 10 ylo yhi\n0 10 zlo zhi\n\n'
        'Masses\n\n1 1.0\n\nAtoms\n\n1 1 0.5 1.0 1.5\n'
    )
    mainfile = tmp_path / 'log.lammps'
    mainfile.write_text('LAMMPS (test)\nunits real\natom_style atomic\n')

    archive = EntryArchive()
    LammpsParser().parse(str(mainfile), archive, get_logger(__name__))

    assert archive.data.model_system == []
