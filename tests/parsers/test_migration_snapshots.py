"""Snapshot-based verification for PR #170 migration.

Run at baseline (1c7acd6) and current (test-data-normalization) to verify
migration only adds properties, doesn't modify existing ones.
"""

import pytest
from nomad.datamodel import EntryArchive
from nomad_simulation_parsers.parsers.vasp.parser import VASPParser
from nomad_simulation_parsers.parsers.abinit.parser import AbinitParser
from nomad_simulation_parsers.parsers.gpaw.parser import GPAWParser
from nomad_simulation_parsers.parsers.octopus.parser import OctopusParser
from nomad_simulation_parsers.parsers.ams.parser import AMSParser
from nomad_simulation_parsers.parsers.exciting.parser import ExcitingParser
from nomad_simulation_parsers.parsers.crystal.parser import CrystalParser
from nomad_simulation_parsers.parsers.fhiaims.parser import FHIAimsParser
from nomad_simulation_parsers.parsers.gromacs.parser import GromacsParser
from nomad_simulation_parsers.parsers.h5md.parser import H5MDParser
from nomad_simulation_parsers.parsers.lammps.parser import LammpsParser
from nomad_simulation_parsers.parsers.phonopy.parser import PhonopyParser
from nomad_simulation_parsers.parsers.quantumespresso.parser import QuantumEspressoParser
from nomad_simulation_parsers.parsers.wannier90.parser import Wannier90Parser
from nomad_simulation_parsers.parsers.yambo.parser import YamboParser


def serialize_archive_data(archive):
    """Convert archive.data to dict for snapshot comparison."""
    if not archive.data:
        return None
    return archive.data.m_to_dict(with_meta=False)


# ============= VASP =============


def test_vasp_vasprun_snapshot(snapshot):
    """Snapshot of VASP vasprun.xml parsing output."""
    parser = VASPParser()
    archive = EntryArchive()
    parser.parse('tests/data/vasp/AgAc_relax/vasprun.xml.relax', archive, None)
    assert serialize_archive_data(archive) == snapshot


def test_vasp_outcar_snapshot(snapshot):
    """Snapshot of VASP OUTCAR parsing output."""
    parser = VASPParser()
    archive = EntryArchive()
    parser.parse('tests/data/vasp/AgAc_relax/OUTCAR', archive, None)
    assert serialize_archive_data(archive) == snapshot


# ============= ABINIT =============


def test_abinit_snapshot(snapshot):
    """Snapshot of ABINIT parsing output."""
    parser = AbinitParser()
    archive = EntryArchive()
    parser.parse('tests/data/abinit/Fe/abinit.out', archive, None)
    assert serialize_archive_data(archive) == snapshot


# ============= GPAW =============


def test_gpaw_snapshot(snapshot):
    """Snapshot of GPAW parsing output."""
    parser = GPAWParser()
    archive = EntryArchive()
    parser.parse('tests/data/gpaw/Fe2.gpw', archive, None)
    assert serialize_archive_data(archive) == snapshot


# ============= Octopus =============


def test_octopus_snapshot(snapshot):
    """Snapshot of Octopus parsing output."""
    parser = OctopusParser()
    archive = EntryArchive()
    parser.parse('tests/data/octopus/Fe_spinpol/inp', archive, None)
    assert serialize_archive_data(archive) == snapshot


# ============= AMS =============


def test_ams_snapshot(snapshot):
    """Snapshot of AMS parsing output."""
    parser = AMSParser()
    archive = EntryArchive()
    parser.parse('tests/data/ams/scf/ams.log', archive, None)
    assert serialize_archive_data(archive) == snapshot


# ============= Exciting =============


def test_exciting_snapshot(snapshot):
    """Snapshot of Exciting parsing output."""
    parser = ExcitingParser()
    archive = EntryArchive()
    parser.parse('tests/data/exciting/C_minimal/INFO.OUT', archive, None)
    assert serialize_archive_data(archive) == snapshot


# ============= Crystal =============


def test_crystal_snapshot(snapshot):
    """Snapshot of Crystal parsing output."""
    parser = CrystalParser()
    archive = EntryArchive()
    parser.parse('tests/data/crystal/single_point/dft/output.out', archive, None)
    assert serialize_archive_data(archive) == snapshot


# ============= FHI-aims =============


def test_fhiaims_snapshot(snapshot):
    """Snapshot of FHI-aims parsing output."""
    parser = FHIAimsParser()
    archive = EntryArchive()
    parser.parse('tests/data/fhiaims/Si_geomopt/out.out', archive, None)
    assert serialize_archive_data(archive) == snapshot


# ============= GROMACS =============


def test_gromacs_snapshot(snapshot):
    """Snapshot of GROMACS parsing output."""
    parser = GromacsParser()
    archive = EntryArchive()
    parser.parse('tests/data/gromacs/fe_test/md.log', archive, None)
    assert serialize_archive_data(archive) == snapshot


# ============= H5MD =============


def test_h5md_snapshot(snapshot):
    """Snapshot of H5MD parsing output."""
    parser = H5MDParser()
    archive = EntryArchive()
    parser.parse('tests/data/h5md/test_traj_openmm_reduced-SOL_5frames_07-10-25.h5', archive, None)
    assert serialize_archive_data(archive) == snapshot


# ============= LAMMPS =============


def test_lammps_snapshot(snapshot):
    """Snapshot of LAMMPS parsing output."""
    parser = LammpsParser()
    archive = EntryArchive()
    parser.parse('tests/data/lammps/1_xyz_files/log.lammps', archive, None)
    assert serialize_archive_data(archive) == snapshot


# ============= Phonopy =============


def test_phonopy_snapshot(snapshot):
    """Snapshot of Phonopy parsing output."""
    parser = PhonopyParser()
    archive = EntryArchive()
    parser.parse('tests/data/phonopy/vasp/phonopy.yaml', archive, None)
    assert serialize_archive_data(archive) == snapshot


# ============= Quantum Espresso =============


def test_quantumespresso_snapshot(snapshot):
    """Snapshot of Quantum Espresso parsing output."""
    parser = QuantumEspressoParser()
    archive = EntryArchive()
    parser.parse('tests/data/quantumespresso/pwscf/TiO2_opt/pw.out', archive, None)
    assert serialize_archive_data(archive) == snapshot


# ============= Wannier90 =============


def test_wannier90_snapshot(snapshot):
    """Snapshot of Wannier90 parsing output."""
    parser = Wannier90Parser()
    archive = EntryArchive()
    parser.parse('tests/data/wannier90/lco_mlwf/lco.wout', archive, None)
    assert serialize_archive_data(archive) == snapshot
