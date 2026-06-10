"""Snapshot-based verification for PR #170 migration.

Run at baseline (1c7acd6) and current (test-data-normalization) to verify
migration only adds properties, doesn't modify existing ones.
"""

import pytest
from nomad.datamodel import EntryArchive
from nomad_simulation_parsers.parsers import (
    ams_parser,
    abinit_parser,
    exciting_parser,
    gpaw_parser,
    octopus_parser,
    vasp_parser,
)


def serialize_archive_data(archive):
    """Convert archive.data to dict for snapshot comparison."""
    if not archive.data:
        return None
    return archive.data.m_to_dict(with_meta=False)


# ============= VASP =============


def test_vasp_vasprun_snapshot(snapshot):
    """Snapshot of VASP vasprun.xml parsing output."""
    parser = vasp_parser.VASPParser()
    archive = EntryArchive()
    parser.parse('tests/data/vasp/AgAc_relax/vasprun.xml.relax', archive, None)
    assert serialize_archive_data(archive) == snapshot


def test_vasp_outcar_snapshot(snapshot):
    """Snapshot of VASP OUTCAR parsing output."""
    parser = vasp_parser.VASPParser()
    archive = EntryArchive()
    parser.parse('tests/data/vasp/AgAc_relax/OUTCAR', archive, None)
    assert serialize_archive_data(archive) == snapshot


# ============= ABINIT =============


def test_abinit_snapshot(snapshot):
    """Snapshot of ABINIT parsing output."""
    parser = abinit_parser.ABINITParser()
    archive = EntryArchive()
    parser.parse('tests/data/abinit/Fe/abinit.out', archive, None)
    assert serialize_archive_data(archive) == snapshot


# ============= GPAW =============


def test_gpaw_snapshot(snapshot):
    """Snapshot of GPAW parsing output."""
    parser = gpaw_parser.GPAWParser()
    archive = EntryArchive()
    parser.parse('tests/data/gpaw/Fe2.gpw', archive, None)
    assert serialize_archive_data(archive) == snapshot


# ============= Octopus =============


def test_octopus_snapshot(snapshot):
    """Snapshot of Octopus parsing output."""
    parser = octopus_parser.OctopusParser()
    archive = EntryArchive()
    parser.parse('tests/data/octopus/Fe_spinpol/inp', archive, None)
    assert serialize_archive_data(archive) == snapshot


# ============= AMS =============


def test_ams_snapshot(snapshot):
    """Snapshot of AMS parsing output."""
    parser = ams_parser.AMSParser()
    archive = EntryArchive()
    parser.parse('tests/data/ams/scf/ams.log', archive, None)
    assert serialize_archive_data(archive) == snapshot


# ============= Exciting =============


def test_exciting_snapshot(snapshot):
    """Snapshot of Exciting parsing output."""
    parser = exciting_parser.ExcitingParser()
    archive = EntryArchive()
    parser.parse('tests/data/exciting/C_minimal/INFO.OUT', archive, None)
    assert serialize_archive_data(archive) == snapshot
