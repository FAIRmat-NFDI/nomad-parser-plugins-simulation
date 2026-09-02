#!/usr/bin/env python3
"""Generate a JSON serialization of each parser's ``archive.data`` for the
additivity check.

This is not a pytest test: it is run directly by ``check_additivity.sh`` (and
manually) to snapshot a ref. Snapshots are written to a plain JSON file and
never committed -- they are regenerated from scratch for both the target and
source refs and compared by ``compare_snapshots.py``. A JSON *dump* (not a
hash) is used so the comparison can report *where* the archives differ.

Usage::

    generate_snapshots.py <out.json> [parser ...]

With no parser names, every parser in ``PARSERS`` is snapshotted; otherwise only
the named parser directories (e.g. ``orca vasp``) are.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

from nomad.datamodel import EntryArchive
from nomad.utils import get_logger

from nomad_simulation_parsers.parsers.abinit.parser import AbinitParser
from nomad_simulation_parsers.parsers.ams.parser import AMSParser
from nomad_simulation_parsers.parsers.crystal.parser import CrystalParser
from nomad_simulation_parsers.parsers.exciting.parser import ExcitingParser
from nomad_simulation_parsers.parsers.fhiaims.parser import FHIAimsParser
from nomad_simulation_parsers.parsers.gpaw.parser import GPAWParser
from nomad_simulation_parsers.parsers.gromacs.parser import GromacsParser
from nomad_simulation_parsers.parsers.h5md.parser import H5MDParser
from nomad_simulation_parsers.parsers.lammps.parser import LammpsParser
from nomad_simulation_parsers.parsers.octopus.parser import OctopusParser
from nomad_simulation_parsers.parsers.orca.parser import OrcaParser
from nomad_simulation_parsers.parsers.phonopy.parser import PhonopyParser
from nomad_simulation_parsers.parsers.quantumespresso.parser import (
    QuantumEspressoParser,
)
from nomad_simulation_parsers.parsers.vasp.parser import VASPParser
from nomad_simulation_parsers.parsers.wannier90.parser import Wannier90Parser

LOGGER = get_logger(__name__)

# Parser directory -> representative mainfiles. Keyed by the parser directory so
# selection matches the names `check_additivity.sh` derives from a diff.
PARSERS: dict[str, list[tuple[type, str]]] = {
    'vasp': [
        (VASPParser, 'tests/data/vasp/AgAc_relax/vasprun.xml.relax'),
        (VASPParser, 'tests/data/vasp/AgAc_relax/OUTCAR'),
    ],
    'abinit': [(AbinitParser, 'tests/data/abinit/Fe/abinit.out')],
    'gpaw': [(GPAWParser, 'tests/data/gpaw/Fe2.gpw')],
    'octopus': [(OctopusParser, 'tests/data/octopus/Fe_spinpol/inp')],
    'ams': [(AMSParser, 'tests/data/ams/scf/ams.log')],
    'exciting': [(ExcitingParser, 'tests/data/exciting/C_minimal/INFO.OUT')],
    'crystal': [(CrystalParser, 'tests/data/crystal/single_point/dft/output.out')],
    'fhiaims': [(FHIAimsParser, 'tests/data/fhiaims/Si_geomopt/out.out')],
    'gromacs': [(GromacsParser, 'tests/data/gromacs/fe_test/md.log')],
    'h5md': [
        (
            H5MDParser,
            'tests/data/h5md/test_traj_openmm_reduced-SOL_5frames_07-10-25.h5',
        )
    ],
    'lammps': [(LammpsParser, 'tests/data/lammps/1_xyz_files/log.lammps')],
    'phonopy': [(PhonopyParser, 'tests/data/phonopy/vasp/phonopy.yaml')],
    'quantumespresso': [
        (QuantumEspressoParser, 'tests/data/quantumespresso/pwscf/TiO2_opt/pw.out')
    ],
    'wannier90': [(Wannier90Parser, 'tests/data/wannier90/lco_mlwf/lco.wout')],
    'orca': [(OrcaParser, 'tests/data/orca/RI_MP2_water.out')],
}


def serialize_archive_data(archive: EntryArchive):
    """Convert archive.data to a JSON-serializable dict for comparison."""
    if not archive.data:
        return None
    return archive.data.m_to_dict(with_meta=False)


def generate(parser_dirs: list[str]) -> dict:
    snapshots = {}
    for parser_dir in parser_dirs:
        for parser_class, mainfile in PARSERS[parser_dir]:
            archive = EntryArchive()
            parser_class().parse(mainfile, archive, LOGGER)
            key = f'{parser_dir}:{Path(mainfile).name}'
            snapshots[key] = serialize_archive_data(archive)
    return snapshots


def main() -> int:
    if len(sys.argv) < 2:
        print(f'usage: {sys.argv[0]} <out.json> [parser ...]', file=sys.stderr)
        return 2
    out = sys.argv[1]
    requested = sys.argv[2:] or sorted(PARSERS)
    unknown = [name for name in requested if name not in PARSERS]
    if unknown:
        print(f'unknown parsers: {", ".join(unknown)}', file=sys.stderr)
        return 2
    snapshots = generate(requested)
    Path(out).write_text(json.dumps(snapshots, sort_keys=True, indent=1, default=str))
    print(f'wrote {len(snapshots)} snapshot(s) for {len(requested)} parser(s): {out}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
