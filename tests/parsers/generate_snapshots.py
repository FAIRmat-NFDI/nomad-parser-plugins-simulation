#!/usr/bin/env python3
"""Generate a JSON serialization of each parser's ``archive.data`` for the
additivity check.

This is not a pytest test: it is run directly by ``check_additivity.sh`` (and
manually) to snapshot a ref. Snapshots are written to a plain JSON file and
never committed -- they are regenerated from scratch for both the target and
source refs and compared by ``compare_snapshots.py``. A JSON *dump* (not a
hash) is used so the comparison can report *where* the archives differ.

The output is ``{"provenance": {...}, "snapshots": {key: archive}}``. The
provenance block records the interpreter, serializer, harness, ``uv.lock`` and
fixture hashes (see ``provenance``) so a stale environment cannot masquerade as
an identical run; ``compare_snapshots.py`` reports it but never diffs it.

Usage::

    generate_snapshots.py <out.json> [parser ...]

With no parser names, every parser in ``PARSERS`` is snapshotted; otherwise only
the named parser directories (e.g. ``orca vasp``) are. A registered parser whose
fixture is not present in the checkout (yambo until one is added, or a large
fixture) is skipped rather than failing.
"""

from __future__ import annotations

import hashlib
import importlib
import json
import platform
import sys
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

from nomad.datamodel import EntryArchive
from nomad.utils import get_logger

LOGGER = get_logger(__name__)

# Parser directory -> [(dotted parser class, mainfile)]. The classes are given as
# `module:ClassName` strings and imported lazily in generate() for the requested
# parsers only, so snapshotting one parser never drags in another's optional
# dependencies (e.g. MDAnalysis via gromacs/h5md). Keyed by the parser directory
# so selection matches the names `check_additivity.sh` derives from a diff.
PARSERS: dict[str, list[tuple[str, str]]] = {
    'vasp': [
        (
            'nomad_simulation_parsers.parsers.vasp.parser:VASPParser',
            'tests/data/vasp/AgAc_relax/vasprun.xml.relax',
        ),
        (
            'nomad_simulation_parsers.parsers.vasp.parser:VASPParser',
            'tests/data/vasp/AgAc_relax/OUTCAR',
        ),
    ],
    'abinit': [
        (
            'nomad_simulation_parsers.parsers.abinit.parser:AbinitParser',
            'tests/data/abinit/Fe/abinit.out',
        )
    ],
    'gpaw': [
        (
            'nomad_simulation_parsers.parsers.gpaw.parser:GPAWParser',
            'tests/data/gpaw/Fe2.gpw',
        )
    ],
    'octopus': [
        (
            'nomad_simulation_parsers.parsers.octopus.parser:OctopusParser',
            'tests/data/octopus/Fe_spinpol/inp',
        )
    ],
    'ams': [
        (
            'nomad_simulation_parsers.parsers.ams.parser:AMSParser',
            'tests/data/ams/scf/ams.log',
        )
    ],
    'exciting': [
        (
            'nomad_simulation_parsers.parsers.exciting.parser:ExcitingParser',
            'tests/data/exciting/C_minimal/INFO.OUT',
        )
    ],
    'crystal': [
        (
            'nomad_simulation_parsers.parsers.crystal.parser:CrystalParser',
            'tests/data/crystal/single_point/dft/output.out',
        )
    ],
    'fhiaims': [
        (
            'nomad_simulation_parsers.parsers.fhiaims.parser:FHIAimsParser',
            'tests/data/fhiaims/Si_geomopt/out.out',
        )
    ],
    'gromacs': [
        (
            'nomad_simulation_parsers.parsers.gromacs.parser:GromacsParser',
            'tests/data/gromacs/fe_test/md.log',
        )
    ],
    'h5md': [
        (
            'nomad_simulation_parsers.parsers.h5md.parser:H5MDParser',
            'tests/data/h5md/test_traj_openmm_reduced-SOL_5frames_07-10-25.h5',
        )
    ],
    'lammps': [
        (
            'nomad_simulation_parsers.parsers.lammps.parser:LammpsParser',
            'tests/data/lammps/1_xyz_files/log.lammps',
        )
    ],
    'phonopy': [
        (
            'nomad_simulation_parsers.parsers.phonopy.parser:PhonopyParser',
            'tests/data/phonopy/vasp/phonopy.yaml',
        )
    ],
    'quantumespresso': [
        (
            'nomad_simulation_parsers.parsers.quantumespresso.parser'
            ':QuantumEspressoParser',
            'tests/data/quantumespresso/pwscf/TiO2_opt/pw.out',
        )
    ],
    'wannier90': [
        (
            'nomad_simulation_parsers.parsers.wannier90.parser:Wannier90Parser',
            'tests/data/wannier90/lco_mlwf/lco.wout',
        )
    ],
    'orca': [
        (
            'nomad_simulation_parsers.parsers.orca.parser:OrcaParser',
            'tests/data/orca/RI_MP2_water.out',
        )
    ],
    'lobster': [
        (
            'nomad_simulation_parsers.parsers.lobster.parser:LobsterParser',
            'tests/data/lobster/NaCl/lobsterout',
        )
    ],
    # yambo is registered so a future fixture is picked up with no further wiring,
    # but it has no end-to-end fixture in the checkout yet (its tests are
    # unit-only). generate() skips any mainfile that is not present, so this stays
    # dormant until tests/data/yambo/ exists -- replace the placeholder path then.
    'yambo': [
        (
            'nomad_simulation_parsers.parsers.yambo.parser:YamboParser',
            'tests/data/yambo/example/r-example',
        )
    ],
}


def load_parser(dotted: str) -> type:
    """Import a ``module:ClassName`` string into the parser class it names."""
    module_name, _, class_name = dotted.partition(':')
    return getattr(importlib.import_module(module_name), class_name)


def serialize_archive_data(archive: EntryArchive):
    """Convert archive.data to a JSON-serializable dict for comparison."""
    if not archive.data:
        return None
    return archive.data.m_to_dict(with_meta=False)


def _sha256(path: Path) -> str | None:
    """SHA-256 of a file, or None if it does not exist."""
    if not path.is_file():
        return None
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b''):
            digest.update(chunk)
    return digest.hexdigest()


def _package_version(name: str) -> str | None:
    try:
        return version(name)
    except PackageNotFoundError:
        return None


def provenance(mainfiles: list[str]) -> dict:
    """Record what produced this snapshot, so a stale environment or a changed
    fixture cannot masquerade as an identical run.

    Includes the interpreter version, the serializer (``nomad-lab``) version, the
    harness hash (this generator), the ``uv.lock`` hash when present, and a hash
    per fixture that was snapshotted. ``compare_snapshots.py`` reports these for
    context but never diffs them: two refs may legitimately differ here.
    """
    return {
        'python': platform.python_version(),
        'nomad_lab': _package_version('nomad-lab'),
        'harness_sha256': _sha256(Path(__file__)),
        'uv_lock_sha256': _sha256(Path('uv.lock')),
        'fixtures': {mainfile: _sha256(Path(mainfile)) for mainfile in mainfiles},
    }


def generate(parser_dirs: list[str]) -> tuple[dict, list[str]]:
    snapshots = {}
    mainfiles: list[str] = []
    for parser_dir in parser_dirs:
        for dotted, mainfile in PARSERS[parser_dir]:
            if not Path(mainfile).is_file():
                # A registered parser whose fixture is absent from the checkout
                # (yambo until one is added, or a large fixture) is skipped, not
                # a failure -- it simply does not contribute a snapshot.
                print(f'skip {parser_dir}: fixture not present ({mainfile})')
                continue
            try:
                parser_class = load_parser(dotted)
            except ImportError as error:
                # The parser's optional dependencies are not installed in this
                # environment; skip it rather than failing the whole run.
                print(f'skip {parser_dir}: {dotted} unavailable ({error})')
                continue
            archive = EntryArchive()
            parser_class().parse(mainfile, archive, LOGGER)
            key = f'{parser_dir}:{Path(mainfile).name}'
            snapshots[key] = serialize_archive_data(archive)
            mainfiles.append(mainfile)
    return snapshots, mainfiles


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
    snapshots, mainfiles = generate(requested)
    document = {'provenance': provenance(mainfiles), 'snapshots': snapshots}
    Path(out).write_text(json.dumps(document, sort_keys=True, indent=1, default=str))
    print(f'wrote {len(snapshots)} snapshot(s) for {len(requested)} parser(s): {out}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
