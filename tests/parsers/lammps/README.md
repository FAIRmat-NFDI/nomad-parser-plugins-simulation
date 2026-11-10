# LAMMPS parser test matrix

The [fixture guide](../../data/lammps/README.md) lists every bundled LAMMPS
directory, its parser mainfile, expected outcome, and intentionally absent
dependencies. It also distinguishes optional NOMAD Central diagnostics from
versioned test data.

| Feature | Recognition | Reader | Archive writer | Integration | Pipeline |
| --- | --- | --- | --- | --- | --- |
| Mainfile recognition | `test_recognition.py` | not applicable | not applicable | not applicable | not applicable |
| LAMMPS commands, units, and auxiliary-file selection | not applicable | `test_reader.py` | `test_units.py` | `test_units.py` | not applicable |
| MD run selection and workflow-method settings | not applicable | `test_reader.py`, `test_md_workflow_method.py` | `test_mapping.py` | `test_integration.py` | representative normalization in `test_pipeline.py` |
| Archive-writer unit conversion | not applicable | not applicable | `test_archive_writer.py` | not applicable | not applicable |
| Trajectory frames, positions, velocities, forces, and cells | not applicable | `test_reader.py` | not applicable | `test_integration.py` | `test_pipeline.py` |
| Input structure when a declared dump is absent | not applicable | `test_reader.py` | not applicable | `test_input_structure.py` | not applicable |
| Molecular topology and particle identity | not applicable | not applicable | not applicable | `test_integration.py` | `test_pipeline.py` |
| Archive and workflow validation and serialization | not applicable | not applicable | not applicable | `test_integration.py` | not applicable |

The fixtures use the existing LAMMPS regression examples under
`tests/data/lammps`. Reader tests use small inline trajectories and injected
parser results, while archive integration tests cover methyl-naphthalene,
hexane, methane XYZ/DCD, and an optional external polymer-melt minimization
example. The
`1_xyz_files` and `2_xyz_files` examples use reduced LJ units: they remain
reader fixtures, but archive parsing rejects them before assigning physical
quantities. `test_units.py` checks that rejection and the `real`/`metal`
velocity scales; the methane XYZ archive covers a physical-unit XYZ input.

LAMMPS populates a `MolecularDynamics` workflow and method for a selected MD
run. `test_md_workflow_method.py` checks command-order selection and direct-parser edge
cases; `test_mapping.py` maps controlled settings to schema quantities without
reading a log. MD integration fixtures compose the shared workflow contract.
The method tests also cover canonical-sampling thermostats, Berendsen pressure
control, anisotropic/tilted Nose–Hoover pressure, and Brownian metadata using
small synthetic logs. These are not yet positive Central-corpus fixtures.
The minimization fixture intentionally does not require an MD workflow; a
geometry-optimization workflow is separate follow-up work. Thermodynamic
`outputs` remain unimplemented.
`test_input_structure.py` checks a source-backed, non-representative input
system for data-only logs, redirected-log `read_data`, missing dumps, and
non-angstrom physical units. It does not claim a final minimized or MD frame.
