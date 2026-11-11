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
| Thermodynamic records and output mapping | not applicable | `test_reader.py`, `test_outputs.py` | not applicable | `test_outputs.py` | log-only round trip and normalization in `test_outputs.py` |
| Molecular topology and particle identity | not applicable | not applicable | not applicable | `test_integration.py` | `test_pipeline.py` |
| Source-backed force-field settings and typed interaction coefficients | not applicable | `test_reader.py` | `test_force_field.py` | hexane fixture in `test_force_field.py` | `test_pipeline.py` |
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
geometry-optimization workflow is separate follow-up work. `test_outputs.py`
now verifies ordered custom/multi thermo records, temperature, energy,
pressure, volume, density and enthalpy units, elapsed CPU seconds, and
force-bearing dump frames, real-fixture output counts, sparse links to
genuinely matching dump steps, log-only serialization, and non-MD
`WorkflowOutputs`. Dump-backed systems record source `time_step` values for
exact, unique output matching. For DCD/XYZ, the reader reports frame ordinals
rather than source LAMMPS steps, so the parser leaves those links unset. The
shared normalizer can still apply its positional fallback when the two series
have equal length; those inferred links are not source-step verified. These
typed thermodynamic subsections and
conservative normalization of sparse output/system references require the
`lammps-thermo-output-properties` branch in `nomad-simulations`. Unrecognized
thermo columns remain explicitly warned about; extra force-only dump frames
without a matching thermo step are not yet promoted to separate outputs.
`test_force_field.py` checks the bundled hexane log/data/dump for separate
Lennard-Jones, harmonic bond/angle, and OPLS dihedral types. Its interaction
rows use zero-based trajectory indices, not raw LAMMPS atom IDs. Small
source-like tests cover data-file coefficients, explicit log overrides,
pair-specific cutoffs, and unsupported hybrid styles. For supported atom
styles, it also checks effective masses and partial charges against exact
particle references, retaining distinct charge groups within one LAMMPS
type. The mapped subset does not include other bonded/pair styles, wildcard
coefficients, or unsupported mixed-style parameterizations.
`test_input_structure.py` checks a source-backed, non-representative input
system for data-only logs, redirected-log `read_data`, missing dumps, and
non-angstrom physical units. It does not claim a final minimized or MD frame.
