# LAMMPS parser test matrix

| Feature | Recognition | Reader | Archive writer | Integration | Pipeline |
| --- | --- | --- | --- | --- | --- |
| Mainfile recognition | `test_recognition.py` | not applicable | not applicable | not applicable | not applicable |
| LAMMPS commands, units, and auxiliary-file selection | not applicable | `test_reader.py` | `test_units.py` | `test_units.py` | not applicable |
| MD run selection and method settings | not applicable | `test_reader.py`, `test_method.py` | `test_mapping.py` | `test_integration.py` | representative normalization in `test_pipeline.py` |
| Archive-writer unit conversion | not applicable | not applicable | `test_archive_writer.py` | not applicable | not applicable |
| Trajectory frames, positions, velocities, forces, and cells | not applicable | `test_reader.py` | not applicable | `test_integration.py` | `test_pipeline.py` |
| Molecular topology and particle identity | not applicable | not applicable | not applicable | `test_integration.py` | `test_pipeline.py` |
| Archive and workflow validation and serialization | not applicable | not applicable | not applicable | `test_integration.py` | not applicable |

The fixtures use the existing LAMMPS regression examples under
`tests/data/lammps`. Reader tests use small inline trajectories and injected
parser results, while archive integration tests cover methyl-naphthalene,
hexane, methane XYZ/DCD, and polymer-melt minimization examples. The
`1_xyz_files` and `2_xyz_files` examples use reduced LJ units: they remain
reader fixtures, but archive parsing rejects them before assigning physical
quantities. `test_units.py` checks that rejection and the `real`/`metal`
velocity scales; the methane XYZ archive covers a physical-unit XYZ input.

LAMMPS populates a `MolecularDynamics` workflow and method for a selected MD
run. `test_method.py` checks command-order selection and direct-parser edge
cases; `test_mapping.py` maps controlled settings to schema quantities without
reading a log. MD integration fixtures compose the shared workflow contract.
The minimization fixture intentionally does not require an MD workflow; a
geometry-optimization workflow is separate follow-up work. Thermodynamic
`outputs` remain unimplemented.
