# LAMMPS parser test matrix

| Feature | Recognition | Reader | Archive writer | Integration | Pipeline |
| --- | --- | --- | --- | --- | --- |
| Mainfile recognition | `test_recognition.py` | not applicable | not applicable | not applicable | not applicable |
| LAMMPS commands, units, and auxiliary-file selection | not applicable | `test_reader.py` | not applicable | not applicable | not applicable |
| Archive-writer unit conversion | not applicable | not applicable | `test_archive_writer.py` | not applicable | not applicable |
| Trajectory frames, positions, velocities, forces, and cells | not applicable | `test_reader.py` | not applicable | `test_integration.py` | `test_pipeline.py` |
| Molecular topology and particle identity | not applicable | not applicable | not applicable | `test_integration.py` | `test_pipeline.py` |
| Archive validation and serialization | not applicable | not applicable | not applicable | `test_integration.py` | not applicable |

The fixtures use the existing LAMMPS regression examples under
`tests/data/lammps`. Reader tests use small inline trajectories and injected
parser results, while archive integration tests cover methyl-naphthalene,
1-xyz, hexane, methane XYZ/DCD, and polymer-melt minimization examples.

LAMMPS currently does not populate thermodynamic `outputs` or a workflow, so
those contracts are intentionally not listed.
