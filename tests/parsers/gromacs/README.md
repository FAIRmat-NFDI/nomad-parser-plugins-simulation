# GROMACS parser test matrix

| Feature | Recognition | Reader and mapping | Archive construction | Integration |
| --- | --- | --- | --- | --- |
| Mainfile recognition | `test_recognition.py` | not applicable | not applicable | not applicable |
| Trajectory configurations, labels, velocities, and cells | not applicable | `test_reader.py` | `test_mapping.py` | `test_integration.py` |
| GROMACS input-parameter transformations | not applicable | `test_reader.py` | `test_mapping.py` | not applicable |
| Thermodynamic outputs and energies | not applicable | `test_reader.py` | `test_mapping.py` | `test_integration.py` |
| Topology bonds and force-field contributions | not applicable | not applicable | `test_mapping.py` | `test_integration.py` |
| Molecular hierarchy and particle-index validity | not applicable | not applicable | not applicable | `test_integration.py` |
| Integrator enum mappings (`water_AA_ENUM_tests`) | not applicable | not applicable | not applicable | `test_integration.py` |
| Free-energy and XVG data | not applicable | `test_reader.py` | not applicable | `test_integration.py` |
| Archive validation and serialization | not applicable | not applicable | not applicable | `test_integration.py` |
| NOMAD normalization pipeline | not applicable | not applicable | not applicable | `test_pipeline.py` |

The fixture set is intentionally composed of the GROMACS regression examples and
larger molecular-dynamics inputs under `tests/data/gromacs`. Unit tests use
small stubs for parser boundaries; integration tests parse the corresponding
GROMACS log, topology, trajectory, energy, and XVG files together. The archive tests also cover the polymer minimization and water integrator regression fixtures.
