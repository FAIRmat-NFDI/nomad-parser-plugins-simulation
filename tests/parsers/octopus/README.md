# octopus parser test matrix

| Feature | Reader | Mapping | Integration |
| --- | --- | --- | --- |
| Mainfile recognition | `test_recognition.py` | not applicable | not applicable |
| Input and auxiliary-file parsing | `test_reader.py` | not yet isolated | `Si_scf`, `Fe_spinpol` |
| Structure, lattice, and periodicity | not yet isolated | not yet isolated | `test_integration.py` |
| Exchange-correlation functional mapping | not yet isolated | `test_mapping.py` | `Si_scf`, `Fe_spinpol` |
| Energies and forces | not yet isolated | `test_mapping.py` | `test_integration.py` |
| Eigenvalues, occupations, and band values | `test_reader.py` | `test_mapping.py` | `Si_scf` |
| Spin-polarised output | not yet isolated | not yet isolated | `Fe_spinpol` |
| Single-point workflow | source inputs not isolated | not yet isolated | `Si_scf`, `Fe_spinpol` |
| Reader error handling and line-ending invariance | not yet isolated | not applicable | not applicable |
| NOMAD normalization | not applicable | not applicable | `Si_scf` pipeline test |

The complete calculation directories in `tests/data/octopus/` exercise the
multi-file parser (`stdout.txt`, `inp`, `exec/`, and `static/`). Rows marked
**not yet isolated** identify lower-layer coverage to add when those parsing or
mapping rules change.

The integration tests assert program metadata, DFT/XC mapping, silicon and
iron geometry, lattice vectors, energies, forces, occupations, and band
values. K-mesh and per-k-point coordinates are asserted at the `InfoParser`
reader level in `test_reader.py`; SCF timing remains uncovered because it is
not currently exposed by the simulation archive schema.
