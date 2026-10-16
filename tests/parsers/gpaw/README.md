# GPAW parser test matrix

| Feature | Reader | Mapping | Integration |
| --- | --- | --- | --- |
| Mainfile recognition | `test_recognition.py` | not applicable | `Fe2.gpw`, `H2.gpw` |
| Program metadata and units | `.gpw`/`.gpw2` fixtures in `test_reader.py` | controlled source in `test_mapping.py` | `Fe2.gpw` |
| Structure, atomic species, and boundary conditions | `test_reader.py` | not yet isolated | `Fe2.gpw`, `H2.gpw` |
| DFT exchange-correlation method | not yet isolated | controlled source in `test_mapping.py` | `Fe2.gpw` |
| Energies, forces, and SCF state | not yet isolated | `test_mapping.py` | `Fe2.gpw` |
| Eigenvalues, occupations, and gaps | not yet isolated | `test_mapping.py` | `Fe2.gpw`, `H2.gpw` |
| Band structure | not yet isolated | controlled source in `test_mapping.py` | not present in current fixtures |
| Single-point workflow | not yet isolated | `test_mapping.py` | `Fe2.gpw`, `H2.gpw` |
| `.gpw` and `.gpw2` reader compatibility | `test_reader.py` | not applicable | `Fe2.gpw`, `Si_pw.gpw2`, `Si_lcao.gpw2` |
| NOMAD normalization | not applicable | not applicable | `Fe2.gpw` pipeline test |

The mapping tests use a small controlled source object so conversion helpers can
be checked without coupling each assertion to a complete binary fixture.
Complete `.gpw` and `.gpw2` files in `tests/data/gpaw/` remain end-to-end
inputs. Rows marked **not yet isolated** identify lower-layer coverage to add
when those parsing or mapping rules change.
