# ABINIT parser test matrix

| Feature | Reader | Mapping | Integration |
| --- | --- | --- | --- |
| Mainfile recognition | `test_recognition.py` | not applicable | Fe |
| Program metadata and datetime | minimal ABINIT output in `test_reader.py` | controlled source | Fe, Si |
| Structure and atomic species | input variables in `test_reader.py` | `TestOutMapping` | Fe, Si, ZrO₂ GW |
| DFT exchange-correlation method | not yet isolated | `TestOutMapping` | Fe, Si |
| Energies, forces, and SCF history | minimal SCF rows in `test_reader.py` | `TestOutMapping` | Fe, H₂, Si |
| Eigenvalues, occupations, and gaps | results block in `test_reader.py` | `TestOutMapping` | Fe, ZrO₂ GW |
| Density of states | DOS table in `test_reader.py` | `TestDosMapping` | Fe |
| Single-point workflow | input variables not yet isolated | `TestWorkflowMapping` | Fe, Si, ZrO₂ GW |
| Geometry-optimization workflow | input variables not yet isolated | `TestWorkflowMapping` | H₂ |
| Partial and compressed output handling | truncated output and gzip source in `test_reader.py` | not applicable | truncated-output scenario |
| NOMAD normalization | not applicable | not applicable | Fe pipeline test |

The minimal reader sources are embedded in `test_reader.py` because each
fragment demonstrates one ABINIT syntax feature and is small enough to review
inline. Controlled sources in `test_mapping.py` isolate mapping behavior from
file extraction. Complete calculation fixtures in `tests/data/abinit/` remain
end-to-end inputs. Rows marked **not yet isolated** identify lower-layer
coverage that can be added when those parsing rules change.
