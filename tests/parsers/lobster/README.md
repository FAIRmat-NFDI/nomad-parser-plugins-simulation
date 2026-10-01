# LOBSTER parser test matrix

| Feature | Reader | Mapping | Integration |
| --- | --- | --- | --- |
| Mainfile recognition | `test_recognition.py` | not applicable | Fe, NaCl, HfV₂, QE/Ni, Si |
| Program metadata, source code, and basis | minimal `lobsterout` in `test_reader.py` | `TestLobsterMainfileParserMapping` | Fe, NaCl, HfV₂, QE/Ni, Si |
| Charge parsing | minimal `CHARGE.lobster` in `test_reader.py` | `TestLobsterCHARGEParserMapping` | Fe, NaCl |
| Structure and atomic species | not yet isolated | not yet isolated | Fe, NaCl, HfV₂, QE/Ni, Si |
| Bond-resolved COHP/COOP/COBI data | minimal `ICOHPLIST.lobster` and `COHPCAR.lobster` in `test_reader.py` | `TestLobsterICOXPLISTParserMapping`, `TestLobsterCOXPCARParserMapping` | Fe, NaCl, HfV₂, Si |
| Density of states | minimal `DOSCAR.lobster` in `test_reader.py` | `TestLobsterDOSCARParserMapping` | Fe, NaCl |
| Spin-polarized output | not yet isolated | not yet isolated | Fe, Si |
| Non-spin-polarized output | not yet isolated | not yet isolated | NaCl, HfV₂ |
| Quantum Espresso structure and metadata | not yet isolated | not yet isolated | QE/Ni |
| Partial, backup-structure, and compressed output handling | not yet isolated | not applicable | HfV₂ backup structure, compressed Si |
| Unimplemented schema fields | not applicable | not applicable | TODO assertions in `test_integration.py` |
| Workflow mapping | not applicable | not applicable | not applicable (LOBSTER has no workflow section) |
| NOMAD normalization | not applicable | not applicable | Fe pipeline test |

The minimal reader sources are embedded in `test_reader.py` because each
fragment demonstrates one LOBSTER syntax feature and is small enough to review
inline. Controlled sources in `test_mapping.py` isolate helper behavior from
file extraction. Complete calculation fixtures in `tests/data/lobster/` remain
end-to-end inputs, and the original scientific assertions are integrated in
`test_integration.py`. Rows marked **not yet isolated** identify lower-layer
coverage that can be added when those parsing or mapping rules change.
