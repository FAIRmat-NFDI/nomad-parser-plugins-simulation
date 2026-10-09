# phonopy parser test matrix

| Feature | Recognition | Calculator | Mapping | Archive construction | Integration |
| --- | --- | --- | --- | --- | --- |
| Mainfile recognition | `test_recognition.py` | not applicable | not applicable | not applicable | not applicable |
| K-path generation and calculator setup | not applicable | `test_calculator.py` | not applicable | not applicable | not applicable |
| Non-canonical hexagonal band path | not applicable | not applicable | not applicable | not applicable | `test_integration.py` |
| Unit and supercell structures | not applicable | not applicable | `test_mapping.py` | `test_archive_writer.py` | `test_integration.py` |
| Force-constant sidecar loading | not applicable | not applicable | not applicable | `test_archive_writer.py` | not applicable |
| Band structures, DOS, and thermodynamics | not applicable | `test_calculator.py` | `test_mapping.py` | `test_archive_writer.py` | `test_integration.py` |
| Archive validation and serialization | not applicable | not applicable | `test_mapping.py` | not applicable | `test_integration.py` |
| NOMAD normalization | not applicable | not applicable | not applicable | not applicable | `test_pipeline.py` |

The complete fixture is intentionally kept as an end-to-end input because
phonopy reads the structure and force constants together. The focused tests use
small phonopy objects or controlled property results at the calculator, mapping,
and archive-construction boundaries.
