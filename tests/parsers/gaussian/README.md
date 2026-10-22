# Gaussian parser test matrix

| Feature | Reader | Mapping | Integration |
| --- | --- | --- | --- |
| Mainfile recognition | `test_recognition.py` | not applicable | not applicable |
| Program metadata and route section | inline mock | `test_maps_program_metadata` | `Al_scf` |
| Atomic positions and species | inline mock | `test_maps_atomic_numbers_and_positions` | all fixtures |
| DFT functional | route section | `test_maps_only_dft_routes` | `Al_scf`, `Al_multistep`, `CHO_freq` |
| SCF energy | inline mock | `test_maps_energy_and_forces` | `Al_scf` |
| Geometry optimization and forces | inline mock | `test_maps_energy_and_forces` | `Al_multistep`, `CHO_freq` |
| Frequencies | inline mock | not yet mapped to the public schema | `CHO_freq` |
| MP and coupled-cluster energies | inline mock | DFT exclusion | `NO_mp` |
| NOMAD normalization | not applicable | not applicable | `Al_scf` pipeline |

The complete calculation fixtures in `tests/data/gaussian/` are copied from the
legacy `electronicparsers` Gaussian tests and are used only by integration and
pipeline tests. Reader tests use minimal inline mock text, while integration
tests assert the migrated public simulation schema.
