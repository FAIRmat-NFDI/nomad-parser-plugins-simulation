# Wannier90 parser test matrix

| Feature | Recognition | Reader | Mapping | Integration | Pipeline |
| --- | --- | --- | --- | --- | --- |
| `.wout` recognition and DFT child selection | `test_recognition.py` |  |  | `lco.wout` |  |
| `.wout` version, structure, lattice, k-mesh, and run settings |  | `test_reader.py` | `test_mapping.py` | `lco.wout` |  |
| `.win` projections, Fermi energy, branches, and orbital states |  | `test_reader.py` | `test_mapping.py` (`TestWInTextParser`) | `lco.win` |  |
| Wannier orbital quantum-number mapping |  |  | `test_mapping.py` (`TestWInTextParser`) |  |  |
| `_hr.dat` degeneracies and complex hopping rows |  | `test_reader.py` | `test_mapping.py` | `lco_hr.dat` |  |
| DOS and band array mapping |  |  | `test_mapping.py` | `lco-dos.dat`, `lco_band.dat` |  |
| structure, exact raw positions/lattice, bands, method, and workflow |  |  |  | `test_integration.py` |  |
| NOMAD normalization compatibility |  |  |  |  | `test_pipeline.py` |

The integration fixture is the existing La-Cu-O maximally localized Wannier
calculation. The `lco_mlwf_archive` fixture parses the complete calculation.
Reader tests use small inline mock texts with the
format-specific readers, while mapping tests call conversion helpers with
controlled arrays and use the shared `approx`/`assert_approx` helpers.

Integration tests compare positions and lattice vectors against `lco.wout` and
representative first/last band energies against `lco_band.dat`. The current
archive contract does not yet populate `particle_states` or an explicit band-gap
section for Wannier90.
