# exciting parser test matrix

| Feature | Reader | Mapping | Integration |
| --- | --- | --- | --- |
| Mainfile recognition | `test_recognition.py` | not applicable | `C_minimal` |
| Program metadata and initialization lattice | minimal `INFO.OUT` in `test_reader.py` (version hash) | not yet isolated | `C_minimal` |
| Structure and atomic species | cell volume, atom count and species mass in `test_reader.py` | not yet isolated | `C_minimal`, `C_gs` |
| DFT exchange-correlation method | not yet isolated | `test_maps_xc_type_to_functional` | `C_minimal`, `PbI_hybrids` |
| Spin treatment | minimal `INFO.OUT` in `test_reader.py` | `test_maps_spin_treatment` | `C_minimal`, `PbI_hybrids` |
| Smearing kind and width | `test_extracts_method_and_system_settings` | `test_maps_smearing_kind` | `C_minimal`, `PbI_hybrids` |
| k-point mesh | `test_extracts_method_and_system_settings` | `test_maps_k_mesh` | `C_minimal`, `PbI_hybrids` |
| APW basis-set cutoff (`rgkmax`) | `test_extracts_method_and_system_settings` | not applicable (direct path) | `C_minimal` |
| Energies and SCF convergence history | charge distance in `test_reader.py` | `test_maps_scf_convergence_quantities_with_units` | `C_minimal`, `C_gs` |
| Eigenvalues, occupations, and gaps | minimal `EIGVAL.OUT` in `test_reader.py` | `test_maps_eigenvalues_to_band_gaps` | `C_minimal` |
| Band structure | not yet isolated | `test_maps_k_path` | `C_minimal` |
| Density of states | not yet isolated | `test_maps_energy_axis_and_values` | `CeO_dos` |
| Single-point workflow | source inputs not isolated | not yet isolated | `C_minimal`, `C_gs`, `CeO_dos` |
| Geometry-optimization workflow | source inputs not isolated | not yet isolated | `GaO_sodium`, `GaO_strucopt` |
| Reader error handling and line-ending invariance | truncated/malformed and LF/CRLF `INFO.OUT` in `test_reader.py` | not applicable | not applicable |
| NOMAD normalization | not applicable | not applicable | `C_minimal` pipeline test |

The minimal reader source is embedded in `test_reader.py` because each fragment
demonstrates one syntax feature and is small enough to review inline. Complete
calculation fixtures in `tests/data/exciting/` remain end-to-end inputs. Rows
marked **not yet isolated** identify the lower-layer coverage still needed when
those parsing or mapping rules change.
