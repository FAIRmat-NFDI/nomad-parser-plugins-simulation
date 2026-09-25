# Quantum ESPRESSO parser test guide

The suite follows the parser-testing pattern used by the other simulation
plugins: readers test source parsing, mapping tests exercise parser-to-schema
transformations, and integration tests validate complete `EntryArchive`
objects.

| Feature | Reader | Mapping | Integration |
| --- | --- | --- | --- |
| Mainfile recognition | `test_recognition.py` | not applicable | text |
| Program metadata and module dispatch | `test_reader.py` | `test_mapping.py` | text, XML, EPW, phonon, XSpectra, GIPAW |
| Structure, lattice, periodicity, and XC | `test_reader.py` | `test_mapping.py` | PWSCF text/XML `TiO2_opt` |
| Energies, forces, SCF, eigenvalues, and bands | `test_reader.py` | `test_mapping.py` | PWSCF text/XML `TiO2_opt` |
| GIPAW mapper methods | `test_reader.py` | `test_mapping.py` | GIPAW text/XML fixtures |
| PWSCF mapper methods | `test_reader.py` | `test_mapping.py` | PWSCF text `TiO2_opt` |
| Density of states | `test_reader.py` | not yet isolated | `W_dos` |
| Geometry optimization workflow | not yet isolated | `test_mapping.py` | PWSCF text/XML `TiO2_opt` |
| Additional QE module readers | EPW, phonon, XSpectra, GIPAW | not applicable | dedicated integration classes |
| Reader error handling and line-ending invariance | not yet isolated | not applicable | not applicable |
| NOMAD normalization | not applicable | not applicable | text `TiO2_opt` pipeline test |

## Integration fixtures

`test_integration.py` uses dedicated fixtures from `conftest.py` for PWSCF text
and XML, DOS, EPW, phonon, XSpectra, and all GIPAW NMR, EFG, hyperfine, and
delta-g text/XML cases. Every module class inherits the common
`QuantumEspressoIntegrationSuite`; fixtures that do not provide particle-state
identities skip only the identity-specific contract.

`test_pipeline.py` covers normalization of the PWSCF text archive.
