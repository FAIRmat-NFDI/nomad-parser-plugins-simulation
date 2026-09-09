# FHI-aims parser test matrix

| Feature | Reader | Mapping | Integration |
| --- | --- | --- | --- |
| Mainfile recognition | `test_recognition.py` | not applicable | `Si_geomopt` |
| Geometry and control inputs | `test_reader.py` | not yet isolated | `Si_geomopt` |
| Structure and periodic cell | representative output reader | not yet isolated | `Si_geomopt` |
| SCF convergence criteria | representative output reader | `test_mapping.py` | `Si_geomopt` |
| K-space mesh and offset | `test_reader.py` | `test_mapping.py` | `Si_geomopt` |
| Eigenvalues and band structures | representative output reader | `test_mapping.py` | `Si_geomopt` |
| Density of states | not yet isolated | `test_mapping.py` | `ClNa_dos` |
| Geometry optimization workflow | representative output reader | not yet isolated | `Si_geomopt` |
| Negative and line-ending handling | minimal `aims.out` in `test_reader.py` | not applicable | not applicable |
| NOMAD normalization | not applicable | not applicable | `Si_geomopt` pipeline test |

The complete calculation fixtures in `tests/data/fhiaims/` remain end-to-end
inputs. Rows marked **not yet isolated** identify lower-layer coverage to add when
those parsing or mapping rules change.

The integration fixtures are declared in `conftest.py` and each is mapped to an
explicit integration class in `test_integration.py`. This covers the legacy
electronicparsers cases (spin-polarized SCF and bands, silicon DOS/bands, DOS, MD,
hybrid, DFT+U, GW, and native tiers).

Where the new `nomad-simulations` schema has an equivalent quantity, the fixture
classes also assert the mapped method, k-mesh, DOS, band, eigenvalue, trajectory,
and workflow values. Legacy electronicparsers-only quantities are not asserted
against unrelated fields.
