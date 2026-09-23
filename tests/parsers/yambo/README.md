# YAMBO parser test matrix

| Feature | Recognition | Reader | Mapping | Integration | Pipeline |
| --- | --- | --- | --- | --- | --- |
| YAMBO mainfile recognition | `test_recognition.py` |  |  | `hBN/r-10b_1Ry_HF_and_locXC_gw0_em1d_ppa` |  |
| Mainfile band-edge and temperature values |  | minimal report fragment in `test_reader.py` |  | hBN report |  |
| Input key/value parsing |  | `test_reader.py` |  |  |  |
| Eigenvalues and highest occupied energy |  |  | `test_mapping.py` | hBN report |  |
| Band-gap mapping |  |  | `test_mapping.py` | hBN report |  |
| NetCDF structure and k-space helpers |  |  | `test_mapping.py` | hBN databases |  |
| GW/QP output archive |  |  |  | `hBN` integration class |  |
| Lifetime output archive |  |  |  | `Aluminum` integration class |  |
| QP/PPA output archive |  |  |  | `LiF` integration class |  |
| GW output archive |  |  |  | `GaSb` integration class |  |
| Minimal setup database archive |  |  |  | `CH4_db_minimal` integration class |  |
| NOMAD normalization compatibility |  |  |  |  | `test_pipeline.py` |

Small reader fragments isolate YAMBO syntax from the complete report fixtures.
The hBN calculation is the representative end-to-end archive because its
report and NetCDF databases are checked in together. The other reports cover
QP, GW, and lifetime output variants. Mapping tests use the public mapper
methods directly so changes to schema transformation are localized.

The integration assertions cover the current `Simulation` schema: program
versions, K-point data, band gaps, eigenvalue shapes and representative values,
and highest-occupied energies. Custom `x_yambo_*` method and calculation
sections without current schema equivalents are not included.

Numeric comparisons follow the units reported by YAMBO: NetCDF lattice vectors
are compared in bohr (`[a.u.]`), while report energies, band gaps, and
eigenvalues are compared in eV. The shared `approx` and `assert_approx` helpers
are used for scalar and array comparisons.
