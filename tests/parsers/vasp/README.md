# VASP parser test matrix

| Feature | Recognition | Reader | Mapping | Integration | Pipeline |
| --- | --- | --- | --- | --- | --- |
| OUTCAR recognition | `test_recognition.py` |  |  | `AgAc_relax/OUTCAR` |  |
| vasprun.xml recognition | `test_recognition.py` |  |  | `AgAc_relax/vasprun.xml.relax` |  |
| OUTCAR header, parameters, and structure |  | `test_reader.py` |  | `AgAc_relax/OUTCAR` |  |
| vasprun.xml source sections |  | `test_reader.py` |  | `AgAc_relax/vasprun.xml.relax` |  |
| XC functional mapping |  |  | `test_mapping.py` | both archives |  |
| SCF, electronic outputs, and workflow |  |  | `test_mapping.py` | both archives |  |
| Mg static and band vasprun files |  |  |  | `test_integration.py` |  |
| Silicon DOS and band vasprun files |  |  |  | `test_integration.py` |  |
| GW, gamma-only, hybrid, and meta-GGA inputs |  |  |  | `test_integration.py` |  |
| Alternative pseudopotentials |  |  |  | `test_integration.py` |  |
| DFT+U inputs with and without INCAR |  |  |  | `test_integration.py` |  |
| Boolean INCAR spellings |  |  |  | `test_integration.py` |  |
| Malformed time and broken XML handling |  |  |  | `test_integration.py` |  |
| NOMAD normalization compatibility |  |  |  |  | `test_pipeline.py` |

The integration archive assertions include the following numerical settings and
functional checks:

- Silicon GW verifies the 6x6x6 KMesh, representative k-point coordinates,
  weights, and offset.
- Gamma-only, hybrid, and meta-GGA archives verify their KMesh values.
- Hybrid verifies the `HSE06` functional; meta-GGA verifies the mapped `PBE`
  functional exposed by the fixture.

Numerical arrays in the parser tests use the shared `assert_approx` helper from
`tests.parsers.common` rather than direct NumPy assertion calls.

The reader and mapping tests use small, explicit source fragments where possible;
the integration and pipeline tests retain the representative VASP relaxation
fixtures. Changes to a parser rule should update the smallest applicable layer
and add an integration assertion when the public archive contract changes.

The legacy electronicparsers tests also contained helper-level band-path tests
for APIs that are not present in the current VASP parser, and a Cu3PS4 fixture
that is not available in this package's test data. Those cases are intentionally
not copied verbatim; the available band and compressed-vasprun inputs are covered
through the current archive contract instead.
