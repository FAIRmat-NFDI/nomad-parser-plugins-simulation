# ORCA parser test matrix

| Feature | Reader | Mapping | Integration |
| --- | --- | --- | --- |
| Mainfile recognition | `test_recognition.py` | not applicable | not applicable |
| Structure and atomic species | `CARTESIAN COORDINATES` fragment in `test_reader.py` | not yet isolated (`get_atoms` is a root aggregator) | `RI_MP2_water` |
| SCF settings (charge, multiplicity, reference) | SCF-settings fragment in `test_reader.py` | `TestDFTMapping` (reference form) | `RI_MP2_water`, `dft-print-MOs` |
| DFT exchange-correlation | not yet isolated | `TestDFTMapping` | `dft-print-MOs` |
| HF reference | not yet isolated | not yet isolated (`get_hf_methods`) | `dlpno-coupled-cluster` |
| MP2 / perturbation (SCS, local correlation) | not yet isolated | `TestPerturbationMapping` | `RI_MP2_water`, `dlpno-coupled-cluster` |
| Coupled cluster (order, (T), DLPNO, thresholds) | not yet isolated | `TestCoupledClusterMapping` | `dlpno-coupled-cluster` |
| Multireference CASSCF / CASCI (active space, states) | not yet isolated | `TestMultireferenceMapping` | `CoPc_CASCI_QD` (nightly) |
| Multireference PT (NEVPT2 name) | not yet isolated | `TestMultireferenceMapping` | `CoPc_CASCI_QD` (nightly) |
| Scalar relativity (DKH) | not yet isolated | `TestRelativityMapping` | `CoPc_CASCI_QD` (nightly) |
| Basis sets (role assignment per method) | not yet isolated | `TestBasisSetMapping` | `RI_MP2_water`, `CoPc_CASCI_QD` (nightly) |
| Orbital energies and occupations | `ORBITAL ENERGIES` fragment in `test_reader.py` | not yet isolated | `RI_MP2_water` |
| Molecular-orbital coefficients (HDF5) | not yet isolated | not yet isolated | `dft-print-MOs` |
| Derived HOMO-LUMO gap | not applicable | not applicable | `RI_MP2_water` pipeline |
| HF -> CC serial workflow | source inputs not isolated | not yet isolated | `dlpno-coupled-cluster` |
| Reader error handling and line-ending invariance | truncated / malformed / LF-CRLF fragments in `test_reader.py` | not applicable | not applicable |
| NOMAD normalization | not applicable | not applicable | `RI_MP2_water` pipeline |

The minimal reader source is embedded in `test_reader.py` because each fragment
demonstrates one block and is small enough to review inline. Complete
calculation fixtures in `tests/data/orca/` remain end-to-end inputs. Rows marked
**not yet isolated** identify the lower-layer coverage still needed when those
parsing or mapping rules change.

Notes:

- `CoPc_CASCI_QD.out` is ~19 MB; its integration class is marked
  `@pytest.mark.large_fixture` and runs only in the nightly gate.
- `dft-print-MOs.out` stores the MO coefficient matrix through the HDF5 backend,
  so `TestDFTPrintMOs` skips the two generic serialization round-trips (they
  detach from the upload-backed context).
- `CoPc_CASSCF_SA.out` and `CoPc_CASSCF_SS.out` (~8.8 MB each) are present in
  `tests/data/orca/` but not yet referenced by a test.
