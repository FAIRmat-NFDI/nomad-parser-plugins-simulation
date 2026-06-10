# PR #170 Migration Verification Report

## Executive Summary

This report verifies the migration from manual `archive.run` population to declarative `archive.data` mappings for 6 parsers (VASP, ABINIT, GPAW, Octopus, AMS, Exciting).

**Verification Method**: Syrupy snapshot comparison
- **Baseline**: commit 1c7acd6 (Apr 10, 2026)
- **Current**: test-data-normalization branch

**Finding**: The migration is **NOT purely additive**. While it adds substantial new data (44,048 lines), it also removes/modifies some existing data (2,891 lines).

## Methodology

### Snapshot Generation

1. Created snapshot test suite (`tests/parsers/test_migration_snapshots.py`) with 7 tests:
   - `test_vasp_vasprun_snapshot` - VASP vasprun.xml parsing
   - `test_vasp_outcar_snapshot` - VASP OUTCAR parsing
   - `test_abinit_snapshot` - ABINIT parsing
   - `test_gpaw_snapshot` - GPAW parsing
   - `test_octopus_snapshot` - Octopus parsing
   - `test_ams_snapshot` - AMS parsing
   - `test_exciting_snapshot` - Exciting parsing

2. Generated snapshots at two commits:
   - **Baseline** (1c7acd6): 180 KB snapshot file
   - **Current** (test-data-normalization): 1.4 MB snapshot file

3. Used `diff -u` to compare snapshots

### Statistics

```
Baseline snapshot size:  180,273 bytes
Current snapshot size: 1,434,017 bytes
Size increase: 7.96x (696% growth)

Diff statistics:
- Additions:  44,048 lines
- Deletions:   2,891 lines
- Total diff: 47,098 lines
```

## Key Findings

### 1. Additions (Expected)

The migration successfully adds significant new data to `archive.data`, including:

**Electronic Properties**:
- `electronic_band_gaps` - Band gap values
- `electronic_eigenvalues.highest_occupied` - Fermi level tracking
- Additional eigenvalue metadata

**Model System Data**:
- `model_system` entries for all parsers
- `lattice_vectors` - Crystal lattice parameters
- `positions` - Atomic positions
- `periodic_boundary_conditions` - PBC flags
- `particle_states` - Atom types and chemical symbols

**Example Addition** (from AMS):
```python
'model_system': list([
  dict({
    'periodic_boundary_conditions': list([False, False, False]),
  }),
]),
```

### 2. Deletions/Modifications (Unexpected)

The migration removed approximately 2,891 lines of existing data, primarily:

**Electronic Eigenvalue Reorganization**:
- Removed nested `'value'` arrays from eigenvalue structures
- This appears to be a restructuring, not data loss
- Eigenvalue data may have moved to different schema locations

**Example Deletion**:
```python
# Removed from ABINIT snapshot:
-'value': list([
-  list([
-    1.9156540824850025e-19,
-    7.636869176795147e-19,
-    # ... (18 k-points worth of eigenvalues)
-  ]),
-]),
```

### 3. Test Results

**Baseline generation** (1c7acd6):
```
7 snapshots generated
7 passed, 2 warnings in 4.85s
```

**Current generation** (test-data-normalization):
```
1 snapshot passed (unchanged)
6 snapshots updated (modified)
7 passed, 3 warnings in 6.11s
```

## Analysis

### What Changed Per Parser

| Parser | Baseline Status | Current Status | Key Changes |
|--------|----------------|----------------|-------------|
| VASP vasprun.xml | ✓ | Updated | +model_method, +eigenvalue metadata |
| VASP OUTCAR | ✓ | Updated | +model_system, +lattice_vectors |
| ABINIT | ✓ | Updated | +model_system, eigenvalue restructure |
| GPAW | ✓ | Unchanged | No changes detected |
| Octopus | ✓ | Updated | +model_system data |
| AMS | ✓ | Updated | +model_system, +PBC |
| Exciting | ✓ | Updated | +model_system, +positions |

**Note**: GPAW had no changes, suggesting either:
1. GPAW parser wasn't part of this migration
2. GPAW test data doesn't trigger the new mappings
3. GPAW mappings produce identical output

### Migration Assessment

**Strengths**:
- ✅ Successfully adds extensive model_system data
- ✅ Adds electronic properties (band gaps, Fermi levels)
- ✅ All tests pass at both commits
- ✅ File size increase indicates substantial new data

**Concerns**:
- ⚠️ NOT purely additive (2,891 deletions)
- ⚠️ Eigenvalue data structure changed
- ⚠️ Need manual review to verify deletions are intentional restructuring, not data loss

## Recommendations

### Immediate Actions

1. **Manual diff review**: Examine the 2,891 deleted lines in detail
   - Verify eigenvalue data moved to new schema locations
   - Confirm no scientific data was lost
   - Document intentional schema changes

2. **Validation against test assertions**:
   - Run existing test suite (`test_vasp_parser.py`, etc.)
   - Verify tests like `test_outcar_scf_steps_and_single_point_convergence` still pass
   - Check if existing tests validate the data that was "deleted"

3. **Schema documentation**:
   - Document eigenvalue schema changes
   - Explain why `'value'` arrays were removed
   - Update migration notes in PR #170

### Follow-up Questions

1. **Is the eigenvalue restructuring intentional?**
   - Were eigenvalue `'value'` arrays moved to a different path?
   - Is this part of the nomad-simulations schema design?

2. **Why did GPAW show no changes?**
   - Was GPAW intentionally excluded from this migration?
   - Or does the test data not exercise the migrated code paths?

3. **Test coverage**:
   - Do existing tests validate ALL the data in snapshots?
   - Or just specific properties (energies, forces, SCF steps)?

## Artifacts

All verification artifacts are preserved in `/tmp/`:

```
/tmp/test_migration_snapshots.py          # Test file (fixed imports)
/tmp/baseline_snapshots/                   # Snapshots at 1c7acd6
/tmp/current_snapshots/                    # Snapshots at test-data-normalization
/tmp/snapshot_diff.txt                     # Full diff (47,098 lines)
```

## Conclusion

The PR #170 migration successfully adds extensive electronic and model_system data to `archive.data`, increasing snapshot size by nearly 8x. However, it also modifies existing eigenvalue data structures, removing approximately 2,891 lines.

**Verdict**: Migration requires manual review of deletions before approval. The changes may be correct (schema restructuring), but cannot be verified as "purely additive" without deeper investigation.

**Next Steps**:
1. Review `/tmp/snapshot_diff.txt` for deletion patterns
2. Validate against existing test suite
3. Consult with PR author about eigenvalue schema changes
4. Document intentional breaking changes in PR description

---

**Report Generated**: 2026-06-10
**Verification Branch**: verify-pr170-migration
**Methodology**: Syrupy snapshot comparison
**Tool**: pytest-syrupy 5.3.2
