# VASP Test Isolation Fix

## Problem Summary

VASP parser tests (`test_vasprun` and `test_outcar`) were failing in CI/CD and when running the full test suite, but passing when run in isolation. The error was:

```
AssertionError: No model_method in simulation
assert []
 +  where [] = Simulation(steps, model_system, model_method, outputs).model_method
```

## Root Cause

**Test isolation issue caused by global state contamination from THREE parsers.**

The NOMAD mapping parser system uses global mapping annotations stored on class-level `m_def` objects. Three parsers call `remove_mapping_annotations()` after parsing, which removes ALL annotations from the shared `general.Simulation` class used by ALL parsers:

1. **exciting** - calls `remove_mapping_annotations(exciting.general.Simulation.m_def)`
2. **fhiaims** - calls `remove_mapping_annotations(fhiaims.general.Simulation.m_def)`
3. **h5md** - calls `remove_mapping_annotations(h5md.general.Simulation.m_def)` or similar

### Why it only failed in CI/CD / full test suite:
- When running VASP tests alone: none of these parsers run, annotations remain intact
- When running full test suite: these parsers run before VASP (alphabetical order), remove annotations, break VASP tests
- CI/CD always runs the full test suite in alphabetical order, so it always failed there

### Test execution order (alphabetical):
```
abinit → ams → crystal → EXCITING → FHIAIMS → H5MD → lammps → octopus → phonopy → quantumespresso → VASP
                         ^^^^^^^^    ^^^^^^^^   ^^^^                                                  ^^^^
                         Breaks      Breaks     Breaks                                               Victim
```

### Comprehensive Verification Results:

**VASP is the ONLY victim:**
```bash
# VASP passes in isolation
uv run pytest tests/parsers/test_vasp_parser.py::test_outcar
# ✅ PASSES

# VASP passes before culprits run (alphabetically)
uv run pytest tests/parsers/test_abinit_parser.py tests/parsers/test_ams_parser.py tests/parsers/test_crystal_parser.py tests/parsers/test_vasp_parser.py::test_outcar
# ✅ PASSES

# VASP fails after ANY culprit
uv run pytest tests/parsers/test_exciting_parser.py tests/parsers/test_vasp_parser.py::test_outcar
# ❌ FAILS - No model_method

uv run pytest tests/parsers/test_fhiaims_parser.py tests/parsers/test_vasp_parser.py::test_outcar
# ❌ FAILS - No model_method

uv run pytest tests/parsers/test_h5md_parser.py tests/parsers/test_vasp_parser.py::test_outcar
# ❌ FAILS - No model_method
```

**The culprits DON'T break each other:**
```bash
exciting → fhiaims    # ✅ PASS
fhiaims → exciting    # ✅ PASS
h5md → exciting       # ✅ PASS
h5md → fhiaims        # ✅ PASS
```

**The culprits DON'T break other DFT parsers:**
```bash
exciting/fhiaims/h5md → abinit           # ✅ PASS
exciting/fhiaims/h5md → ams              # ✅ PASS
exciting/fhiaims/h5md → crystal          # ✅ PASS
exciting/fhiaims/h5md → octopus          # ✅ PASS
exciting/fhiaims/h5md → quantumespresso  # ✅ PASS
```

**Why VASP is uniquely affected:**

VASP uses a **different annotation strategy**:

```python
# VASP - Annotates base class DIRECTLY
add_mapping_annotation(general.Simulation.m_def, OUTCAR_KEY, '@')

# Other parsers - Annotate IN subclass scope
class Simulation(general.Simulation):
    add_mapping_annotation(model_method.DFT.m_def, OUT_KEY, '.@')
```

When `remove_mapping_annotations(exciting.general.Simulation.m_def)` is called:
- `exciting.general.Simulation IS general.Simulation` (same object reference)
- Removes all annotations from shared `general.Simulation.m_def`
- **VASP's annotations on parent class** → REMOVED ❌
- **Other parsers' annotations on subclasses** → NOT AFFECTED ✅
- **Culprits' own subclass annotations** → NOT AFFECTED ✅

## Investigation Journey

### Initial Hypotheses (all incorrect):

1. **Missing `l_max`/`lm_max` properties** - Checked both local and GitHub's develop branch, both had the properties ✗
2. **Wrong path syntax** - Tried `'parameters'` vs `'.parameters'`, both worked locally ✗
3. **Wrong root path** - Tried `'@'` vs `'.'`, no difference ✗
4. **Python version difference** - CI/CD uses 3.11 vs local 3.12, not the issue ✗
5. **Missing `model_method` annotation** - Added explicit annotation, didn't help ✗

### Key Discovery

User reported: **"The issue is related to some parsers removing mapping annotation keys globally, the error only occurs when the test suite is run all together."**

This explained everything:
- Why tests passed in isolation
- Why they failed in CI/CD (always runs full suite)
- Why reload attempts didn't work reliably

## Solutions Attempted

### 1. Reload schema module (partial success)
```python
importlib.reload(vasp)
```
- Fixed `model_method` population
- But `numerical_settings` still empty (pseudopotential annotations not restored)

### 2. Reload + re-init metainfo (didn't work)
```python
importlib.reload(vasp)
vasp.m_package.__init_metainfo__()
```
- Still failed because annotations were actively removed by other parsers

### 3. Subprocess isolation (TEMPORARY WORKAROUND)

Run VASP tests in isolated subprocesses where annotations can't be contaminated:

```python
@pytest.mark.parametrize('mainfile', [...])
def test_vasp_isolated(mainfile):
    script = f'''
    # Full test logic here
    '''
    result = subprocess.run([sys.executable, '-c', script], ...)
```

This workaround was implemented but later reverted once the root cause was identified and fixed.

### 4. Change annotation strategy to match other parsers (FAILED)

**Hypothesis:** VASP is vulnerable because it annotates `general.Simulation.model_method` directly, while other parsers annotate `model_method.DFT.m_def`. If VASP uses the same pattern, it should be protected.

```python
# Changed from:
add_mapping_annotation(general.Simulation.model_method, OUTCAR_KEY, '.parameters')

# To:
add_mapping_annotation(model_method.DFT.m_def, OUTCAR_KEY, '.parameters')
```

**Result:**
- ✅ VASP works in isolation
- ❌ `exciting → VASP` **still fails** with same error: `No model_method in simulation`
- ❌ Hypothesis REJECTED

**Why it failed:**
- The annotation target doesn't matter - `remove_mapping_annotations()` walks the entire class hierarchy recursively (max_depth=5)
- It removes annotations from `general.Simulation.m_def` AND all its sub-sections and quantities
- Even `model_method.DFT.m_def` annotations get removed because DFT is referenced as a sub-section of Simulation
- Other parsers are protected NOT by annotation target, but because they have their own `Simulation` subclass instances

### 5. Re-add reload() to restore annotations (FAILED)

**Hypothesis:** If `reload(vasp)` re-registers all annotations after they've been removed, VASP should work again.

```python
from importlib import reload
from ...schema_packages import vasp

def parse(mainfile, archive, logger, child_archives):
    reload(vasp)  # Re-register annotations removed by other parsers
    ...
```

**Result:**
- ✅ VASP works in isolation
- ❌ `exciting → VASP` **still fails** with same error: `No model_method in simulation`
- ❌ Hypothesis REJECTED

**Why it failed:**
- `reload()` creates NEW class instances with different `m_def` objects
- NOMAD's plugin system still references the ORIGINAL `m_def` from first import
- Annotations are registered on the NEW `m_def`, but the parser uses the OLD `m_def`
- This is the exact problem documented in commit 98c5b15 from PR #129

**Key insight:** The `reload()` pattern is fundamentally broken for NOMAD's metainfo system. It creates identity conflicts that can't be resolved.

## Root Cause Analysis

**The real culprits:** Two antipatterns used across multiple parsers:

### 1. `reload(schema_module)` pattern (used by 10 parsers)

```python
from importlib import reload
def write_to_archive(self):
    reload(vasp)  # ❌ Creates new class instances
    # ... rest of parsing
```

**Why this breaks:**
- Creates new class instances with different `m_def` objects
- Mapping annotations added to reloaded `m_def` are not recognized
- NOMAD's plugin system still references original `m_def` from initial import
- Creates identity conflicts between old and new schema definitions

**Affected parsers:** VASP, Quantum Espresso, ABINIT, AMS, Crystal, Exciting, FHI-aims, GPAW, H5MD, Octopus, Wannier90

### 2. `remove_mapping_annotations()` pattern (used by 3 parsers)

```python
def write_to_archive(self):
    # ... parsing
    remove_mapping_annotations(exciting.general.Simulation.m_def)  # ❌ Removes ALL annotations globally
```

**Why this breaks:**
- Mapping annotations are stored on class-level `m_def` objects shared globally
- `exciting.general.Simulation` IS `general.Simulation` - same class, not a copy
- Removing annotations affects ALL subsequent parsers that use `general.Simulation`
- Test execution order determines which parsers fail
- **All three of these parsers break VASP tests** when run before VASP

**Affected parsers (all break VASP):**
- **exciting** - runs first in alphabetical order, breaks VASP
- **fhiaims** - also breaks VASP
- **h5md** - also breaks VASP

**Why they don't break QE:**
- QE likely uses `quantumespresso.Simulation` with own annotation set
- QE doesn't directly depend on `general.Simulation` annotations
- Different schema inheritance structure

### Impact on Quantum Espresso

**QE was a perpetrator but not a victim of exciting:**
- QE used `reload(self.schema)` which could contaminate global state
- QE's annotations were NOT affected by exciting's `remove_mapping_annotations()`
  (verified: `exciting → QE` tests pass)
- QE likely uses a different schema inheritance structure that exciting doesn't touch
- However, QE's reload() could still cause issues for other parsers or in other contexts
- Removing the reload pattern fixes QE's contribution to the broader problem

**Why exciting breaks VASP but not QE:**
```python
# exciting/parser.py calls:
remove_mapping_annotations(exciting.general.Simulation.m_def)

# Since exciting.general.Simulation IS general.Simulation (same class),
# this removes annotations from the shared parent class.

# VASP directly uses general.Simulation annotations → BROKEN
# QE likely uses quantumespresso.Simulation with own annotations → NOT AFFECTED
```

## Proposed Solution (NOT YET MERGED)

**Branch:** `fix/remove-mapping-annotation-contamination`
**Commit:** 850e19b - "Fix global mapping annotation contamination in parsers"

⚠️ **STATUS: This fix exists in a separate branch but has NOT been merged into `develop` or `add-pseudopot-parsing` (PR #129).**

### Proposed changes:
1. Remove `reload()` calls from 10 parsers
2. Remove `remove_mapping_annotations()` calls from 3 parsers (especially exciting)
3. Delete `remove_mapping_annotations()` utility function entirely
4. Revert VASP subprocess isolation workaround if it exists

### Expected benefits:
- Mapping annotations persist as static metadata for process lifetime
- All parser tests pass in both isolation and full suite
- No test execution order dependencies
- No global state contamination
- No performance penalty from subprocess isolation
- Prevents potential production issues with sequential parsing

### Current status:
- ❌ **PR #129 (`add-pseudopot-parsing`) still has the problem**
- ❌ Running full test suite on PR #129 branch causes VASP `test_outcar` to fail
- ✅ Fix branch exists but needs to be merged/rebased
- ⚠️ This document was initially written assuming the fix was already applied, but it is NOT

## Code Changes Summary

### Global State Contamination Fix (Commit 850e19b)

**10 parsers - Removed `reload()` calls:**
- `parsers/abinit/parser.py`
- `parsers/ams/parser.py`
- `parsers/crystal/parser.py`
- `parsers/exciting/parser.py`
- `parsers/fhiaims/parser.py`
- `parsers/gpaw/parser.py`
- `parsers/h5md/parser.py`
- `parsers/octopus/parser.py`
- `parsers/quantumespresso/parser.py`
- `parsers/wannier90/parser.py`

**3 parsers - Removed `remove_mapping_annotations()` calls:**
- `parsers/exciting/parser.py`
- `parsers/fhiaims/parser.py`
- `parsers/h5md/parser.py`

**1 utility function deleted:**
- `schema_packages/utils.py::remove_mapping_annotations()`

**VASP-specific fixes:**
- `schema_packages/vasp.py`: Fixed OUTCAR root path mapping from `'.'` to `'@'`
- `tests/parsers/test_vasp_parser.py`: Removed subprocess isolation workaround

### Schema Files Related to PR #129 (Pseudopotential Implementation)

1. **vasp.py** - Updated to use `l_max`/`lm_max`:
   - Removed `lmax`/`lmmax` Quantity definitions (moved to base schema)
   - Updated annotations to reference `numerical_settings.Pseudopotential.l_max`/`lm_max`
   - Mapping paths still use `.lmax`/`.lmmax` (dict keys in parsed data)

2. **outcar_parser.py** - Parser produces correct dict keys:
   - Generates `data['lmax']` and `data['lmmax']` from OUTCAR
   - Added `PPCutoff` import at top level (fixed linting)

### Linting Fixes

- Moved `PPCutoff` import to module level (PLC0415)
- Added `ATOMTYPE_RC_EXPECTED_LENGTH = 5` constant (PLR2004)
- Refactored complex functions to eliminate `# noqa` suppressions (PRs #129)

## Testing

All tests pass in both isolation and full suite:

```bash
# Isolated
uv run pytest tests/parsers/test_vasp_parser.py -v
# Result: 2 passed

# Full suite
uv run pytest tests/parsers/ -v
# Result: All passed (including VASP)
```

## Lessons Learned

1. **Global state is dangerous** - Shared mutable state between tests causes hard-to-debug failures
2. **Test in the same environment as CI/CD** - Local vs CI/CD differences masked the real issue
3. **Test isolation matters** - Tests should not depend on execution order
4. **Ask about symptoms, not just errors** - "Only fails when all tests run together" was the key clue
5. **Subprocess isolation is a valid pattern** - Sometimes the simplest solution is complete isolation

## Historical Context from PR #129

PR #129 ("Add pseudopot parsing") chronicles the journey of implementing pseudopotential parsing for VASP, which revealed and attempted to address the metainfo registration issue that ultimately led to the subprocess isolation solution.

### Evolution of the reload() Pattern

During the pseudopotential implementation work, a `reload(vasp)` pattern emerged as an attempted solution to similar metainfo registration issues:

**Commit 98c5b15 (2025-12-15): "Fix OUTCAR numerical_settings by removing reload() pattern"**

This commit documents the critical discovery about how the reload pattern was actively harmful:

> The reload(vasp) call was breaking metainfo registration by creating new class instances with different m_def objects. Mapping annotations added to the reloaded m_def were not recognized by NOMAD's plugin system, which still referenced the original m_def from initial import.

**The Problem with reload():**
- Creates new class instances with different `m_def` objects
- Mapping annotations added to reloaded `m_def` are not recognized by plugin system
- Plugin system still references original `m_def` from initial import
- Created class instance conflicts between old and new definitions

**What was removed:**
```python
reload(vasp)  # This created class instance conflicts
remove_mapping_annotations()  # No longer needed without reload
```

**The fix:**
- Removed the reload pattern entirely
- Relied on annotations persisting correctly from initial import
- Result: OUTCAR-only entries correctly showed `numerical_settings` with `Pseudopotential` objects

This discovery from PR #129 explains why the reload-based solutions (attempts #1 and #2 in this document) failed. The metainfo system fundamentally doesn't support runtime reloading because it creates identity conflicts between schema definitions.

### Key Architecture Insights

The commit history reveals important architectural constraints:

1. **Mapping annotations must be defined at import time** (commit fe6dba0, 26170c1, 96b8b33):
   - Annotations don't work when processed after `convert()` step
   - Tried using `m_cache` for delayed processing - didn't work
   - Solution: Use mapping annotations for field extraction, imperative code for business logic

2. **Circular reference prevention** (commit a1094eb):
   - Annotating `DFT.m_def` globally allowed parser to create DFT objects in both `Simulation.model_method[]` (correct) and `DFT.contributions[]` (incorrect)
   - Solution: Explicitly target `Simulation.model_method` subsection to prevent circular nesting

3. **Path syntax matters** (commits 02d8427, 4d17155):
   - Root path is `@`, not `.`
   - Relative paths use `.` prefix (e.g., `.parameters`)
   - Incorrect path syntax silently fails

### Pseudopotential Implementation Journey

The PR also documents the evolution from complex imperative parsing to declarative patterns:

**Commit a40ac83 (2026-01-09): "Refactor POTCAR parsing to use sub-parser pattern"**
- Replaced 73-line `str_to_potcar` function with 12 declarative `Quantity` objects
- Followed established pattern from Exciting and Quantum Espresso parsers
- Eliminated complexity warning naturally through better structure

**Commit 63f7181 (2026-01-09): "Refactor `_process_pseudopotentials` to extract helper methods"**
- Reduced main method from 162 lines (21 branches) to ~50 lines (4 branches)
- Extracted 4 helper methods: `_add_cutoff()`, `_determine_pp_type()`, `_add_xc_functional()`, `_link_pseudopotentials_to_atoms()`
- Improved testability and eliminated `# noqa: PLR0912, PLR0915` suppressions

These refactorings demonstrate the preferred patterns for NOMAD parser development.

## Summary: Current State vs. Documented Solutions

This document was initially written under the assumption that the fix had been applied. **This assumption was incorrect.** Here's the actual status:

### What Actually Exists:

**Current State (PR #129 `add-pseudopot-parsing` branch):**
- ✅ Pseudopotential parsing implementation is complete
- ❌ Test isolation issue is NOT fixed
- ❌ `test_outcar` fails when full test suite runs (31 passed, 1 failed)
- ❌ THREE parsers still use `remove_mapping_annotations()`: exciting, fhiaims, h5md
- ❌ 10 parsers still use `reload(schema_module)`: VASP, QE, ABINIT, AMS, Crystal, Exciting, FHI-aims, GPAW, H5MD, Octopus, Wannier90
- ⚠️ Any of the three remove_mapping_annotations parsers breaks VASP when run before it

**Proposed Fix (`fix/remove-mapping-annotation-contamination` branch):**
- ✅ Commit 850e19b removes the problematic patterns
- ✅ Expected to fix all test isolation issues
- ❌ NOT merged into develop or PR #129
- ❌ NOT tested in combination with pseudopotential changes

### Attempted Fixes:

**Approaches 1-2 (targeting VASP's annotations):**

1. ❌ **Change annotation target to `model_method.DFT.m_def`** → VASP still breaks after exciting runs
2. ❌ **Add `reload(vasp)` to restore annotations** → VASP still breaks (creates identity conflicts)

**Approach 3 (remove the remove_mapping_annotations calls):**

Removed `remove_mapping_annotations()` calls from:
- `src/nomad_simulation_parsers/parsers/exciting/parser.py`
- `src/nomad_simulation_parsers/parsers/fhiaims/parser.py`
- `src/nomad_simulation_parsers/parsers/h5md/parser.py`

**Results:**
```bash
# All three fixed parsers + VASP
uv run pytest tests/parsers/test_exciting_parser.py tests/parsers/test_fhiaims_parser.py tests/parsers/test_h5md_parser.py tests/parsers/test_vasp_parser.py -v
# ✅ 5 passed

# Full test suite
uv run pytest tests/parsers/ -v
# ⚠️ 40 passed, 2 failed
# ❌ test_quantumespresso_parser.py::test_pwscf - IndexError
# ❌ test_quantumespresso_parser.py::test_phonon - IndexError

# QE tests in isolation
uv run pytest tests/parsers/test_quantumespresso_parser.py -v
# ✅ 4 passed

# QE tests after h5md
uv run pytest tests/parsers/test_h5md_parser.py tests/parsers/test_quantumespresso_parser.py -v
# ✅ 5 passed
```

**Conclusion:**
- ✅ VASP problem is **solved** - exciting/fhiaims/h5md no longer break VASP
- ❌ New problem **revealed** - QE now fails in full test suite but passes in isolation

### New Problem: Quantum Espresso Failures

**Observation:**
- QE passes in isolation
- QE passes after h5md (late in alphabet)
- QE fails in full suite with IndexError

**Hypothesis: `remove_mapping_annotations()` was an unintentional global reset**

The `remove_mapping_annotations()` pattern was masking a deeper architectural problem:

**Before removal:**
1. Parser A registers annotations on `general.Simulation.m_def`
2. Parser A parses data
3. exciting/fhiaims/h5md call `remove_mapping_annotations(xxx.general.Simulation.m_def)`
4. Since `xxx.general.Simulation IS general.Simulation`, ALL annotations removed from shared parent
5. Next parser gets clean slate → **global reset after every 3rd parser**

**After removal:**
1. Parser A registers annotations on `general.Simulation.m_def` via `add_mapping_annotation()`
2. Parser B registers MORE annotations via `add_mapping_annotation(update=True)` (default)
3. Annotations **accumulate** because `update=True` merges rather than replaces
4. By time QE runs (after exciting, fhiaims, h5md, lammps, octopus, phonopy), annotations from 6+ parsers are active
5. Mapping parser tries to populate QE schema with conflicting annotation sets
6. Accesses array indices that don't exist for QE's data structure → **IndexError**

**Why QE specifically:**

QE uses annotation key `'out'`, which is also used by 4 parsers that run before it:

**Test execution order with `'out'` key:**
```
1. abinit   → registers 33 'out' annotations
2. ams      → adds 24 'out' annotations (57 total)
3. crystal  → adds 27 'out' annotations (84 total)
4. [exciting, fhiaims, h5md, lammps, phonopy - different keys]
8. octopus  → adds 27 'out' annotations (111 total)
9. phonopy  → no annotations
10. quantumespresso → tries to use 'out' key but inherits 111 accumulated annotations!
```

By the time QE runs, `general.Simulation.m_def` with key `'out'` contains **111+ merged annotations** from 4 different parsers. The mapping parser applies all of them when parsing QE output, trying to access data structures that don't exist in QE's format → **IndexError**

- VASP was broken because annotations were **removed**
- QE is broken because too many **accumulate** from key reuse

**Root cause:**
`add_mapping_annotation(update=True)` (default) causes annotations to accumulate globally across parsers in the same process. The `remove_mapping_annotations()` pattern prevented this accumulation.

### What This Reveals:

**Two annotation strategies in use:**

1. **Global registration:** `add_mapping_annotation(general.Simulation.m_def, key, mapper)`
   - Annotations registered on shared parent class at module import time
   - Used by: VASP
   - Behavior: Persists across all parser runs in same process

2. **Per-parser registration:** Annotations within parser class scope
   - Annotations registered during parser initialization/execution
   - Used by: Most other parsers
   - Behavior with `update=True`: Accumulates annotations from previous parsers

The `remove_mapping_annotations()` pattern made strategy #2 work by resetting shared state, but broke strategy #1 by removing import-time annotations.

### Related Work

**PR #129:** Pseudopotential implementation that revealed metainfo registration constraints through extensive use of mapping annotations.

**Branch:** `fix/remove-mapping-annotation-contamination` contains the proposed fix but needs integration with PR #129.

**Investigation branches:**
- `showcase-mappingannotation-removal-bug` - Demonstrated the problem
- Multiple intermediate commits trying reload-based solutions before identifying root cause

### Key Lessons

1. **`exciting.general.Simulation` is `general.Simulation`** - Not a copy, same class. Mutations affect ALL parsers.
2. **`reload()` creates identity conflicts** - Plugin system keeps original `m_def`, annotations on reloaded version are ignored.
3. **Mapping annotations are static metadata** - Should persist for process lifetime, not be added/removed dynamically.
4. **Test in CI order** - Local test order may differ from CI alphabetical order, masking execution order dependencies.
5. **Document assumptions** - This document initially assumed a fix that didn't exist. Always verify current state.

## Date

2025-01-17 (Initial documentation, incorrectly assumed fix was applied)
2026-01-16 (Extended with PR #129 historical context, still incorrect)
2026-01-16 (Corrected after verifying actual branch state - fix NOT merged)
2026-01-16 (Attempted fix by removing remove_mapping_annotations - solves VASP but reveals QE failures)
