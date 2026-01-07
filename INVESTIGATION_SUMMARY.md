# VASP Parser Investigation Summary

**Date:** 2026-01-07
**Context:** Investigation of VASP test failures and mapping annotation system

---

## Table of Contents

1. [Initial Problem](#initial-problem)
2. [Root Cause Analysis](#root-cause-analysis)
3. [The Real Fix](#the-real-fix)
4. [The False Trail](#the-false-trail)
5. [Test Enhancements](#test-enhancements)
6. [Tags Ordering Issue](#tags-ordering-issue)
7. [Branches and PRs](#branches-and-prs)
8. [Recommendations](#recommendations)

---

## Initial Problem

### Symptoms

VASP parser test `test_outcar` exhibited inconsistent behavior:
- ✅ **Passed** when run in isolation: `pytest tests/parsers/test_vasp_parser.py`
- ❌ **Failed** when run in full suite: `pytest tests/parsers/`

**Error message:**
```
AssertionError: No model_method in simulation
assert []
 +  where [] = Simulation(steps, model_system, model_method, outputs).model_method
```

### Investigation Timeline

1. **VASP_TEST_ISOLATION_FIX.md** documented a subprocess isolation workaround
2. Tests were modified to run in isolated subprocesses to avoid contamination
3. This masked the real issue but didn't fix the root cause

---

## Root Cause Analysis

### The Mapping Parser System

NOMAD's mapping parser uses **global annotations** stored on metainfo class definitions:

```python
# Annotations stored here:
Section.m_def.m_annotations['mapping'][parser_key] = MapperAnnotation(...)

# Example:
general.Simulation.m_def.m_annotations['mapping']['outcar'] = Mapper(mapper='@')
```

**Critical insight:** `m_def` objects are **class-level**, shared globally across all parser instances in the same Python process.

### The Contamination Pattern

Three components caused global state contamination:

#### 1. `reload(schema_module)` Pattern

**What it did:**
```python
def write_to_archive(self):
    reload(exciting)  # Creates NEW class instances with NEW m_def objects
    # ... parsing ...
```

**Why it broke:**
- First import: Classes get `m_def` objects registered in NOMAD's plugin system
- `reload()`: Creates NEW class instances with NEW `m_def` objects
- But NOMAD's plugin system still references ORIGINAL `m_def` from first import
- Annotations on reloaded `m_def` are never found during parsing

#### 2. `remove_mapping_annotations()` Function

**What it did:**
```python
def remove_mapping_annotations(property: Section, max_depth: int = 5):
    """Remove mapping annotations from section and all subsections recursively."""
    property.m_annotations.pop('mapping', None)
    # ... recursive removal of all subsections ...
```

**Why it broke:**
- Removes annotations from **shared global** `m_def` objects
- When Exciting parser removes annotations, it removes them for **ALL** parsers
- VASP parser runs after → annotations missing → parsing fails

#### 3. Parsers Affected

**Used `reload()` + `remove_mapping_annotations()`:**
- Exciting parser
- FHI-aims parser
- H5MD parser

**Used `reload()` only:**
- QuantumESPRESSO parser
- AMS parser
- Abinit parser
- GPAW parser
- Octopus parser
- Wannier90 parser
- Crystal parser

### Why Tests Failed in Full Suite

1. Test suite runs parsers sequentially in same process
2. Exciting/FHI-aims/H5MD parser runs → removes annotations globally
3. VASP parser runs next → annotations missing → `model_method` empty → test fails
4. In isolation: VASP runs alone with fresh annotations → test passes

### Production Impact

**Affected scenarios:**
- Sequential parsing in Celery workers (different file types in same worker)
- Long-running workers (first contaminating parse breaks all subsequent parses)
- Development servers with hot-reload

**Not affected:**
- Multi-process deployments (each process has separate memory)
- Fresh worker per task (annotations reset)

---

## The Real Fix

**Branch:** `fix/remove-mapping-annotation-contamination`
**PR:** #137 → `develop`

### Changes Made

1. **Removed `reload()` from 10 parsers:**
   - quantumespresso/parser.py
   - ams/parser.py
   - abinit/parser.py
   - gpaw/parser.py
   - octopus/parser.py
   - wannier90/parser.py
   - crystal/parser.py
   - exciting/parser.py
   - fhiaims/parser.py
   - h5md/parser.py

2. **Removed `remove_mapping_annotations()` calls from 3 parsers:**
   - exciting/parser.py
   - fhiaims/parser.py
   - h5md/parser.py

3. **Deleted utility function:**
   - `schema_packages/utils.py` - removed `remove_mapping_annotations()` function entirely

4. **Reverted VASP subprocess workaround:**
   - `tests/parsers/test_vasp_parser.py` - removed subprocess isolation, restored simple tests

### Results

- All 42 parser tests pass in both isolation and full suite
- No subprocess workarounds needed
- Annotations correctly persist as static metadata for process lifetime

### Correct Architecture

**Annotations should:**
- ✅ Be loaded once at module import
- ✅ Persist for the process lifetime
- ✅ Be reused across all parses
- ❌ Never be reloaded
- ❌ Never be removed

**Rationale:** Annotations are **metadata** (how to parse), not **data** (what was parsed).

---

## The False Trail

### Branch: `add-pseudopot-parsing`

Between commits `c17c206` and `4d17155`, we made 8 commits attempting to "fix" the VASP parser:

| Commit | Description | What It Changed |
|--------|-------------|-----------------|
| 98c5b15 | Fix OUTCAR numerical_settings by removing reload() | Removed reload/remove_annotations + added test assertions |
| 2f907a0 | Add missing model_method annotation | Added explicit `model_method` mapping |
| db3b052 | Fix formatting error | Formatting fixes |
| 42b5d71 | Remove logging from debugging | Removed debug logs |
| 74e1a7e | Remove incorrect model_method annotation | Removed explicit mapping that interfered |
| 02d8427 | Fix OUTCAR DFT mapping path syntax | Changed `'parameters'` → `'.parameters'` |
| 852201d | Test failing error further | Debug changes |
| 4d17155 | Fix OUTCAR root path '.' → '@' | Changed root path |

### The Discovery

**Experiment at commit c17c206 (before all "fixes"):**

```bash
$ git checkout c17c206
$ uv run nomad parse --show-archive OUTCAR > test.outcar-old.json
$ git checkout 4d17155  # After all "fixes"
$ uv run nomad parse --show-archive OUTCAR > test.outcar-new.json
$ diff test.outcar-old.json test.outcar-new.json
# Result: Files are IDENTICAL (except timestamps and tag ordering)
```

**Conclusion:** The parser was already working correctly at c17c206. All 8 commits were chasing a ghost.

### Mapping Annotation Values Through Commits

**At c17c206 (working state):**
```python
add_mapping_annotation(general.Simulation.m_def, OUTCAR_KEY, '@')  # ✅ Correct
add_mapping_annotation(model_method.DFT.m_def, OUTCAR_KEY, 'parameters')  # Works
```

**The "fixes" alternated between:**
- Root path: `'@'` vs `'.'` (both times changed, it was already `'@'`)
- DFT path: `'parameters'` vs `'.parameters'` (both apparently work)
- Adding/removing explicit `model_method` mapping (not needed)

**Why tests appeared to fail:** Global annotation contamination from other parsers, not VASP's mapping configuration.

### Actual Issue

The test failures were due to:
1. **Other parsers** removing annotations globally
2. **Test execution order** determining which parser ran first
3. **VASP running after** a parser that removed annotations

Not due to VASP's own mapping configuration.

---

## Test Enhancements

### Commit 98c5b15 Changes

**Before (minimal tests):**
```python
def test_vasprun():
    parser = VASPParser()
    archive = EntryArchive()
    parser.parse('tests/data/vasp/AgAc_relax/vasprun.xml.relax', archive, LOGGER)

def test_outcar():
    parser = VASPParser()
    archive = EntryArchive()
    parser.parse('tests/data/vasp/AgAc_relax/OUTCAR', archive, LOGGER)
```

**After (enhanced tests):**
```python
def test_outcar():
    parser = VASPParser()
    archive = EntryArchive()
    parser.parse('tests/data/vasp/AgAc_relax/OUTCAR', archive, LOGGER)

    # Validation assertions
    simulation = archive.data
    assert simulation is not None, 'No simulation data in archive'
    assert simulation.model_method, 'No model_method in simulation'

    method = simulation.model_method[0]
    assert hasattr(method, 'numerical_settings'), (
        'ModelMethod missing numerical_settings attribute'
    )
    assert method.numerical_settings is not None, 'numerical_settings is None'
    assert len(method.numerical_settings) > 0, 'numerical_settings is empty'

    # Check for Pseudopotentials
    pseudopotentials = [
        ns for ns in method.numerical_settings
        if type(ns).__name__ == 'Pseudopotential'
    ]
    assert len(pseudopotentials) > 0, 'No Pseudopotential objects in numerical_settings'

    # Verify OUTCAR-specific fields
    for pp in pseudopotentials:
        assert pp.name is not None
        assert pp.n_valence_electrons is not None

    # ALSO ADDED: Extensive LOGGER.info() debugging output
```

### Analysis

**Valuable additions:**
- ✅ Assertions verifying simulation structure
- ✅ Validation of model_method population
- ✅ Verification of numerical_settings
- ✅ Checks for Pseudopotential objects
- ✅ Field-level validation

**Not valuable (per user preference):**
- ❌ Extensive `LOGGER.info()` debugging output
- ❌ Structure inspection logging
- ❌ Pseudopotential detail logging

**Status:** These enhancements are on the `add-pseudopot-parsing` branch but not in the clean fix (PR #137).

---

## Tags Ordering Issue

### Observation

When running the same parse twice, output JSON differs only in tag array ordering:

```diff
601,605c601,602
<       "Perdew, Burke & Ernzerhof",
<       "TotalForce",
<       "PBE",
<       "wavefunction",
<       "XC_GGA_C_PBE",
---
>       "XC_GGA_X_PBE",
>       "s2p4",
```

### Cause

**Non-deterministic ordering due to:**
1. Tags collected from dictionaries/sets with no guaranteed iteration order
2. Python's hash randomization between runs (`PYTHONHASHSEED` varies)
3. No explicit sorting applied to tags array before JSON serialization

### Impact

**Functionality:** None - all data present, just different order
**Testing:** Makes file comparison difficult, diffs are noisy
**Correctness:** Not a bug - tags are complete and correct

### Current State

- All tags present in both versions
- Content is identical
- Only presentation differs

### Potential Fix (if desired)

Add sorting when generating tags array to ensure deterministic output:
```python
archive.metadata.tags = sorted(set(tags))
```

Location TBD - need to investigate where tags are collected and serialized.

---

## Branches and PRs

### Branch: `fix/remove-mapping-annotation-contamination`

**Status:** ✅ Merged to `develop` via PR #137
**Commits:** 1
**Files changed:** 13 files, +20 -109 lines

**Purpose:** Fix global annotation contamination by removing broken patterns

**Results:**
- All 42 parser tests pass
- No subprocess isolation needed
- Clean, correct architecture

### Branch: `add-pseudopot-parsing`

**Status:** 🔄 Active, needs cleanup
**Commits:** 25 from `develop`
**Current HEAD:** `4d17155`

**Purpose:** Add pseudopotential parsing for VASP

**Core work (valuable):**
- Commits 1bf3455 → c17c206: Pseudopotential extraction implementation
- PPCutoff subsections, l_max/lm_max support
- Parsing from vasprun.xml and OUTCAR
- XC functional mapping, type detection

**False trail (can be dropped):**
- Commits 98c5b15 → 4d17155: Attempts to "fix" already-working parser
- 8 commits that made no functional difference
- Parser was correct at c17c206

**Current state at c17c206:**
- Parsing works correctly
- Output identical to final state (4d17155)
- Tests pass in isolation

**Conflicts with other branch:**
- Has the old `reload()` and `remove_mapping_annotations()` patterns
- Will conflict when rebasing onto `develop` (which now has PR #137)
- These files now have opposite patterns in different branches

---

## Recommendations

### For `add-pseudopot-parsing` Branch

**Option A: Reset and rebase (cleanest)**
```bash
git checkout add-pseudopot-parsing
git reset --hard c17c206  # Drop the 8 false-trail commits
git rebase develop        # Incorporate the real fix from PR #137
# Result: Clean branch with just pseudopotential work
```

**Option B: Interactive rebase (selective)**
```bash
git checkout add-pseudopot-parsing
git rebase -i c17c206
# Mark commits 98c5b15-4d17155 as 'drop'
# Keep test assertions from 98c5b15 if desired (sans LOGGER.info)
git rebase develop
```

**Option C: Cherry-pick test assertions only**
```bash
git checkout add-pseudopot-parsing
git reset --hard c17c206
# Manually add just the assertion logic from 98c5b15
# Skip all LOGGER.info() debugging
git rebase develop
```

### For Future Development

**Architecture rules:**
1. ✅ Annotations loaded once at import, persist forever
2. ❌ Never `reload()` schema modules
3. ❌ Never `remove_mapping_annotations()`
4. ✅ Trust the metainfo system's class-level design

**Testing rules:**
1. ✅ Validate structure with assertions
2. ✅ Check critical fields are populated
3. ❌ Don't add extensive debugging logs to tests
4. ✅ Keep tests simple and focused

**Debugging approach:**
1. Use `nomad parse --show-archive` to inspect output
2. Compare JSON diffs for validation
3. Keep debugging output out of committed tests

### For Tags Ordering

**If deterministic output desired:**
1. Locate where tags are collected (likely in archive metadata generation)
2. Add sorting: `tags = sorted(set(tags))`
3. Minimal change, ensures consistent output
4. Makes diffs clean and file comparison reliable

**If not critical:**
- Current behavior is functionally correct
- Just cosmetic ordering difference
- Can live with it for now

---

## Lessons Learned

1. **Global state is dangerous:** Shared mutable state between tests causes hard-to-debug failures
2. **Test in the same environment as CI/CD:** Local vs CI/CD differences masked the real issue
3. **Test isolation matters:** Tests should not depend on execution order
4. **Understand the architecture before fixing:** We spent 8 commits "fixing" something that wasn't broken
5. **Validate assumptions:** Running the parser directly proved it was already working
6. **Subprocess isolation is a smell:** It's a workaround, not a fix - investigate the root cause
7. **Simple diffs reveal truth:** Comparing actual output cuts through confusion
8. **Annotations are metadata:** They define behavior, not store results - they should persist

---

## Files Referenced

**Documentation:**
- `VASP_TEST_ISOLATION_FIX.md` - Original workaround documentation
- `INVESTIGATION_SUMMARY.md` - This file

**Code locations:**
- `src/nomad_simulation_parsers/schema_packages/utils.py` - Contained problematic functions
- `src/nomad_simulation_parsers/schema_packages/vasp.py` - Mapping annotations
- `src/nomad_simulation_parsers/parsers/vasp/parser.py` - VASP parser implementation
- `tests/parsers/test_vasp_parser.py` - Test file

**Parsers affected:**
- `parsers/exciting/parser.py`
- `parsers/fhiaims/parser.py`
- `parsers/h5md/parser.py`
- `parsers/quantumespresso/parser.py`
- `parsers/ams/parser.py`
- `parsers/abinit/parser.py`
- `parsers/gpaw/parser.py`
- `parsers/octopus/parser.py`
- `parsers/wannier90/parser.py`
- `parsers/crystal/parser.py`

---

## Next Steps

1. **Decide on branch cleanup strategy** for `add-pseudopot-parsing`
2. **Rebase onto `develop`** to incorporate PR #137 fixes
3. **Remove false-trail commits** (98c5b15 → 4d17155)
4. **Optionally add test assertions** (without debug logging)
5. **Verify pseudopotential parsing** still works after rebase
6. **Create PR** for pseudopotential feature

---

**Document Version:** 1.0
**Last Updated:** 2026-01-07
**Authors:** Investigation by Claude and Nathan Daelman
