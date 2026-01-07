# VASP Test Isolation Fix

## Problem Summary

VASP parser tests (`test_vasprun` and `test_outcar`) were failing in CI/CD and when running the full test suite, but passing when run in isolation. The error was:

```
AssertionError: No model_method in simulation
assert []
 +  where [] = Simulation(steps, model_system, model_method, outputs).model_method
```

## Root Cause

**Test isolation issue caused by global state contamination.**

The NOMAD mapping parser system uses global mapping annotations that are registered when schema modules are imported. Other parsers in the test suite were removing these mapping annotations during their test teardown, which caused VASP tests to fail when they ran afterward.

### Why it only failed in CI/CD / full test suite:
- When running VASP tests alone: annotations are registered on import and remain intact
- When running full test suite: other parsers remove annotations, breaking VASP tests
- CI/CD always runs the full test suite, so it always failed there

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

### 3. Subprocess isolation (FINAL SOLUTION ✓)

Run VASP tests in isolated subprocesses where annotations can't be contaminated:

```python
@pytest.mark.parametrize('mainfile', [...])
def test_vasp_isolated(mainfile):
    script = f'''
    # Full test logic here
    '''
    result = subprocess.run([sys.executable, '-c', script], ...)
```

## Final Solution

**File:** `tests/parsers/test_vasp_parser.py`

- Combined `test_vasprun` and `test_outcar` into parametrized `test_vasp_isolated`
- Each test runs in fresh subprocess via `subprocess.run()`
- Subprocess has clean Python environment with fresh annotation state
- Tests are now immune to global state contamination

### Benefits:
- ✅ Tests pass reliably in full suite
- ✅ Tests pass in CI/CD
- ✅ No dependencies on test execution order
- ✅ Clean isolation without complex fixtures

### Tradeoffs:
- ⚠️ Slightly slower (subprocess overhead ~15s total vs ~4s)
- ⚠️ Workaround rather than fixing root cause

## Proper Long-term Fix

The real issue is that **other parsers should not remove VASP's mapping annotations**. The proper fix would be to:

1. Identify which parsers are removing annotations globally
2. Fix those parsers to only remove their own annotations
3. Or prevent annotation removal entirely (annotations should persist)

This would require investigating all parser teardown logic and the mapping parser framework.

## Code Changes Summary

### Schema Files (nomad-parser-plugins-simulation)

1. **vasp.py** - Updated to use `l_max`/`lm_max`:
   - Removed `lmax`/`lmmax` Quantity definitions (moved to base schema)
   - Updated annotations to reference `numerical_settings.Pseudopotential.l_max`/`lm_max`
   - Mapping paths still use `.lmax`/`.lmmax` (dict keys in parsed data)

2. **outcar_parser.py** - Parser produces correct dict keys:
   - Generates `data['lmax']` and `data['lmmax']` from OUTCAR
   - Added `PPCutoff` import at top level (fixed linting)

3. **test_vasp_parser.py** - Subprocess isolation:
   - Replaced two separate tests with one parametrized test
   - Runs each test file in isolated subprocess
   - Includes all assertions for OUTCAR validation

### Linting Fixes

- Moved `PPCutoff` import to module level (PLC0415)
- Added `ATOMTYPE_RC_EXPECTED_LENGTH = 5` constant (PLR2004)
- Added `# noqa` suppressions for complexity warnings (PLR0915, PLR0912)

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

## Date

2025-01-17
