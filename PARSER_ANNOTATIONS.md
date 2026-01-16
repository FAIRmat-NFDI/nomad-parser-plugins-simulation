# Parser Annotation Summary

This document lists all mapping annotations defined by each parser schema.

## Annotation Keys

- Annotation keys are the identifiers used to select which parser schema to use (e.g., `XML_KEY='xml'`, `OUTCAR_KEY='outcar'`)
- Each parser uses 1-5 different keys for different file formats

### **CRITICAL: Key Collisions**

**Multiple parsers use the same key values:**

| Key Value | Parsers Using It | Count |
|-----------|-----------------|-------|
| `'out'` | abinit, ams, crystal, octopus | **4** |
| `'info'` | exciting, octopus | **2** |
| `'dos'` | abinit, wannier90 | **2** |

**Impact:**
- When multiple parsers use the same key (e.g., `'out'`), they register annotations with that key on the **same shared classes**
- Since `add_mapping_annotation(update=True)` is the default, annotations **accumulate** rather than replace
- This means `general.Simulation.m_def` with key `'out'` contains merged annotations from abinit, ams, crystal, AND octopus
- The mapping parser tries to apply ALL accumulated annotations when parsing, even if they're from different parsers

**Why Quantum Espresso Fails:**

QE also uses key `'out'`! Test execution order:
```
1. abinit   - registers 'out' annotations
2. ams      - adds 'out' annotations (merged with abinit)
3. crystal  - adds 'out' annotations (merged with abinit + ams)
4. exciting - different keys
5. fhiaims  - different keys
6. h5md     - different keys
7. lammps   - different keys
8. octopus  - adds 'out' annotations (merged with abinit + ams + crystal)
9. phonopy  - no annotations
10. quantumespresso - tries to use 'out' key but gets accumulated annotations from 4 previous parsers!
```

By the time QE runs, the `'out'` key contains merged annotations from:
- `abinit.py`: 33 annotations
- `ams.py`: 24 annotations
- `crystal.py`: 27 annotations
- `octopus.py`: 27 annotations
- **Total: ~111 accumulated annotations** before QE adds its own

The mapping parser tries to apply all 111+ annotations when parsing QE output, accessing data structures that don't match QE's format → **IndexError**

## Counts by Parser

| Parser | Total Annotations | Keys Used |
|--------|-------------------|-----------|
| h5md | 114 | H5_KEY |
| vasp | 69 | XML_KEY, XML2_KEY, OUTCAR_KEY |
| wannier90 | 59 | OUT_KEY, BAND_KEY |
| fhiaims | 53 | OUT_KEY, AIMS_KEY |
| exciting | 45 | INFO_KEY, INPUT_XML_KEY, EIGVAL_KEY, BANDSTRUCTURE_XML_KEY, DOS_XML_KEY |
| abinit | 33 | OUT_KEY, DOS_KEY |
| octopus | 27 | STDOUT_KEY, STATIC_KEY |
| crystal | 27 | OUT_KEY |
| gpaw | 25 | GPAW_KEY |
| ams | 24 | LOG_KEY |
| phonopy | 0 | (no annotations) |
| quantumespresso | 0 | (annotations in submodules) |

## Annotation Targets by Parser

### Key Observation

Parsers annotate at **different levels**:

1. **Module-level annotations on `general.Simulation.m_def`**
   - VASP: Lines 25-27 (before class definition)
   - Exciting: Lines 23-27 (before class definition)
   - FHI-aims: Lines (before class definition)
   - All other parsers: Line after class definition

2. **Class-level annotations within `Simulation` class**
   - All parsers define these inside their `class Simulation(general.Simulation):` block

3. **Annotations on other schema classes**
   - `model_method.DFT.m_def`
   - `model_system.ModelSystem`
   - `outputs.Outputs`
   - etc.

The key difference is **when** module-level `general.Simulation.m_def` annotations are registered:
- **Before class definition** → Registered at module import time, persist globally
- **After class definition** → Still registered at module import time, but after subclass is defined

## Impact on Test Isolation

When parsers run in sequence (test suite):

1. **Module import** → All `add_mapping_annotation()` calls execute, annotations registered on shared classes
2. **Parser A runs** → Uses Parser A's annotations
3. **Parser B imports (if not already)** → Annotations from B **added to** existing annotations from A (because `update=True` default)
4. **Parser B runs** → May see annotations from both A and B on shared classes

The `remove_mapping_annotations()` pattern was resetting this state, but broke parsers that relied on import-time annotations persisting.

## Detailed Listings

### VASP (69 annotations, 3 keys)

**Module-level (before class):**
- `general.Simulation.m_def` + `XML_KEY` → `'modeling'`
- `general.Simulation.m_def` + `XML2_KEY` → `'modeling'`
- `general.Simulation.m_def` + `OUTCAR_KEY` → `'@'`

**Class Simulation:**
- `general.Simulation.program` + XML_KEY, OUTCAR_KEY
- `general.Simulation.model_method` + XML_KEY, OUTCAR_KEY
- `general.Simulation.model_system` + XML_KEY, OUTCAR_KEY
- `general.Simulation.outputs` + XML_KEY, XML2_KEY, OUTCAR_KEY

**Other classes:**
- Program, ModelSystem, AtomsState, BasisSetContainer, BasisSet
- DFT, XCFunctional, Smearing
- KSpace, KSpaceMesh, FrequencyMesh
- Pseudopotential (14 annotations for cutoffs, functionals, etc.)
- Outputs, TotalEnergy, TotalForce, various output properties

### EXCITING (45 annotations, 5 keys)

**Module-level (before class):**
- `general.Simulation.m_def` + INFO_KEY → `'@'`
- `general.Simulation.m_def` + INPUT_XML_KEY → `'@'`
- `general.Simulation.m_def` + EIGVAL_KEY → `'@'`
- `general.Simulation.m_def` + BANDSTRUCTURE_XML_KEY → `'@'`
- `general.Simulation.m_def` + DOS_XML_KEY → `'@'`

**Class Simulation:**
- `general.Simulation.program` + INFO_KEY
- Various model_method, model_system, outputs annotations

### FHI-AIMS (53 annotations, 2 keys)

**Module-level:**
- `general.Simulation.m_def` + OUT_KEY, AIMS_KEY

**Class Simulation:**
- Extensive annotations for program, model_method, model_system, outputs

### H5MD (114 annotations, 1 key)

**Most annotations of any parser** - molecular dynamics simulations

**Module-level:**
- `general.Simulation.m_def` + H5_KEY

**Extensive subsystems:**
- ModelSystem with hierarchical particle groups
- Workflow with MD method, thermostats, barostats
- Custom outputs for MD observables

### ABINIT (33 annotations, 2 keys)

**Module-level (after class):**
- `general.Simulation.m_def` + OUT_KEY, DOS_KEY

**Class Simulation:**
- Standard DFT workflow annotations

### Other Parsers

- **AMS** (24): Standard DFT with LOG_KEY
- **Crystal** (27): DFT calculations with OUT_KEY
- **GPAW** (25): DFT with GPAW_KEY
- **Octopus** (27): TDDFT with STDOUT_KEY, STATIC_KEY
- **Wannier90** (59): Wannier function analysis with OUT_KEY, BAND_KEY

## Current State (add-pseudopot-parsing branch)

**Parsers still using `reload()`:** 10 parsers
- abinit, ams, crystal, exciting, fhiaims, gpaw, h5md, octopus, quantumespresso, wannier90

**Parsers with `reload()` removed:**
- vasp (removed in PR #129)

**Parsers using `remove_mapping_annotations()`:** 0 parsers
- exciting, fhiaims, h5md (removed in recent commits to fix VASP)

## Summary

1. **VASP, Exciting, FHI-aims**: Register `general.Simulation.m_def` annotations at module import before class definition
2. **Other parsers**: Register annotations after class definition or only within class scope
3. **All parsers**: Use `add_mapping_annotation(update=True)` which accumulates annotations
4. **Problem**: Without `remove_mapping_annotations()`, annotations accumulate globally across test runs
5. **Result**: VASP now works (no longer broken by annotation removal), but QE fails (affected by annotation accumulation)
