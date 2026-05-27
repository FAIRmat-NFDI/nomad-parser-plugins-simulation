# Archive Examples

This directory contains parsed NOMAD archive examples demonstrating the data structure, particularly for SCF convergence data.

## Files

### `fhiaims_si_geomopt_archive.json` (15 KB)
- **Parser**: FHI-aims
- **Test file**: `tests/data/fhiaims/Si_geomopt/out.out`
- **Calculation type**: Geometry optimization
- **Content**:
  - 6 `model_system` sections (initial + 5 optimization steps)
  - 5 `outputs` sections (one per geometry step)
  - SCF convergence data in each `outputs[i].scf_steps`
  - Example: 11 SCF iterations for first geometry step

**Use case**: Reference for SCF convergence data structure in DFT calculations.

### `crystal_mgo_archive.json` (258 B)
- **Parser**: CRYSTAL
- **Test file**: `tests/data/crystal/MgO/output/MgO.out`
- **Calculation type**: Single point
- **Content**:
  - Minimal archive structure
  - No model_system or outputs sections

**Use case**: Example of minimal valid archive.

### `gromacs_water_archive.json` (15 MB)
- **Parser**: GROMACS
- **Test file**: `tests/data/gromacs/water/reference_s.log`
- **Calculation type**: Molecular dynamics
- **Content**:
  - 51 `model_system` sections (MD trajectory frames)
  - 11 `outputs` sections (sampled MD steps)
  - Force field `model_method`

**Use case**: Reference for MD trajectory data structure.

## Generation Details

**Generated**: 2026-05-22

**Branches**:
- `nomad-FAIR`: `mapping-parser-multi-defs` (commit 7ef8221cf - MR 3040 with bug fixes)
- `nomad-parser-plugins-simulation`: `develop`

**Purpose**: These archives serve as reference examples for:
1. Understanding the NOMAD archive data structure
2. Verifying SCF convergence data layout (`outputs[i].scf_steps`)
3. Testing schema changes and normalizers
4. Documentation and training materials

## SCF Data Structure

The key pattern for SCF convergence data:

```python
archive.data.outputs[i].scf_steps  # Single SCFSteps object (not a list)
  ├─ delta_energies_total   # Array: (n_iterations,) in joule
  ├─ delta_density_rms      # Array: (n_iterations,) in coulomb
  ├─ durations              # Array: (n_iterations,) in second
  └─ code_specific_quantities  # Dict with parser-specific data
```

Each array element represents one SCF iteration. All arrays have the same length = number of SCF iterations for that calculation.

## Regeneration

To regenerate these examples:

```bash
# From nomad-parser-plugins-simulation root
uv run python /tmp/generate_archive_examples.py

# Or manually:
from nomad.datamodel import EntryArchive
from nomad_simulation_parsers.parsers.fhiaims.parser import FHIAimsParser

parser = FHIAimsParser()
archive = EntryArchive()
parser.parse('tests/data/fhiaims/Si_geomopt/out.out', archive)

import json
with open('tests/archive_examples/fhiaims_si_geomopt_archive.json', 'w') as f:
    json.dump(json.loads(archive.m_to_json()), f, indent=2)
```

## Notes

- Archives are serialized using `archive.m_to_json()` for full metainfo compliance
- File sizes vary significantly based on calculation type (geometry opt vs MD trajectory)
- GROMACS archive is large due to 51 trajectory frames with full atomic coordinates
- CRYSTAL archive is minimal because the test file doesn't produce full output sections
