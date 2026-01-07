# PPCutoff Implementation Summary

## Changes Made

### 1. Schema Changes (`nomad-schema-plugins-simulations`)

#### `numerical_settings.py`
- **Added `PPCutoff` class** (lines 1035-1093):
  - `cutoff_kind`: MEnum specifying physical expansion type (wavefunction, charge_density, augmentation, response, unavailable)
  - `cutoff_role`: MEnum specifying precision context (recommended, recommended_min, recommended_max, fast, balanced, stringent)
  - `value`: Energy value in joules
  - Follows SSSP nomenclature (Beal et al., arXiv:2504.03962, 2025)

- **Modified `Pseudopotential` class** (line 1200):
  - Replaced single `cutoff` Quantity with `cutoffs` SubSection
  - `cutoffs = SubSection(sub_section=PPCutoff.m_def, repeats=True)`

### 2. Parser Changes (`nomad-parser-plugins-simulation`)

#### `vasp.py` schema
- **Removed from `Pseudopotential` class**:
  - `enmax` Quantity
  - `enmin` Quantity  
  - `cutoff_target` Quantity
- **Kept VASP-specific fields**:
  - `sha256`: POTCAR file hash
  - `l_max`: Maximum angular momentum
  - `lm_max`: Total lm-projection operators

#### `outcar_parser.py`
- **Added import**: `from nomad.units import ureg`
- **Modified `_process_pseudopotentials()` method**:
  - Creates `PPCutoff` for ENMAX:
    - `cutoff_kind='wavefunction'`
    - `cutoff_role='recommended'`
    - Converts eV to joules
  - Creates `PPCutoff` for ENMIN:
    - `cutoff_kind='wavefunction'`
    - `cutoff_role='recommended_min'`
    - Converts eV to joules

## Design Rationale

### Why PPCutoff Subsections?

Pseudopotential files often contain multiple cutoff recommendations:
- VASP POTCAR: ENMAX (recommended) and ENMIN (minimum)
- Quantum ESPRESSO UPF: ecutwfc and ecutrho (4× for charge density)
- SSSP libraries: efficiency and precision tiers

The `PPCutoff` subsection structure captures:
1. **Physical expansion type** (`cutoff_kind`): What this cutoff controls
2. **Precision context** (`cutoff_role`): Recommendation level or tier
3. **Value**: The actual cutoff energy

### VASP Mapping

| POTCAR Field | cutoff_kind | cutoff_role | Notes |
|--------------|-------------|-------------|-------|
| ENMAX        | wavefunction | recommended | Standard precision |
| ENMIN        | wavefunction | recommended_min | ~2/3 of ENMAX |

## Testing Status

- ✅ Code compiles successfully (Python syntax check)
- ✅ Ruff formatting and linting passes
- ⚠️ Runtime testing blocked by environment issue (NumPy 2.x / elasticsearch 7.9.1 incompatibility)
  - This is unrelated to our changes
  - Issue: elasticsearch 7.9.1 uses deprecated `np.float_` removed in NumPy 2.0

## Next Steps

1. Test with actual VASP OUTCAR/POTCAR files once environment is fixed
2. Verify `pp.cutoffs` is properly populated with two PPCutoff instances
3. Extend to other parsers (Quantum ESPRESSO, CASTEP, etc.)
