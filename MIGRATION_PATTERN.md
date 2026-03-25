# Parser Migration Pattern: WorkflowConvergenceTarget

**Last Updated**: 2026-01-29
**Reference Implementation**: Exciting parser (`schema_update_convergence_targets` branch)
**Related PR**: #260 in nomad-simulations

---

## ⚠️ CRITICAL WARNING: Architectural Decision Pending

**STATUS (2026-01-29)**: **DO NOT USE THIS MIGRATION GUIDE YET**

The team is **actively reconsidering** the architectural approach documented in this guide. The current repeatable subsection pattern (described below) may be replaced with an alternative explicit fields approach.

### Why This Matters

During team review, the implementer acknowledged that the current approach has significant complexity:
- JMESPath integration is *"hacky"* and *"ugly"* (implementer's words)
- Unit handling requires workarounds
- Normalization logic is complex
- Multiple layers of indirection in parser implementation
- Integration testing not yet completed

The team is evaluating whether **simplicity should be prioritized over elegance**.

### Impact on Migration Work

**If you are a parser developer**:
- ⚠️ **STOP** any migration work using this guide until the architectural decision is made
- Following this guide now risks wasting 4-16 hours of work per parser
- The alternative approach (if chosen) will have a **simpler migration pattern**

**If the explicit fields approach is adopted**:
- Direct field mappings (like the old API, but improved)
- No JMESPath needed
- No dict-based method returns
- Simpler, more straightforward implementation
- This migration guide will be completely rewritten

**If the repeatable subsection approach continues**:
- This guide remains valid
- Exciting parser proves the pattern works
- Continue following the steps below

### Recommendation

**Wait for architectural decision** before proceeding with any parser migration work.

**See**: `packages/nomad-schema-plugins-simulations/ARCHITECTURAL_DECISION.md` for detailed analysis of both approaches.

**Timeline**: Decision expected within 1-2 team meetings.

---

## Overview

⚠️ **REMINDER**: The pattern described below assumes the repeatable subsection approach continues. Do not proceed until architectural decision is confirmed.

This guide documents the migration pattern for updating parsers to use the new `WorkflowConvergenceTarget` system introduced in nomad-simulations PR #260. This replaces the old convergence tolerance fields in `GeometryOptimizationMethod`.

### What Changed

**OLD API (Removed)**:
```python
GeometryOptimizationMethod:
  - convergence_tolerance_energy_difference
  - convergence_tolerance_force_maximum
  - convergence_tolerance_stress_maximum
  - convergence_tolerance_displacement_maximum
```

**NEW API (Use This)**:
```python
GeometryOptimizationMethod:
  - convergence: [WorkflowConvergenceTarget]  # repeatable subsection
  - single_point_convergence: [WorkflowConvergenceTarget]  # for SCF within geometry opt

WorkflowConvergenceTarget:
  - convergence_parameter_name: MEnum('energy', 'force', 'potential', 'charge', 'density')
  - threshold_type: MEnum('absolute', 'relative', 'maximum', 'rms', 'residuum')
  - convergence_threshold: float
  - threshold_unit: str
```

---

## Migration Pattern (Based on Exciting Parser)

### Step 1: Create Convergence Extraction Methods

Add methods to your parser class that return convergence data as dictionaries:

#### Example 1: Geometry Optimization Convergence

```python
def get_geometry_convergence(self, source: dict[str, Any]) -> list[dict[str, Any]]:
    """
    Extract geometry optimization convergence criteria.

    Returns:
        List of dicts, one per convergence criterion. Each dict must have:
        - convergence_parameter_name: str (e.g., 'force', 'energy')
        - threshold_type: str (e.g., 'maximum', 'absolute')
        - convergence_threshold: float (numerical value)
        - convergence_threshold_unit: str (e.g., 'newton', 'joule')
    """
    structure_optimization = source.get('structure_optimization')
    if structure_optimization is None:
        return []

    convergence = []

    # Force convergence (example)
    force_threshold = structure_optimization.get('force_target')
    if force_threshold is not None:
        convergence.append({
            'convergence_parameter_name': 'force',
            'threshold_type': 'maximum',
            'convergence_threshold': force_threshold.to('newton'),
            'convergence_threshold_unit': 'newton'
        })

    # Energy convergence (example)
    energy_threshold = structure_optimization.get('energy_target')
    if energy_threshold is not None:
        convergence.append({
            'convergence_parameter_name': 'energy',
            'threshold_type': 'absolute',
            'convergence_threshold': energy_threshold.to('joule'),
            'convergence_threshold_unit': 'joule'
        })

    return convergence
```

#### Example 2: Single Point (SCF) Convergence

```python
# Define mapping from code-specific names to WorkflowConvergenceTarget
convergence_threshold_mapping = {
    'x_exciting_effective_potential_convergence': {
        'name': 'potential',
        'type': 'rms',
        'unit': 'joule'
    },
    'x_exciting_energy_convergence': {
        'name': 'energy',
        'type': 'absolute',
        'unit': 'joule'
    },
    'x_exciting_charge_convergence': {
        'name': 'charge',
        'type': 'absolute',
        'unit': 'coulomb'
    },
    'x_exciting_IBS_force_convergence': {
        'name': 'force',
        'type': 'absolute',
        'unit': 'newton'
    }
}

def get_single_point_convergence(self, source: dict[str, Any]) -> list[dict[str, Any]]:
    """
    Extract SCF convergence criteria.

    Returns:
        List of dicts with convergence target information.
    """
    last_iteration = source.get('groundstate', {}).get('scf_iteration', [])[-1]
    convergence_targets = []

    for key, info in convergence_threshold_mapping.items():
        quantity = last_iteration.get(key)
        if quantity is None:
            continue

        convergence_targets.append({
            'convergence_parameter_name': info['name'],
            'threshold_type': info['type'],
            'convergence_threshold': quantity[1].to(info['unit']),
            'convergence_threshold_unit': info['unit']
        })

    return convergence_targets
```

### Step 2: Add Mapping Annotations in Schema File

In your parser's schema file (e.g., `schema_packages/yourparser.py`), add mapping annotations:

#### For Convergence Methods

```python
from nomad_simulations.schema_packages import workflow

# Map the convergence subsections to your parser methods
add_mapping_annotations(
    # Geometry optimization convergence
    (workflow.geometry_optimization.GeometryOptimizationMethod.convergence,
     INFO_KEY, ('get_geometry_convergence', ['.@'])),

    # Single point convergence within geometry optimization
    (workflow.geometry_optimization.GeometryOptimizationMethod.single_point_convergence,
     INFO_KEY, ('get_single_point_convergence', ['.@'])),

    # Single point convergence for standalone SinglePoint workflows
    (workflow.single_point.SinglePointMethod.convergence,
     INFO_KEY, ('get_single_point_convergence', ['.@']))
)
```

#### For WorkflowConvergenceTarget Fields

```python
# Direct mapping of WorkflowConvergenceTarget properties
# These connect the dict keys returned by your methods to the schema fields
add_mapping_annotations(
    # For INFO_KEY context (typical single point)
    (workflow.general.WorkflowConvergenceTarget.convergence_parameter_name,
     INFO_KEY, '.convergence_parameter_name'),
    (workflow.general.WorkflowConvergenceTarget.convergence_threshold,
     INFO_KEY, '.convergence_threshold'),
    (workflow.general.WorkflowConvergenceTarget.threshold_type,
     INFO_KEY, '.threshold_type'),
    (workflow.general.WorkflowConvergenceTarget.threshold_unit,
     INFO_KEY, '.convergence_threshold_unit'),

    # For GEO_OPT_KEY context (geometry optimization)
    (workflow.general.WorkflowConvergenceTarget.convergence_parameter_name,
     GEO_OPT_KEY, '.convergence_parameter_name'),
    (workflow.general.WorkflowConvergenceTarget.convergence_threshold,
     GEO_OPT_KEY, '.convergence_threshold'),
    (workflow.general.WorkflowConvergenceTarget.threshold_type,
     GEO_OPT_KEY, '.threshold_type'),
    (workflow.general.WorkflowConvergenceTarget.threshold_unit,
     GEO_OPT_KEY, '.convergence_threshold_unit')
)
```

### Step 3: Migrate SCFSteps (if applicable)

If your parser populates SCF iteration data, migrate from the old `SCFOutputs` to new `SCFSteps`:

```python
def get_scf_steps(self, source: dict[str, Any]) -> dict[str, Any]:
    """
    Extract SCF iteration data.

    Returns:
        Dict with arrays for each SCF quantity:
        - energies_total: array of total energies at each iteration
        - durations: array of time spent in each iteration
        - delta_energies_total: array of energy changes (optional)
        - delta_potential_rms: array of potential RMS changes (optional)
        - delta_density_rms: array of density RMS changes (optional)
        - delta_force_abs: array of absolute force changes (optional)
    """
    scf_steps = source.get('groundstate', {}).get('scf_iteration', [])

    energies = []
    wall_times = []
    delta_energies = []
    delta_potential = []
    delta_charge = []
    delta_force = []

    def safe_append(source, value_name, unit_conversion, out):
        """Append values only if they exist for the current step"""
        value = source.get(value_name)
        if value is not None:
            out.append(value.to(unit_conversion)[0])

    for idx, step in enumerate(scf_steps):
        energies.append(step.get('energy_total').to('joule'))
        wall_times.append(step.get('time_physical').to('seconds').magnitude)

        # Optional convergence deltas
        safe_append(step, 'x_exciting_energy_convergence', 'joule', delta_energies)
        safe_append(step, 'x_exciting_effective_potential_convergence',
                    'joule', delta_potential)
        safe_append(step, 'x_exciting_charge_convergence', 'coulomb', delta_charge)
        safe_append(step, 'x_exciting_IBS_force_convergence', 'newton', delta_force)

    # Compute durations from cumulative wall times
    durations = []
    for idx, time in enumerate(wall_times):
        duration = time if idx == 0 else time - wall_times[idx-1]
        durations.append(duration)

    # Build output dict (only include non-empty arrays)
    out = {
        'energies_total': energies,
        'durations': durations
    }

    for name, values in zip(
        ['delta_energies_total', 'delta_potential_rms',
         'delta_density_rms', 'delta_force_abs'],
        [delta_energies, delta_potential, delta_charge, delta_force]
    ):
        if len(values) > 0:
            out[name] = values

    return out
```

Add mapping annotations:

```python
add_mapping_annotations(
    (outputs.Outputs.scf_steps, INFO_KEY, ('get_scf_steps', ['@'])),
    (outputs.SCFSteps.durations, INFO_KEY, '.durations'),
    (outputs.SCFSteps.energies_total, INFO_KEY, '.energies_total'),
    (outputs.SCFSteps.delta_energies_total, INFO_KEY, '.delta_energies_total'),
    (outputs.SCFSteps.delta_potential_rms, INFO_KEY, '.delta_potential_rms'),
    (outputs.SCFSteps.delta_density_rms, INFO_KEY, '.delta_density_rms'),
    (outputs.SCFSteps.delta_force_abs, INFO_KEY, '.delta_force_abs')
)
```

---

## Common Patterns by Parser Type

### For Parsers That Only Do Geometry Optimization

```python
# In parser class:
def get_geometry_convergence(self, source):
    # Extract force/energy/displacement thresholds
    return [...]

# In schema file:
add_mapping_annotation(
    workflow.geometry_optimization.GeometryOptimizationMethod.convergence,
    YOUR_KEY, ('get_geometry_convergence', ['.@'])
)
```

### For Parsers That Do Both Geometry Opt + SCF

```python
# In parser class:
def get_geometry_convergence(self, source):
    # Geometry-level convergence (force, displacement)
    return [...]

def get_single_point_convergence(self, source):
    # SCF-level convergence (energy, potential, density)
    return [...]

# In schema file:
add_mapping_annotations(
    (workflow.geometry_optimization.GeometryOptimizationMethod.convergence,
     GEO_OPT_KEY, ('get_geometry_convergence', ['.@'])),
    (workflow.geometry_optimization.GeometryOptimizationMethod.single_point_convergence,
     GEO_OPT_KEY, ('get_single_point_convergence', ['.@']))
)
```

### For Parsers That Only Do Single Point (SCF)

```python
# In parser class:
def get_single_point_convergence(self, source):
    # SCF convergence criteria
    return [...]

# In schema file:
add_mapping_annotation(
    workflow.single_point.SinglePointMethod.convergence,
    INFO_KEY, ('get_single_point_convergence', ['.@'])
)
```

---

## Field Reference

### convergence_parameter_name

**Type**: `MEnum`
**Values**: `'energy'`, `'force'`, `'potential'`, `'charge'`, `'density'`
**Description**: What physical quantity is being converged

**Examples**:
- `'energy'` - Total energy convergence (SCF or geometry opt)
- `'force'` - Force convergence (geometry opt)
- `'potential'` - Effective potential convergence (SCF)
- `'charge'` / `'density'` - Electron density convergence (SCF)

### threshold_type

**Type**: `MEnum`
**Values**: `'absolute'`, `'relative'`, `'maximum'`, `'rms'`, `'residuum'`
**Description**: How the convergence check is performed

**Examples**:
- `'absolute'` - Absolute difference between iterations: `|E_n - E_{n-1}|`
- `'relative'` - Relative difference: `|E_n - E_{n-1}| / |E_n|`
- `'maximum'` - Maximum value across components: `max|F_i,n - F_i,{n-1}|`
- `'rms'` - Root mean square: `sqrt(sum|F_i,n - F_i,{n-1}|^2 / N)`
- `'residuum'` - Difference from initial estimate

### convergence_threshold

**Type**: `float`
**Description**: Numerical threshold value (without units)

**Examples**:
- `1e-6` (for energy in joules)
- `0.001` (for forces in eV/Angstrom converted to newtons)

### threshold_unit

**Type**: `str`
**Description**: Unit of the threshold (pint notation)

**Examples**:
- `'joule'` (for energy)
- `'newton'` (for force)
- `'coulomb'` (for charge/density)
- `'pascal'` (for stress)
- `'meter'` (for displacement)

---

## Complete Example: Abinit Parser Migration

### Before (Broken):

```python
# schema_packages/abinit.py
class GeometryOptimizationMethod(
    workflow.geometry_optimization.GeometryOptimizationMethod
):
    add_mapping_annotation(
        workflow.geometry_optimization.GeometryOptimizationMethod.convergence_tolerance_energy_difference,  # ❌ REMOVED
        OUT_KEY,
        ('get_input_var', [], dict(name='tolmxde', n_dataset=1, default=0.0)),
        unit='hartree',
    )
    add_mapping_annotation(
        workflow.geometry_optimization.GeometryOptimizationMethod.convergence_tolerance_force_maximum,  # ❌ REMOVED
        OUT_KEY,
        ('get_input_var', [], dict(name='tolmxf', n_dataset=1, default=0.0)),
        unit='hartree/bohr',
    )
```

### After (Fixed):

```python
# parsers/abinit/parser.py
class OutParser(TextParser):
    def get_geometry_convergence(self, source: dict[str, Any]) -> list[dict[str, Any]]:
        """Extract geometry optimization convergence from input variables."""
        convergence = []

        # Energy tolerance
        tolmxde = source.get('tolmxde')  # in hartree
        if tolmxde is not None:
            convergence.append({
                'convergence_parameter_name': 'energy',
                'threshold_type': 'absolute',
                'convergence_threshold': tolmxde.to('joule'),
                'convergence_threshold_unit': 'joule'
            })

        # Force tolerance
        tolmxf = source.get('tolmxf')  # in hartree/bohr
        if tolmxf is not None:
            convergence.append({
                'convergence_parameter_name': 'force',
                'threshold_type': 'maximum',
                'convergence_threshold': tolmxf.to('newton'),
                'convergence_threshold_unit': 'newton'
            })

        return convergence

# schema_packages/abinit.py
# Add mapping for convergence method
add_mapping_annotation(
    workflow.geometry_optimization.GeometryOptimizationMethod.convergence,
    OUT_KEY,
    ('get_geometry_convergence', ['.@'])
)

# Add mappings for WorkflowConvergenceTarget fields
add_mapping_annotations(
    (workflow.general.WorkflowConvergenceTarget.convergence_parameter_name,
     OUT_KEY, '.convergence_parameter_name'),
    (workflow.general.WorkflowConvergenceTarget.convergence_threshold,
     OUT_KEY, '.convergence_threshold'),
    (workflow.general.WorkflowConvergenceTarget.threshold_type,
     OUT_KEY, '.threshold_type'),
    (workflow.general.WorkflowConvergenceTarget.threshold_unit,
     OUT_KEY, '.convergence_threshold_unit')
)
```

---

## Testing Your Migration

### 1. Basic Test: Check Data Extraction

```python
def test_convergence_extraction():
    """Test that convergence methods return correct structure."""
    parser = YourParser()
    source = {...}  # Your test data

    convergence = parser.get_geometry_convergence(source)

    assert isinstance(convergence, list)
    assert len(convergence) > 0

    for target in convergence:
        assert 'convergence_parameter_name' in target
        assert 'threshold_type' in target
        assert 'convergence_threshold' in target
        assert 'convergence_threshold_unit' in target
        assert target['convergence_parameter_name'] in ['energy', 'force', 'potential', 'charge', 'density']
        assert target['threshold_type'] in ['absolute', 'relative', 'maximum', 'rms', 'residuum']
```

### 2. Integration Test: Check Schema Population

```python
def test_workflow_convergence_populated(test_archive):
    """Test that WorkflowConvergenceTarget is populated correctly."""
    # Parse with your parser
    parser.parse(test_file, archive=test_archive)

    # Check workflow section
    workflow = test_archive.data.workflow2
    assert workflow is not None
    assert workflow.method is not None
    assert workflow.method.convergence is not None
    assert len(workflow.method.convergence) > 0

    # Check individual targets
    for target in workflow.method.convergence:
        assert target.convergence_parameter_name is not None
        assert target.threshold_type is not None
        assert target.convergence_threshold is not None
        assert target.threshold_unit is not None
```

### 3. Unit Test: Check Unit Conversions

```python
def test_unit_conversions():
    """Test that units are converted correctly."""
    parser = YourParser()
    source = {
        'force_tolerance': 0.001 * ureg('eV/angstrom')  # Example input
    }

    convergence = parser.get_geometry_convergence(source)
    force_target = next(t for t in convergence if t['convergence_parameter_name'] == 'force')

    # Should be converted to SI (newton)
    assert force_target['convergence_threshold_unit'] == 'newton'
    expected_value = (0.001 * ureg('eV/angstrom')).to('newton').magnitude
    assert np.isclose(force_target['convergence_threshold'], expected_value)
```

---

## Troubleshooting

### Issue: "WorkflowConvergenceTarget has no attribute..."

**Cause**: Missing mapping annotation for WorkflowConvergenceTarget fields.

**Fix**: Ensure you have all four field mappings:
```python
add_mapping_annotations(
    (workflow.general.WorkflowConvergenceTarget.convergence_parameter_name, ...),
    (workflow.general.WorkflowConvergenceTarget.convergence_threshold, ...),
    (workflow.general.WorkflowConvergenceTarget.threshold_type, ...),
    (workflow.general.WorkflowConvergenceTarget.threshold_unit, ...)
)
```

### Issue: "convergence is None in workflow.method"

**Cause**: Your `get_geometry_convergence()` method isn't being called or isn't mapped.

**Fix**: Check that:
1. Method exists and returns a list (even if empty)
2. Mapping annotation exists for the convergence subsection
3. Method signature matches: `def get_geometry_convergence(self, source: dict[str, Any]) -> list[dict[str, Any]]`

### Issue: "Unit conversion error"

**Cause**: Threshold value already has units or wrong type.

**Fix**: Ensure threshold is a plain float (no units attached):
```python
# CORRECT:
'convergence_threshold': force_value.to('newton').magnitude  # ✓ float

# WRONG:
'convergence_threshold': force_value.to('newton')  # ✗ still a Quantity
```

Wait, actually reviewing the Exciting parser code, they DO pass quantities:
```python
'convergence_threshold': quantity[1].to(info_['unit']),  # This is a Quantity
```

So the threshold CAN be a Quantity object. The mapping parser will handle the conversion.

### Issue: "Convergence targets empty but data exists"

**Cause**: Logic error in extraction method (e.g., wrong dict keys, missing data checks).

**Fix**: Add defensive checks and logging:
```python
def get_geometry_convergence(self, source: dict[str, Any]) -> list[dict[str, Any]]:
    convergence = []

    force_tolerance = source.get('force_tolerance')
    if force_tolerance is None:
        self.logger.warning('No force tolerance found in source')
        return convergence  # Return empty list, not None

    # ... extract convergence
    return convergence
```

---

## Migration Checklist

For each parser migration:

- [ ] Identify which convergence fields the old parser populated
- [ ] Create `get_geometry_convergence()` method (if geometry opt)
- [ ] Create `get_single_point_convergence()` method (if SCF)
- [ ] Create `get_scf_steps()` method (if SCF iteration data)
- [ ] Add mapping annotations for convergence subsections
- [ ] Add mapping annotations for WorkflowConvergenceTarget fields
- [ ] Add mapping annotations for SCFSteps fields (if applicable)
- [ ] Remove or comment out old convergence tolerance mappings
- [ ] Write/update tests for convergence extraction
- [ ] Verify units are converted to SI
- [ ] Test with real parser output files
- [ ] Document any parser-specific quirks

---

## Reference Files

**Exciting Parser Implementation**:
- Schema: `src/nomad_simulation_parsers/schema_packages/exciting.py`
- Parser: `src/nomad_simulation_parsers/parsers/exciting/parser.py`
- Branch: `schema_update_convergence_targets`

**Schema Definitions**:
- `nomad-simulations/src/nomad_simulations/schema_packages/workflow/general.py` (WorkflowConvergenceTarget)
- `nomad-simulations/src/nomad_simulations/schema_packages/workflow/geometry_optimization.py` (GeometryOptimizationMethod)
- `nomad-simulations/src/nomad_simulations/schema_packages/outputs.py` (SCFSteps)

---

*Migration guide by: Claude Code*
*Based on: Exciting parser (schema_update_convergence_targets branch)*
*Last Updated: 2026-01-29*
