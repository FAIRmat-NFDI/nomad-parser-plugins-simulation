# CLI Parser Testing Framework

Exhaustive first-order validation of all parsers using `nomad parse` CLI command.

## Overview

This framework provides automated testing of all 15 parsers by:
1. Parsing test files via `nomad parse` command
2. Validating archive structure and content
3. Checking workflow types and program identification
4. Verifying required data sections are populated

## Files

- **`parser_test_matrix.yaml`** - Defines all test cases with expected outcomes
- **`test_cli_parsing.py`** - Pytest framework for CLI validation

## Running Tests

### Test All Parsers
```bash
uv run pytest tests/cli/test_cli_parsing.py -v
```

### Test Specific Parser
```bash
uv run pytest tests/cli/test_cli_parsing.py -k "abinit"
```

### Test with CLI Marker
```bash
uv run pytest -m cli
```

## Test Matrix Structure

The `parser_test_matrix.yaml` defines test cases for each parser:

```yaml
parsers:
  parser_name:
    test_cases:
      - name: "Test case description"
        mainfile: "path/to/test/file"
        expected_workflow: "SinglePoint|GeometryOptimization|MolecularDynamics|Phonon"
        required_fields:
          - "data.model_system"
          - "data.outputs"
          - "workflow2.method"  # for workflows
        skip: false  # Optional, defaults to false
```

## Validation Checks

Each test performs first-order validation:

1. **Parser doesn't crash** - Exit code 0
2. **Archive created** - Non-null archive returned
3. **Simulation data present** - `data` section populated
4. **Program identified** - Program name/version extracted
5. **Workflow type correct** - Matches expected workflow (if workflow2 exists)
6. **Required fields present** - All specified fields populated

## Test Results (Initial Run)

| Parser          | Status | Notes                                    |
|-----------------|--------|------------------------------------------|
| abinit          | PASS   | 2/2 test cases passing                   |
| ams             | -      | Not yet tested                           |
| crystal         | -      | Not yet tested                           |
| exciting        | FAIL   | JSON parsing error in archive output     |
| fhiaims         | -      | Not yet tested                           |
| gpaw            | PASS   | 1/1 test case passing                    |
| gromacs         | -      | Not yet tested                           |
| h5md            | -      | Not yet tested                           |
| lammps          | -      | Not yet tested                           |
| octopus         | -      | Not yet tested                           |
| phonopy         | FAIL   | AttributeError: 'PhonopyAtoms' no 'get_cell' |
| quantumespresso | -      | Not yet tested                           |
| vasp            | PASS   | 1/1 test case passing                    |
| wannier90       | -      | Not yet tested                           |
| yambo           | SKIP   | No test data available                   |

## Known Issues

### Exciting Parser
**Error:** JSON parsing fails when extracting archive
**Status:** Needs investigation - may be archive output format issue

### Phonopy Parser
**Error:** `AttributeError: 'PhonopyAtoms' object has no attribute 'get_cell'`
**Status:** Parser code bug - incompatible with current phonopy version

## Adding New Test Cases

1. Add test case to `parser_test_matrix.yaml` under appropriate parser:
   ```yaml
   - name: "Your test description"
     mainfile: "tests/data/parser/path/to/file"
     expected_workflow: "WorkflowType"
     required_fields:
       - "data.model_system"
       - "data.outputs"
   ```

2. Add test data to `tests/data/parser/` if needed

3. Run test:
   ```bash
   uv run pytest tests/cli/test_cli_parsing.py::test_cli_parsing -k "your_test"
   ```

## Design Decisions

**Simple Required Fields**
We validate top-level sections (`data.model_system`, `data.outputs`) rather than deep nested fields. This provides first-order validation without being brittle to schema changes.

**Workflow Validation Optional**
Not all parsers create `workflow2` sections. Workflow type validation only runs if the section exists.

**Timeout Protection**
All CLI parsing has a 60-second timeout to prevent hanging tests.

## Future Enhancements

1. Add deep field validation for critical parsers (VASP, Quantum ESPRESSO)
2. Validate units and numerical values
3. Test error handling (malformed files)
4. Performance benchmarking (large file parsing speed)
5. Memory usage tracking
6. Coverage reporting per parser
7. Integration with real-world samples from `.tests/`

## Related

- Unit tests: `tests/parsers/test_*_parser.py` - More detailed parser-specific tests
- Real-world tests: `.tests/pr-150-convergence/` - Production data validation
- GUI tests: `tests/gui/` - Upload workflow testing (to be implemented)
