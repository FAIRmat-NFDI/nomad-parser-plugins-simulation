"""
Exhaustive CLI parsing tests for all parsers.

Tests all parsers using `nomad parse` command to validate:
- Parser doesn't crash
- Archive structure created
- Program identification
- Workflow type correctness
- Required fields populated

Uses parser_test_matrix.yaml to define test cases.
"""

import json
import subprocess
from pathlib import Path
from typing import Any

import pytest
import yaml


def load_test_matrix() -> dict[str, Any]:
    """Load parser test matrix from YAML file."""
    matrix_path = Path(__file__).parent / 'parser_test_matrix.yaml'
    with open(matrix_path) as f:
        return yaml.safe_load(f)


def generate_test_cases():
    """Generate pytest parameters from test matrix."""
    matrix = load_test_matrix()
    test_cases = []

    for parser_name, parser_config in matrix['parsers'].items():
        for test_case in parser_config['test_cases']:
            # Skip cases marked for skipping (e.g., YAMBO with no test data)
            if test_case.get('skip', False):
                continue

            test_id = f"{parser_name}::{test_case['name']}"
            test_cases.append(
                pytest.param(
                    parser_name,
                    test_case,
                    id=test_id,
                )
            )

    return test_cases


def parse_via_cli(mainfile: str) -> dict[str, Any] | None:
    """
    Parse a file using `nomad parse` CLI and return the archive as dict.

    Args:
        mainfile: Path to the file to parse

    Returns:
        Parsed archive as dictionary, or None if parsing failed
    """
    try:
        result = subprocess.run(
            ['uv', 'run', 'nomad', 'parse', mainfile, '--show-archive'],
            capture_output=True,
            text=True,
            timeout=60,  # 60 second timeout
            check=True,
        )

        # Parse JSON from output
        # The archive is output as JSON after any warning messages
        output_lines = result.stdout.strip().split('\n')

        # Find where JSON starts (first line with '{')
        json_start = None
        for i, line in enumerate(output_lines):
            if line.strip().startswith('{'):
                json_start = i
                break

        if json_start is None:
            return None

        json_str = '\n'.join(output_lines[json_start:])
        return json.loads(json_str)

    except subprocess.CalledProcessError as e:
        pytest.fail(f"nomad parse failed with exit code {e.returncode}:\n{e.stderr}")
    except subprocess.TimeoutExpired:
        pytest.fail(f"nomad parse timed out after 60 seconds")
    except json.JSONDecodeError as e:
        pytest.fail(f"Failed to parse archive JSON: {e}")

    return None


def get_nested_value(data: dict, path: str) -> Any:
    """
    Get value from nested dictionary using dot-notation path.

    Args:
        data: Dictionary to search
        path: Dot-notation path (e.g., 'data.model_system')

    Returns:
        Value at path, or None if not found
    """
    keys = path.split('.')
    value = data

    for key in keys:
        if isinstance(value, dict) and key in value:
            value = value[key]
        else:
            return None

    return value


@pytest.mark.cli
@pytest.mark.parametrize('parser_name,test_case', generate_test_cases())
def test_cli_parsing(parser_name: str, test_case: dict[str, Any]):
    """
    Test parser via CLI for first-order validation.

    Validates:
    1. Parser doesn't crash (exit code 0)
    2. Archive created
    3. Simulation data present
    4. Program identified
    5. Workflow type matches expected
    6. Required fields populated
    """
    mainfile = test_case['mainfile']
    expected_workflow = test_case['expected_workflow']
    required_fields = test_case['required_fields']

    # Parse via CLI
    archive = parse_via_cli(mainfile)

    # Check 1: Archive created
    assert archive is not None, f"Failed to create archive for {mainfile}"

    # Check 2: Simulation data present
    assert 'data' in archive, "Archive missing 'data' section"
    data = archive['data']
    assert data is not None, "data section is None"

    # Check 3: Program identified
    program = get_nested_value(archive, 'data.program')
    assert program is not None, "Program not identified"
    assert 'name' in program or 'version' in program, "Program missing name/version"

    # Check 4: Workflow type correct (if workflow2 exists)
    if 'workflow2' in archive and archive['workflow2']:
        workflow = archive['workflow2']
        if workflow and 'm_def' in workflow:
            workflow_type = workflow['m_def'].split('.')[-1]
            assert (
                workflow_type == expected_workflow
            ), f"Expected workflow {expected_workflow}, got {workflow_type}"

    # Check 5: Required fields present
    for field_path in required_fields:
        value = get_nested_value(archive, field_path)
        assert value is not None, f"Required field '{field_path}' not populated"


@pytest.mark.cli
def test_all_parsers_covered():
    """Verify all parsers have at least one test case."""
    matrix = load_test_matrix()
    parsers = list(matrix['parsers'].keys())

    # Expected parsers (from source code)
    expected_parsers = {
        'abinit',
        'ams',
        'crystal',
        'exciting',
        'fhiaims',
        'gpaw',
        'gromacs',
        'h5md',
        'lammps',
        'octopus',
        'phonopy',
        'quantumespresso',
        'vasp',
        'wannier90',
        'yambo',
    }

    covered_parsers = set(parsers)
    missing = expected_parsers - covered_parsers

    assert (
        not missing
    ), f"Parsers missing from test matrix: {', '.join(sorted(missing))}"


@pytest.mark.cli
def test_matrix_yaml_valid():
    """Verify test matrix YAML is valid and well-formed."""
    matrix = load_test_matrix()

    assert 'parsers' in matrix, "Matrix missing 'parsers' key"
    assert 'global_checks' in matrix, "Matrix missing 'global_checks' key"

    for parser_name, parser_config in matrix['parsers'].items():
        assert (
            'test_cases' in parser_config
        ), f"Parser {parser_name} missing 'test_cases'"

        for i, test_case in enumerate(parser_config['test_cases']):
            assert 'name' in test_case, f"{parser_name} test case {i} missing 'name'"
            assert (
                'mainfile' in test_case
            ), f"{parser_name} test case {i} missing 'mainfile'"
            assert (
                'expected_workflow' in test_case
            ), f"{parser_name} test case {i} missing 'expected_workflow'"
            assert (
                'required_fields' in test_case
            ), f"{parser_name} test case {i} missing 'required_fields'"
