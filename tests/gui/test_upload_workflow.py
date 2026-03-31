"""
GUI upload workflow integration tests.

Tests parser functionality through the NOMAD web interface by uploading
test files and validating the archive display.
"""

from pathlib import Path
from typing import Any

import pytest
import yaml


def load_test_matrix() -> dict[str, Any]:
    """Load parser test matrix from YAML file."""
    matrix_path = Path(__file__).parent.parent / 'cli' / 'parser_test_matrix.yaml'
    with open(matrix_path) as f:
        return yaml.safe_load(f)


def generate_gui_test_cases():
    """Generate pytest parameters from test matrix for GUI tests."""
    matrix = load_test_matrix()
    test_cases = []

    for parser_name, parser_config in matrix['parsers'].items():
        for test_case in parser_config['test_cases']:
            # Skip cases marked for skipping
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


@pytest.mark.gui
@pytest.mark.slow
@pytest.mark.parametrize('parser_name,test_case', generate_gui_test_cases())
def test_upload_workflow(upload_page, parser_name: str, test_case: dict[str, Any]):
    """
    Test file upload and archive display through GUI.

    Validates:
    1. File upload succeeds
    2. Processing completes
    3. Archive entry appears in upload list
    4. Entry detail page loads
    5. Archive data displayed correctly
    6. Program identification shown
    7. Workflow type badge displayed
    """
    mainfile = test_case['mainfile']
    expected_workflow = test_case['expected_workflow']
    file_path = Path(mainfile)

    # Verify file exists
    assert file_path.exists(), f"Test file not found: {mainfile}"

    # Step 1: Upload file
    # Wait for upload button to be visible
    upload_button = upload_page.locator('[data-testid="upload-button"]')
    upload_button.wait_for(state='visible', timeout=10000)

    # Trigger file input
    file_input = upload_page.locator('input[type="file"]')
    file_input.set_input_files(str(file_path.absolute()))

    # Wait for upload to complete (check for success message or entry in list)
    upload_page.wait_for_selector(
        f'text="{file_path.name}"', state='visible', timeout=30000
    )

    # Step 2: Wait for processing to complete
    # Look for status indicator showing processing is done
    # This depends on NOMAD GUI implementation - might need adjustment
    upload_page.wait_for_selector(
        '[data-testid="processing-complete"]', state='visible', timeout=120000
    )

    # Step 3: Verify entry appears in upload list
    entry_row = upload_page.locator(f'text="{file_path.name}"')
    assert entry_row.is_visible(), f"Upload entry not found for {file_path.name}"

    # Step 4: Click on entry to open detail page
    entry_row.click()
    upload_page.wait_for_load_state('networkidle')

    # Step 5: Validate archive sections rendered
    # Check that main archive sections are displayed
    archive_sections = ['Overview', 'Data', 'Workflow']
    for section in archive_sections:
        section_header = upload_page.locator(f'text="{section}"')
        # Section might not always be visible, just check it exists in DOM
        section_header.wait_for(state='attached', timeout=5000)

    # Step 6: Verify program identification
    # Look for program name in the archive display
    program_section = upload_page.locator('[data-testid="program-info"]')
    assert program_section.is_visible(), "Program information not displayed"

    # Step 7: Verify workflow type badge
    # Check for workflow type indicator
    workflow_badge = upload_page.locator(f'text="{expected_workflow}"')
    assert (
        workflow_badge.is_visible()
    ), f"Workflow type '{expected_workflow}' not displayed"


@pytest.mark.gui
def test_gui_server_accessible(nomad_dev_server):
    """Verify NOMAD dev server is accessible."""
    import requests

    response = requests.get(f'{nomad_dev_server}/alive')
    assert response.status_code == 200
    assert response.json().get('alive') is True


@pytest.mark.gui
def test_upload_page_loads(upload_page):
    """Verify upload page loads successfully."""
    assert upload_page.title() != '', 'Upload page title should not be empty'
    assert upload_page.url.endswith('/user/uploads'), 'Should be on uploads page'
