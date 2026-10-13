"""
Pytest fixtures for GUI integration tests.

Manages NOMAD dev server lifecycle and Playwright browser instances.
"""

import subprocess
import time
from pathlib import Path

import pytest
import requests


@pytest.fixture(scope='session')
def nomad_dev_server():
    """
    Start NOMAD dev server for GUI testing.

    Requires:
    - Docker services running (docker compose up -d)
    - nomad.yaml configured for development mode
    - Backend dependencies installed

    Yields the base URL of the NOMAD server.
    """
    # Check if server is already running
    base_url = 'http://localhost:8000'
    try:
        response = requests.get(f'{base_url}/', timeout=2)
        # Any HTTP response (including 404 with JSON) indicates server is running
        if response.status_code in [200, 307, 404]:
            print('NOMAD dev server already running')
            yield base_url
            return
    except requests.exceptions.RequestException:
        pass

    # Start server
    root_dir = Path(__file__).resolve().parents[4]  # Navigate to workspace root
    process = subprocess.Popen(
        ['uv', 'run', 'nomad', 'admin', 'run', 'appworker'],
        cwd=root_dir,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    # Wait for server to be ready (max 60 seconds)
    max_wait = 60
    start_time = time.time()
    while time.time() - start_time < max_wait:
        try:
            response = requests.get(f'{base_url}/', timeout=2)
            # Any HTTP response (including 404 with JSON) indicates server is running
            if response.status_code in [200, 307, 404]:
                print(f'NOMAD dev server started (PID: {process.pid})')
                break
        except requests.exceptions.RequestException:
            time.sleep(1)
    else:
        process.kill()
        pytest.fail('NOMAD dev server failed to start within 60 seconds')

    yield base_url

    # Cleanup
    process.terminate()
    try:
        process.wait(timeout=10)
    except subprocess.TimeoutExpired:
        process.kill()
    print('NOMAD dev server stopped')


@pytest.fixture(scope='session')
def browser_context_args(browser_context_args):
    """Configure Playwright browser context."""
    return {
        **browser_context_args,
        'viewport': {'width': 1920, 'height': 1080},
        'ignore_https_errors': True,
    }


@pytest.fixture
def upload_page(page, nomad_dev_server):
    """
    Navigate to NOMAD upload page.

    Returns Playwright page object ready for upload testing.
    """
    page.goto(f'{nomad_dev_server}/gui/user/uploads')
    return page
