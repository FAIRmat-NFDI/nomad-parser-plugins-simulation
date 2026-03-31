# GUI Upload Workflow Testing Framework

Browser-based integration testing for all parsers using Playwright to validate the complete upload workflow through the NOMAD web interface.

## Overview

This framework provides automated GUI testing by:
1. Starting NOMAD dev server (backend + frontend)
2. Launching automated browser via Playwright
3. Uploading test files through the web interface
4. Validating archive processing and display
5. Checking metadata extraction and visualization

## Prerequisites

### 1. Docker Services

```bash
# Start required services (MongoDB, Elasticsearch, Temporal)
docker compose up -d

# Verify services are running
docker compose ps
```

### 2. NOMAD Configuration

Ensure `nomad.yaml` has development mode enabled:

```yaml
services:
  mode: development
```

### 3. Playwright Dependencies

```bash
# Install GUI testing dependencies (Playwright + pytest-playwright)
uv pip install -e ".[gui]"

# Install browser binaries (~300MB for Chromium)
uv run playwright install chromium

# Optional: Install all browsers (Chromium, Firefox, WebKit ~500MB total)
uv run playwright install
```

## Running Tests

### Run All GUI Tests

```bash
# Requires dev server to be running or fixture will start it
uv run pytest tests/gui/ -v -m gui
```

### Run Specific Parser

```bash
uv run pytest tests/gui/test_upload_workflow.py -k "vasp" -v
```

### Run Basic Smoke Tests

```bash
# Test server accessibility and page load only
uv run pytest tests/gui/test_upload_workflow.py::test_gui_server_accessible -v
uv run pytest tests/gui/test_upload_workflow.py::test_upload_page_loads -v
```

### Skip GUI Tests in CI

```bash
# Run all tests except GUI (useful for CI)
uv run pytest -m "not gui"
```

## Test Structure

### Files

- **`test_upload_workflow.py`** - Main upload and validation tests
- **`conftest.py`** - Pytest fixtures for server and browser management
- **`README.md`** - This documentation

### Test Workflow

Each GUI test performs the following steps:

1. **Navigate to upload page** - `http://localhost:8000/gui/user/uploads`
2. **Upload test file** - Select file via file input
3. **Wait for processing** - Monitor status until complete (max 120s timeout)
4. **Verify entry appears** - Check upload list contains new entry
5. **Open entry detail** - Click to view archive
6. **Validate archive display** - Check sections rendered correctly
7. **Verify program info** - Program name and version shown
8. **Check workflow type** - Workflow badge displays correct type

### Validation Checks

- [ ] File upload succeeds (no error messages)
- [ ] Processing completes (status indicator shows done)
- [ ] Entry appears in upload list
- [ ] Archive detail page loads
- [ ] Archive sections rendered (Overview, Data, Workflow)
- [ ] Program identification displayed
- [ ] Workflow type badge shows correct type

## Server Management

### Automatic Server Startup

The `nomad_dev_server` fixture automatically:
- Checks if server already running (port 8000)
- Starts server if needed via `nomad admin run appworker`
- Waits for server to be ready (max 60s)
- Stops server when tests complete

### Manual Server Management

If you prefer to manage the server manually:

```bash
# Terminal 1: Start backend
uv run nomad admin run appworker

# Terminal 2: Start frontend (if testing full stack)
cd packages/nomad-FAIR/gui
yarn start

# Terminal 3: Run tests
uv run pytest tests/gui/ -v -m gui
```

## Troubleshooting

### Server Won't Start

**Error:** `NOMAD dev server failed to start within 60 seconds`

**Solutions:**
1. Check Docker services running: `docker compose ps`
2. Verify nomad.yaml has `services.mode: development`
3. Check port 8000 not already in use: `lsof -i :8000`
4. Review server logs for errors

### Playwright Not Found

**Error:** `Failed to spawn: playwright`

**Solution:**
```bash
uv sync --extra dev
uv run playwright install
```

### Browser Doesn't Launch

**Error:** Playwright browser fails to start

**Solutions:**
1. Install system dependencies: `uv run playwright install-deps`
2. Try different browser: `pytest --browser chromium` (or firefox/webkit)

### File Upload Fails

**Error:** File input not found or upload doesn't trigger

**Solutions:**
1. Check NOMAD GUI version compatibility
2. Verify test file paths are correct (relative to project root)
3. Update selectors if GUI structure changed:
   - `[data-testid="upload-button"]`
   - `input[type="file"]`

### Test Timeouts

**Error:** `wait_for_selector` timeout after 120s

**Solutions:**
1. Increase timeout for large files: modify `timeout=120000` in test
2. Check parser actually completes processing
3. Verify server has sufficient resources
4. Check for parser bugs causing hangs

## Test Coverage

Reuses CLI test matrix (`tests/cli/parser_test_matrix.yaml`):
- Same 20 test cases
- Same parsers (abinit, ams, crystal, etc.)
- Same test files
- Validates GUI workflow instead of CLI

**Coverage:** 17/20 tests (3 skipped with parser bugs)

## Browser Options

Playwright supports multiple browsers:

```bash
# Chromium (default)
uv run pytest tests/gui/ --browser chromium

# Firefox
uv run pytest tests/gui/ --browser firefox

# WebKit (Safari engine)
uv run pytest tests/gui/ --browser webkit

# Headed mode (see browser window)
uv run pytest tests/gui/ --headed
```

## Debugging

### Visual Debugging

```bash
# Run with headed browser and slow motion
uv run pytest tests/gui/ --headed --slowmo 1000

# Take screenshots on failure (automatic with pytest-playwright)
uv run pytest tests/gui/ --screenshot on-failure

# Record video of test execution
uv run pytest tests/gui/ --video on
```

### Playwright Inspector

```bash
# Launch with debugger
PWDEBUG=1 uv run pytest tests/gui/test_upload_workflow.py::test_upload_workflow[vasp] -s
```

## Design Decisions

**Server Lifecycle:** Session-scoped fixture starts server once for all tests, improving performance.

**Selector Strategy:** Uses `data-testid` attributes where available, falls back to text selectors. May need updates if GUI changes.

**Timeout Values:**
- Page load: 30s (sufficient for upload confirmation)
- Processing: 120s (handles slow parsers)
- Server startup: 60s (allows for cold start)

**Test Isolation:** Each test uploads a file but doesn't clean up uploads. For production testing, add cleanup in fixture teardown.

**Reuse Test Matrix:** GUI tests use same YAML matrix as CLI tests for consistency. Any test case changes apply to both.

## Future Enhancements

1. Add cleanup of test uploads after each test
2. Test error handling (malformed files)
3. Test batch uploads (multiple files)
4. Validate specific archive fields (beyond basic display)
5. Test search and filtering of uploaded entries
6. Test download functionality
7. Performance benchmarking (upload → processing → display time)
8. Cross-browser compatibility matrix
9. Mobile viewport testing
10. Accessibility testing (ARIA labels, keyboard navigation)

## Related

- CLI tests: `tests/cli/` - Command-line integration testing
- Unit tests: `tests/parsers/` - Parser-specific detailed tests
- Test matrix: `tests/cli/parser_test_matrix.yaml` - Shared test case definitions
