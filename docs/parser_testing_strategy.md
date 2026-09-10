## Fundamentals

### Different test types serve different purposes and catch different categories of bugs

1. **Unit Tests** catch logic errors early and cheaply  
2. **Integration Tests** catch interface mismatches and coordination issues  
3. **E2E Tests** catch real-world edge cases and validate user-facing behavior

#### Why the testing pyramid

* **Cost increases up the pyramid**: E2E tests are 10-100x slower and harder to maintain  
* **Feedback speed decreases up the pyramid**: Unit tests give instant feedback, E2E tests take minutes  
* **Debugging difficulty increases up the pyramid**: A failing E2E test could be caused by any of 20 components  
* **Coverage breadth decreases up the pyramid**: 1000 unit tests can cover edge cases that would take 100,000 E2E tests

### Unit Tests

Definition: A **unit test** verifies a single function, method, or class in isolation.

#### Purpose

* **Fast feedback**: Developers get instant feedback during coding (\< 10ms execution)  
* **Pinpoint failures: When a unit test fails, you know exactly which function broke**  
* **Enable refactoring**: Safe to restructure code when comprehensive unit tests exist  
* **Document behavior**: Unit tests serve as executable documentation of what each function should do  
* **Catch logic errors: Off-by-one errors, incorrect conditionals, edge cases in algorithms**

#### When to write

* For any function with conditional logic (if/else, loops)  
* For data transformation functions (parsing, formatting, validation)  
* For calculations and algorithms  
* For utility functions used across the codebase

#### Characteristics

* Have **no external dependencies** (no file I/O, no database, no network)  
* Execute in **\< 10ms**  
* Test **one thing** at a time  
* Use **mocks/stubs** instead of actual dependencies  
* Should be **60-90% of the test suite**

**Example:** Testing a pattern extraction function

````python
def test_pattern_extraction_from_text():
    """Test regex pattern matching for key-value extraction"""
    pattern = r"(\\w\+)\\s\*\=\\s\*([-\\d.]\+)"
    test_line = "temperature = 300.5"
    match = re.match(pattern, test_line)

    assert match.group(1) == "temperature"
    assert float(match.group(2)) == 300.5
````

**Example:** Testing unit conversion logic  

````python
def test_value_conversion_with_units():
    """Test that numerical values are correctly converted with unit objects"""
    converter = UnitConverter()
    converter.set_unit_system("metric")

    raw_value = 100.0
    converted = converter.apply_unit(raw_value, "length")

    assert converted.magnitude == 100.0
    assert converted.units == ureg.meter
````

### Integration Tests

Definition: An **integration test** verifies that multiple units work together correctly through their real interfaces.

#### Purpose

* **Catch interface mismatches**: Unit A might work alone, but pass wrong data type to Unit B  
* **Verify data flow**: Ensure data transformations work correctly through the pipeline  
* **Test realistic scenarios**: Components might work individually but fail when coordinating  
* **Validate configuration**: Test that components are correctly wired together  
* **Catch integration bugs**: Timing issues, resource sharing, state management problems

#### When to write

* When multiple classes must coordinate (e.g., Parser \+ Validator \+ Writer)  
* When testing file discovery and multi-file handling  
* When validating data flows through processing pipelines  
* When testing configuration and dependency injection  
* For testing public APIs that orchestrate multiple components

#### Characteristics

* May involve **controlled file I/O** (temporary files, minimal data)  
* Test **class interactions** and data flow  
* Execute in **10-100ms**  
* Use **real implementations** where practical (but with minimal/synthetic data)  
* Should comprise **20-30% of the test suite**

**Example:** Parser coordination test

````python
@pytest.fixture
def mock_components(mocker):
    """Mock all components to test orchestration logic only"""
    return {
        "reader": mocker.Mock(spec=FileReader),
        "validator": mocker.Mock(spec=DataValidator),
        "writer": mocker.Mock(spec=OutputWriter),
    }

def test_pipeline_orchestration(mock_components):
    """Test how Pipeline coordinates components without I/O"""
    pipeline = DataPipeline()
    pipeline._reader = mock_components["reader"]
    pipeline._validator = mock_components["validator"]
    pipeline._writer = mock_components["writer"]

    # Setup mock behavior
    mock_components["reader"].read_chunk.return_value = {
        "data": np.zeros((10, 3)),
        "metadata": {"timestamp": 0.0},
    }
    mock_components["validator"].is_valid.return_value = True

    # Execute pipeline
    result = pipeline.process_chunk(0)

    # Verify coordination
    mock_components["reader"].read_chunk.assert_called_once_with(0)
    mock_components["validator"].is_valid.assert_called_once()
    mock_components["writer"].write.assert_called_once()
    assert result.success is True
````

**Example:** File discovery and matching  

````python
@pytest.mark.parametrize(
    "main_file, auxiliary_files, expected_match",
    [
        # Strategy 1: Exact prefix match
        ("config.main", ["main_data.txt", "other.txt"], "main_data.txt"),
        # Strategy 2: Token-based similarity
        (
            "config.project_v2",
            ["proj_v2_results.csv", "unrelated.csv"],
            "proj_v2_results.csv",
        ),
        # Strategy 3: Fallback to first file
        ("config.unique", ["file1.txt", "file2.txt"], "file1.txt"),
    ],
)
def test_auxiliary_file_matching(main_file, auxiliary_files, expected_match):
    """Test file matching with various strategies"""
    finder = FileFinder()
    finder.logger = LOGGER

    result = finder.find_best_match(auxiliary_files, main_file)

    assert len(result) == 1
    assert result[0] == expected_match
````

**Example:** Minimal synthetic data  

````python
@pytest.fixture
def minimal_config_file(tmp_path):
    """Generate minimal valid configuration content"""
    content = """# Configuration v1.0
setting1 = value1
setting2 = 42
setting3 = true
"""
    file = tmp_path / "config.txt"
    file.write_text(content)
    return file

def test_config_parsing_with_minimal_data(minimal_config_file):
    """Test parsing with controlled synthetic data"""
    parser = ConfigParser()
    parser.parse(minimal_config_file)

    assert parser.get("setting1") == "value1"
    assert parser.get("setting2") == 42
    assert parser.get("setting3") is True
````

### End-to-End (E2E) Tests

**Definition**: An end-to-end test verifies the complete processing pipeline with realistic production data.

#### Purpose

* **Validate real-world behavior**: Catch edge cases that only appear with real data  
* **Test the full user journey**: From input files to final output  
* **Integration with external systems**: Database, file system, network in realistic configuration  
* **Catch regression bugs**: Ensure new changes don't break existing functionality  
* **Build confidence**: Demonstrate the system works end-to-end before deployment  
* **Document system capabilities**: Show what file formats and scenarios are supported

#### When to write

* For each major file format or input type supported  
* For critical user workflows that must always work  
* When adding support for a new data source  
* To validate bug fixes with the exact data that caused the bug  
* As acceptance tests for new features

#### When **NOT** to write

* For testing edge cases (use unit tests instead)  
* For testing error conditions (use unit/integration tests)  
* For every possible combination (combinatorial explosion)  
* As the first test you write (start with unit tests)

#### Characteristics

* Use **real production data files**  
* Test the **entire workflow** from file → processed output  
* Execute in **100ms-1s** (or more for complex cases)  
* Are **expensive to run** and maintain  
* Should comprise **5-10% of the test suite**  
* **Should be marked for selective execution**

**Example:** Complete pipeline with real data  

````python
@pytest.mark.slow  
@pytest.mark.e2e  
def test_complete_data_processing_pipeline():  
    """E2E test: Process complete dataset with all features"""  
    output_container = OutputContainer()  
    parser = DataParser()

    parser.parse("tests/data/realistic_dataset/main_file.txt", output_container)  
    apply_all_validations(output_container)

    # High-level assertions only - don't verify every detail  
    assert output_container.data is not None  
    assert len(output_container.data.records) > 0  
    assert output_container.data.records[0].values is not None  
    assert output_container.metadata.processing_status == "success"
````

**Example:** Format-specific E2E test  

````python
@pytest.mark.slow
@pytest.mark.e2e
def test_json_format_end_to_end():
    """E2E test for JSON input format"""
    output = OutputContainer()
    parser = DataParser()

    parser.parse("tests/data/formats/sample.json", output)

    # Verify format-specific features
    assert output.data.format == "json"
    assert output.data.has_nested_structures is True  
    assert len(output.data.records) == 100  # Known from fixture
````

### Test organization structure

tests/  
├── unit/                          \# Fast, no I/O (60-90% of tests)  
│   ├── test\_parsers.py           \# Individual parser components  
│   ├── test\_validators.py        \# Validation logic  
│   ├── test\_transformers.py      \# Data transformation functions  
│   ├── test\_utils.py             \# Utility functions  
│   └── test\_converters.py        \# Unit conversion, type conversion  
│  
├── integration/                   \# Minimal file I/O (20-30% of tests)  
│   ├── test\_pipeline\_coordination.py  \# Multi-component workflows  
│   ├── test\_file\_discovery.py    \# File matching and discovery  
│   └── fixtures/                  \# Synthetic test data (\< 1KB each)  
│       ├── minimal\_config.txt  
│       ├── minimal\_data.csv  
│       └── minimal\_schema.json  
│  
└── e2e/                           \# Real data (5-10% of tests)  
    ├── test\_supported\_formats.py  \# One test per format  
    ├── test\_representative\_workflows.py  \# Common user scenarios  
    └── data/                      \# Real production samples  
        ├── json\_format/  
        │   └── sample\_v1.json  
        ├── csv\_format/  
        │   └── sample\_export.csv  
        └── xml\_format/  
            └── sample\_config.xml

### Specific testing patterns

1. **Synthetic Minimal Data**  
   **When to use**: Integration tests that need files but want to avoid I/O overhead

````python
@pytest.fixture
def minimal_data_file(tmp_path):
    """Generate minimal valid data content"""
    content = """# Data File v1.0
timestamp,value,status
0.0,100.0,ok
1.0,101.5,ok
2.0,99.8,ok
"""
    file = tmp_path / "data.csv"
    file.write_text(content)
    return file

@pytest.fixture
def minimal_frame_data():
    """Generate minimal valid frame structure"""
    return {
        "timestamp": 0,
        "count": 5,
        "bounds": {"min": -10.0, "max": 10.0},
        "records": [
            {"id": 1, "type": "A", "x": 0.0, "y": 0.0, "z": 0.0},
            {"id": 2, "type": "A", "x": 1.0, "y": 1.0, "z": 1.0},
            {"id": 3, "type": "B", "x": 2.0, "y": 2.0, "z": 2.0},
            {"id": 4, "type": "B", "x": 3.0, "y": 3.0, "z": 3.0},
            {"id": 5, "type": "B", "x": 4.0, "y": 4.0, "z": 4.0},
        ],
    }
````

2. **Mock objects and test doubles**
   **When to use**: Unit tests for complex orchestration without external dependencies

````python
class StubDataReader:
    """Test double for DataReader to avoid file I/O"""

    def __init__(self, records, metadata):
        self._records = records
        self._metadata = metadata

    def get_records(self, idx):
        return self._records[idx]

    def get_metadata(self, idx):
        return self._metadata[idx]

    def count(self):
        return len(self._records)

def test_data_processing_with_stub():
    """Test processing logic without real files"""
    stub_reader = StubDataReader(
        records=[np.zeros((100, 3))], metadata=[{"timestamp": 0.0, "source": "test"}]
    )

    processor = DataProcessor()
    processor.data_source = stub_reader

    # Test logic without file I/O
    processed = processor.process_batch(0)
    assert processed.shape == (100, 3)
    assert processed.metadata["source"] == "test"
````

3. **Parameterized tests**
   **When to use:** Testing the same logic with multiple input variations

````python
@pytest.mark.parametrize(
    "input_format,expected_type",
    [
        ("json", DataType.STRUCTURED),
        ("csv", DataType.TABULAR),
        ("xml", DataType.HIERARCHICAL),
        ("txt", DataType.UNSTRUCTURED),
    ],
)
def test_format_detection(input_format, expected_type):
    """Test format detection without separate test functions"""
    detector = FormatDetector()
    result = detector.detect(f"sample.{input_format}")
    assert result == expected_type

@pytest.mark.parametrize(
    "value,unit,expected",
    [
        (100, "meters", 100.0 * ureg.meter),
        (1000, "millimeters", 1.0 * ureg.meter),
        (1, "kilometers", 1000.0 * ureg.meter),
    ],
)
def test_unit_normalization(value, unit, expected):
    """Test all unit conversions in one test function"""
    converter = UnitConverter()
    result = converter.normalize(value, unit, target="meters")
    assert result == pytest.approx(expected)
````

4. **Schema validation**  
   **When to use:** Ensuring output conforms to expected structure

````python
def test_output_schema_validation():
    """Verify output conforms to expected schema"""
    output = ProcessedData(
        records=np.zeros((10, 3)), count=10, metadata={"source": "test"}
    )

    # Schema validation
    assert hasattr(output, "records")
    assert hasattr(output, "count")
    assert hasattr(output, "metadata")
    assert output.records.shape == (10, 3)
    assert output.count == 10
    assert isinstance(output.metadata, dict)

def test_nested_schema_validation():
    """Verify nested structures conform to schema"""
    container = DataContainer()
    container.header = Header(version="1.0", author="test")
    container.body = Body(records=[Record(id=1, value=100)])

    # Validate structure
    assert container.header.version == "1.0"
    assert len(container.body.records) == 1
    assert container.body.records[0].id == 1
````

### Performance Optimization

* **Fixture management: reuse extensive setup**  
  **Rationale:** Expensive setup operations (parsing large files, database connections) should be reused across tests

````python
@pytest.fixture(scope="session")  # Reuse across entire test session
def reference_dataset():
    """Parse once, use many times for validation"""
    container = DataContainer()
    processor = DataProcessor()
    processor.process("tests/fixtures/reference.dat", container)
    return container

@pytest.fixture
def dataset_copy(reference_dataset):
    """Fast deepcopy for test isolation"""
    return copy.deepcopy(reference_dataset)

def test_with_reused_data(dataset_copy):
    """Test can modify copy without affecting other tests"""
    dataset_copy.records[0].value = 999
    # Other tests still see original data
````

**Fixture scope levels:**

* **scope="function"** (default): New fixture for each test  
  * **scope="class"**: One fixture per test class  
  * **scope="module"**: One fixture per test file  
  * **scope="session"**: One fixture for entire test run

* **Selective test execution**  
  **Rationale**: Different development stages need different test coverage/speed trade-offs

````python
# pytest.ini or pyproject.toml
[tool.pytest.ini_options]
markers = [
    "unit: Fast unit tests (< 10ms each)",
    "integration: Tests requiring file I/O (< 100ms each)",
    "slow: End-to-end tests with real data (> 100ms each)",
    "e2e: Complete pipeline tests",
    "requires_network: Tests needing network access",
]
````

#### Execution strategies

````bash
  # Only fast tests during active development (~1-2 seconds)    
  pytest -m "not slow" -v  

  # Unit tests only (~1-2 seconds)    
  pytest tests/unit -v  

  # Integration tests (~10-30 seconds)    
  pytest tests/integration -v  

  # Everything except network tests (~1-5 minutes)    
  pytest -m "not requires_network" -v  

  # Full suite including slow tests (~5-15 minutes)    
  pytest -v  

  # Run until first failure (for debugging)    
  pytest -x --tb=short  

  # Run last failed tests only    
  pytest --lf  
````

#### Continuous testing workflow

* **Development Workflow:**

````bash
  # 1. During active development (watch mode - requires pytest-watch)  
  pytest-watch tests/unit  # Re-run on file changes  

  # 2. Before committing changes  
  pytest tests/unit -x --tb=short  # Stop on first failure, short traceback  

  # 3. Before pushing  
  pytest tests/unit tests/integration -v  

  # 4. Before creating PR  
  pytest -m "not slow" --cov --cov-report=term-missing
  ````

* **CI/CD Pipeline:**

````yaml
  # .github/workflows/tests.yml (or similar CI config)
  jobs:

   unit-tests:

     runs-on: ubuntu-latest

     steps:

       - name: Run unit tests

         run: pytest tests/unit -v --cov --cov-report=xml

         timeout-minutes: 5  # Fail if tests take too long

    integration-tests:

     runs-on: ubuntu-latest

     needs: unit-tests  # Only run if unit tests pass

     steps:

       - name: Run integration tests

         run: pytest tests/integration -v

         timeout-minutes: 10

    e2e-tests:

     runs-on: ubuntu-latest

     needs: integration-tests

     steps:

       - name: Run E2E tests

         run: pytest tests/e2e -v

         timeout-minutes: 15
````

* **Parallel Execution**  
  **Rationale:** Utilize multiple CPU cores to speed up test execution

````bash
  # Install pytest-xdist  

  pip install pytest-xdist  

  # Run tests in parallel  

  pytest -n auto  # Use all available CPU cores  
  pytest -n 4     # Use 4 worker processes  

  # Parallel execution with coverage (requires pytest-cov)  

  pytest -n auto --cov --cov-report=html
````

#### Important considerations

* Tests must be independent (no shared state)  
  * File-based tests need unique temporary directories  
  * Database tests need separate test databases  
  * \~30-50% speedup for CPU-bound tests

### Testing schedule?

| When | What | Expected Time | Purpose |
| :---- | :---- | :---- | :---- |
| On save (watch mode) | Unit tests for changed file | \< 1 second | Immediate feedback during coding |
| Before commit | All unit tests | 1-5 seconds | Ensure no unit-level regressions |
| On PR | Unit \+ Integration | 10-60 seconds | Validate changes don't break integration |
| Before merge | Full suite (fast)  | 1-5 minutes | Comprehensive validation minus slow tests |
| On main branch | Full suite including E2E | 5-15 minutes | Complete validation before deployment |
| Nightly | Full suite \+ slow E2E | 5-15 minutes | Catch issues with production-like data |

### Key principles Summary

1. **Isolation**  
   * Each test should be completely independent  
   * Tests that depend on execution order or shared state are brittle and hard to debug  
   * Use fixtures for setup, clean up after each test, avoid global state

````python
# BAD: Tests depend on execution order
def test_step1():
    global state
    state = initialize()

def test_step2():  # Fails if run alone
    process(state)

# GOOD: Each test is self-contained
def test_step1(initialized_state):
    result = process(initialized_state)
    assert result.is_valid()

def test_step2(initialized_state):
    result = transform(initialized_state)
    assert result.is_transformed()
````

2. **Speed**

   * Tests should be as fast as possible for their category  
   * Slow tests discourage running them frequently, reducing their value  
   * Avoid unnecessary I/O, use mocks, keep test data minimal

    **Target execution times:**

    * Unit tests: \< 10ms each  
    * Integration tests: 10-100ms each  
    * E2E tests: 100ms-1s each

3. **Clarity**

   * Test names and assertions should clearly describe intent  
   * Tests serve as documentation and must be understood 6 months later  
   * Use descriptive names, clear assertions, avoid testing implementation details

````python
   # BAD: Test name describes implementation
   def test_regex_pattern_matches():

      assert parser._position_regex.match("x y z")

   # GOOD: Test name describes behavior

   def test_extracts_coordinates_from_position_data():

      """Test that 3D coordinates are correctly extracted"""

      parser = CoordinateParser()

      result = parser.extract_coordinates("x=1.0 y=2.0 z=3.0")

      assert result == (1.0, 2.0, 3.0)
````

4. **Behavior over Implementation**

   * Test **what** the code does, not **how** it does it  
   * Implementation can change; behavior should remain stable  
   * Test public APIs, not private methods or internal state

````python
# BAD: Testing implementation details
def test_internal_cache_structure():
    parser = DataParser()
    parser._parse_file("data.txt")
    assert len(parser._cache) == 5  # Breaks if caching strategy changes

# GOOD: Testing observable behavior
def test_parsing_returns_expected_records():
    parser = DataParser()
    records = parser.parse("data.txt")
    assert len(records) == 5
    assert records[0].id == 1
````

5. **Maintainability**

   * Tests should be easy to update when requirements change  
   * Brittle tests become technical debt and get disabled  
   * Use fixtures, parametrization, helper functions; avoid duplication

````python
# BAD: Duplicated setup code
def test_case_1():
    parser = DataParser()
    parser.set_format("json")
    parser.set_encoding("utf-8")
    # ... test logic

def test_case_2():
    parser = DataParser()
    parser.set_format("json")
    parser.set_encoding("utf-8")
    # ... test logic

# GOOD: Shared setup via fixture
@pytest.fixture
def configured_parser():
    parser = DataParser()
    parser.set_format("json")
    parser.set_encoding("utf-8")
    return parser

def test_case_1(configured_parser):
    # ... test logic

def test_case_2(configured_parser):
    # ... test logic
````

6. **Selective Execution**

   * Tests should be categorized for selective running  
   * Not all tests need to run all the time  
   * Use markers to enable running subsets based on context

````python
# Mark expensive tests
@pytest.mark.slow
@pytest.mark.e2e
def test_full_pipeline_with_large_dataset():
    """This test takes 5 seconds - skip during development"""
    pass

# Mark tests requiring external resources
@pytest.mark.requires_network
def test_api_integration():
    """Requires network access - skip in offline mode"""
    pass

# Mark tests for specific features
@pytest.mark.feature_x
def test_new_feature_x():
    """Test for feature X - run when working on that feature"""
    pass
````

## Common testing anti-patterns to avoid

1. **Testing implementation details**

````python
   # BAD: Test breaks when refactoring internal structure
   def test_internal_method():
      obj._internal_helper_method()  # Don't test private methods

   # GOOD: Test public behavior
   def test_public_method():
      result = obj.process_data(input)
      assert result == expected_output
````

2. **Overly coupled tests**

````python
   # BAD: Test depends on exact internal state
   def test_processing():
      processor.process()
      assert len(processor._internal_cache) == 3
      assert processor._state == "processed"
      
   # GOOD: Test observable outcomes
   def test_processing():
      result = processor.process()
      assert result.success is True
      assert len(result.records) == 3
````

3. **Mystery guest pattern**

````python
   # BAD: Test depends on external file without making it clear
   def test_parse():
      parser.parse("data.txt")  # Where is this file? What's in it?

   # GOOD: Clear test data source
   def test_parse(tmp_path):
      test_file = tmp_path / "data.txt"
      test_file.write_text("id,value\\n1,100\\n2,200")
      result = parser.parse(test_file)
      assert len(result) == 2
````

4. **Asserting too much**

````python
   # BAD: Test is fragile and hard to maintain
   def test_output():
      result = process()
      assert result.field1 == value1
      assert result.field2 == value2
      assert result.field3 == value3
      # ... 20 more assertions

   # GOOD: Test one behavior per test
   def test_output_has_required_fields():
      result = process()
      assert hasattr(result, 'field1')
      assert hasattr(result, 'field2')
      
   def test_output_field1_value():
      result = process()
      assert result.field1 == expected_value1
````
