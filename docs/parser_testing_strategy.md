## Checklists

Two lists to run through before a pull request is opened, one for a parser change and one for a schema change. Each item names the protocol layer it belongs to and links to the section that explains it; the rest of this document is the reasoning behind the items, and [`how_to/test_parsers.md`](how_to/test_parsers.md) is the protocol they follow.

### Parser: adding or changing a mapped quantity

The items follow the order of the work: the text is read first, then mapped, with any shared helper tested on its own; the last group is run once before the pull request is opened. Each item shows a minimal pair from this repository: first the implementation line under test (marked `# implementation`), then the test that pins it (marked `# test`), both trimmed to the lines that make the point. Where no parser has the test yet, the item says so and shows the shape.

**Text processing**

1. **Extraction.** In `test_reader.py`, write a few lines of output text that contain the quantity, run the reader on them, and check that the regex extracted the right value. If the quantity is assembled from several blocks of the output, write one such test per block ([layer 2](how_to/test_parsers.md#layer-2-source-extraction)). Only extraction is tested: the text goes into a file under `tmp_path`, the reader turns it into its dictionary, and the test compares entries of that dictionary with the numbers and strings visible in the text. The schema is not involved. Whether the value then lands in the right schema field is checked in item 4, on the archive the whole parser produces. ORCA:

````python
# implementation: parsers/orca/text_parser.py
Quantity('hf_type', r'Hartree-Fock type\s*HFTyp\s*\.+\s*(\w+)', convert=False),
Quantity('total_charge', rf'Total Charge\s*Charge\s*\.+\s*({re_float})', dtype=float),

# test: tests/parsers/orca/test_reader.py
@pytest.mark.unit
class TestOrcaReader:
    def test_extracts_scf_settings(self, tmp_path):
        source = _read(  # writes the fragment to tmp_path and runs OutReader on it
            tmp_path,
            SINGLE_POINT_SCF + 'SCF SETTINGS\n------------\n'
            'Hartree-Fock type      HFTyp           .... RHF\n'
            'Total Charge           Charge          .... 0\n'
            '------------\n',
        )
        scf = _scf_settings(source)
        assert scf['hf_type'] == 'RHF'
        assert scf['total_charge'] == 0.0
````

   H5MD and phonopy read HDF5 and YAML through libraries and have no reader layer of their own; for them this item is empty and the work starts at item 2.

**Mapping**

2. **Mapped value.** In `test_mapping.py`, call the transformer function directly on a small dictionary that stands in for the reader's output, and check the value, unit and shape it returns ([layer 3](how_to/test_parsers.md#layer-3-mapping-contract)). Each input dictionary with its expected result is one row of a `pytest.mark.parametrize` table, with an `ids` entry naming the case ([parameterized tests](#specific-testing-patterns)). Add rows for the quantity being absent, empty, zero or `False` wherever the source can produce them, because these are the inputs that mapping code most often drops silently. VASP, a pure function with a table of cases:

````python
# implementation: parsers/vasp/common.py
def functional_key_from_params(parameters: dict[str, Any]) -> str | None:

# test: tests/parsers/vasp/test_mapping.py
    @pytest.mark.parametrize(
        ('parameters', 'expected'),
        [
            ({'GGA': 'PE'}, 'PBE'),
            ({'LHFCALC': True, 'HFSCREEN': 0.2}, 'HSE06'),
            ({'LHFCALC': True, 'AEXX': 0.25, 'HFSCREEN': 0.5}, None),
        ],
    )
    def test_maps_xc_tags_to_canonical_functionals(self, parameters, expected):
        assert functional_key_from_params(parameters) == expected
````

   The absent case, ORCA. The transformer is the function named in the mapping annotation of the `DFT` section, `('get_dft_methods', ['.@'])`, and it runs on every ORCA output because the SCF-settings block it reads is written for every calculation. Without a functional in that block it must return nothing, otherwise a Hartree-Fock run would get a DFT method section:

````python
# implementation: parsers/orca/parser.py
    def get_dft(self, scf_settings: dict[str, Any]) -> dict[str, Any]:
        xc = self._build_xc(scf_settings)
        if not xc:
            return {}

# test: tests/parsers/orca/test_mapping.py
    def test_absent_xc_yields_no_method(self):
        assert OutParser().get_dft({}) == {}
````

   LAMMPS and phonopy have no `test_mapping.py` yet; a change to their mapping adds the file.

3. **Failure visibility.** Give the transformer a value it cannot map, for example a string where a number is expected, and check that something visible happens: an exception, or a log entry at error level ([assertion strength](#axes)). If the quantity just comes out empty, or the bad value is passed on unchanged, and nothing is reported, the failure is silent, and that is what this test exists to catch. No parser has such a test yet. The ORCA exact-exchange fraction shows why the item exists: the implementation passes whatever the reader captured straight through, so a non-numeric value surfaces only later, at schema assignment, if at all. The test states what should happen instead:

````python
# implementation: parsers/orca/parser.py, _build_xc
        hf_fraction = self._to_scalar(scf_settings.get('fraction_hf_exchange'))
        if hf_fraction is not None:
            xc['global_exact_exchange'] = hf_fraction

# test (proposed): tests/parsers/orca/test_mapping.py
    def test_non_numeric_exchange_fraction_is_rejected(self):
        source = {'exchange_functional': 'B88', 'fraction_hf_exchange': 'abc'}
        with pytest.raises(ValueError):
            OutParser().get_dft(source)
````

   At the entry point the same idea uses `monkeypatch` to replace the reader's output with a malformed subtree and `caplog` to check that `parse` reports it at error level rather than producing a thinner archive silently.

4. **Integration.** The integration test runs the whole parser on a real output file and checks a few values in the resulting archive. Add the new quantity there only if it identifies the calculation, meaning a wrong file, a wrong unit or a wrong ordering would change it; if it does not, the mapping test of item 2 is the right place ([layer 4](how_to/test_parsers.md#layer-4-parser-integration), [scope](#scope-different-test-types-serve-different-purposes-and-catch-different-categories-of-bugs)). Most integration classes check only lengths and shapes so far; where a class has no identifying values yet, this item adds the first ones. VASP, where the whole entry point runs on one fixture and the chemical symbols and the SCF step count identify it:

````python
# implementation: parsers/vasp/parser.py
class VASPParser(MatchingParser):
    def parse(self, mainfile: str, archive: 'EntryArchive', logger: 'BoundLogger', ...) -> None:

# test: tests/parsers/vasp/test_integration.py
class TestAgacRelaxOutcarArchive(VASPIntegrationSuite):
    archive_fixture = 'agac_relax_outcar_archive'  # VASPParser().parse on OUTCAR

    @pytest.mark.integration
    def test_archive_contract(self, archive):
        symbols = [s.chemical_symbol for s in archive.data.model_system[0].particle_states]
        assert symbols == ['Ac', 'Ag']
        assert len(archive.data.outputs[0].scf_steps.energies_total) == 11
````

**Helper functions**

5. **Shared helpers.** A function that is added to `parsers/utils`, or moved there because a second parser needs it, gets its own tests in `tests/parsers/test_general_utils.py` or a file next to it: a `parametrize` table of inputs and expected outputs, written out as literals, with the empty input as one row ([unit tests](#unit-tests)). A helper that only one parser uses stays in that parser and is covered by its mapping tests (item 2). The mapping-table helper, with the expected table written out rather than derived:

````python
# implementation: parsers/utils/general.py
def create_mapping_table(file_parser: Any, archive_parser: Any, ...) -> list[dict[str, Any]]:

# test: tests/parsers/test_general_utils.py
def test_create_mapping_table_compares_file_quantities_with_mapper_paths():
    file_parser = _FileParser([_Quantity('n_states'), _Quantity('unused')])
    archive_parser = _FileParser([])
    archive_parser.mapper = _Mapper(['n_states'])

    assert create_mapping_table(file_parser, archive_parser) == [
        {'quantity': 'n_states', 'mapped': True, 'mapping': ['n_states']},
        {'quantity': 'unused', 'mapped': False, 'mapping': []},
    ]
````

**Before the pull request**

6. **Markers.** Mark reader, mapping and helper tests `unit`, the integration class `integration`, and only the shared normalization suite `pipeline`; add `large_fixture` to a class whose input file lives in the fixture bundle. The markers decide in which CI stage a test runs ([markers and CI](how_to/test_parsers.md#pytest-markers-and-ci), [selective execution](#performance-optimization)). The fast gate is `pytest -m 'unit and not large_fixture'`.

````toml
# definition: pyproject.toml
markers = [
    "unit: recognition, source-reader, utility, and mapping-contract tests",
    "integration: tests of the public parser entry point",
    "pipeline: parser and NOMAD normalization compatibility tests",
    "large_fixture: test requires a fixture excluded from the default CI checkout",
]
````

````python
# use: tests/parsers/orca/test_integration.py
@pytest.mark.large_fixture
class TestCASCIQuantumDots(OrcaParserIntegrationSuite):
    # ~19 MB mainfile; `large_fixture` moves this class to the nightly gate.
    archive_fixture = 'casci_qd_archive'
````

7. **Fixture.** Prefer a few lines of made-up output text, as in item 1. If a real output file is added, keep it small, note where it came from, which program version wrote it and under which licence it may be redistributed (in the fixture's `conftest.py` entry or a `README` in its `tests/data` directory), and put it in the fixture bundle when it is large ([fixtures and manifests](how_to/test_parsers.md#fixtures-and-expected-data-manifests), [synthetic minimal data](#specific-testing-patterns)). ORCA, where small files stay in `tests/data` and the bundle root comes from the environment:

````python
# fixtures: tests/parsers/orca/conftest.py
DATA_DIR = Path(__file__).resolve().parents[2] / 'data' / 'orca'
LARGE_FIXTURE_DIR = (
    Path(os.environ['NOMAD_SIM_PARSERS_LARGE_FIXTURE_ROOT']) / 'orca'
    if 'NOMAD_SIM_PARSERS_LARGE_FIXTURE_ROOT' in os.environ
    else DATA_DIR
)

@pytest.fixture(scope='module')
def ri_mp2_water_archive() -> EntryArchive:
    return parse_orca('RI_MP2_water.out')  # small, in tests/data

@pytest.fixture(scope='module')
def casscf_state_averaged_archive() -> EntryArchive:
    return parse_large_orca('CoPc_CASSCF_SA.out')  # in the bundle
````

8. **Output check.** Run the additive-output check against `develop`. It parses the fixtures with both revisions and compares the archives value by value: new values are expected, every changed or removed value has to be explained in the pull request ([additive-output check](how_to/test_parsers.md#additive-output-check), [the role of snapshots](#the-role-of-snapshots)).

````bash
tests/parsers/check_additivity.sh develop HEAD orca
````

9. **Matrix.** Add the quantity to the parser's feature matrix, the table that lists per supported feature which reader, mapping and integration test covers it, in the parser's tests or documentation; no parser has started this table yet, so the first change creates it ([test ownership and minimum coverage](how_to/test_parsers.md#test-ownership-and-minimum-coverage)).

### Schema: adding or changing normalization logic

1. **Is a test due?** Yes unless the logic is trivial. The test is worth writing when it states the behaviour in terms simpler than the code: a table of input and expected output, a law the result obeys. A test that only re-reads a default, an IRI or a type (a guard) does not meet that bar and is not added ([assertion strength](#axes), [testing the framework instead of your code](#common-testing-anti-patterns-to-avoid)).
2. **Cases.** Write the cases as rows of one `pytest.mark.parametrize` table with an `ids` entry per row, and write each expected value out as a literal, never as a formula that repeats the rule under test. Include rows for the empty input, the unset quantity and the undetermined result ([parameterized tests](#specific-testing-patterns), [oracle](#axes)).
3. **Negative.** For every input that must produce a warning or an exception, one test that checks the captured log or uses `pytest.raises` ([unit tests, when to write](#unit-tests)).
4. **Visibility.** Where the code catches exceptions, one test that forces a failure inside the protected block and checks that it is logged at error level or raised again ([assertion strength](#axes)).
5. **Docs.** Regenerate the schema documentation and commit it with the change.
6. **Recurring constraints.** When the same kind of constraint (a shape, a length, a bounded range) needs a test for the third time, do not write the third test; propose a mechanism that enforces it everywhere, in the metainfo definition, in the bounded types, or as a contract ([maintainability](#key-principles-summary)).

## Fundamentals

> This document is conceptual background: why the different kinds of parser test exist and what each is good at. The concrete, NOMAD-specific protocol that puts these ideas into practice — layer layout, fixtures, markers, and CI stages — lives in [`how_to/test_parsers.md`](how_to/test_parsers.md). Prefer that document when actually writing tests; this one explains the reasoning behind it.

Tests vary along several independent dimensions. This document develops one of them in depth — the **scope** of code each test exercises, running from a single function up to the complete parser pipeline. This axis is the most immediately useful to us today: we currently write almost exclusively unit tests, so laying out the full scope ladder makes the missing higher tiers explicit and helps surface parser behaviour that is currently going untested. The protocol refines the generic unit/integration/E2E pyramid into five parser-specific layers — recognition, source extraction, mapping contract, parser integration, and NOMAD-pipeline compatibility — so that a failure points at the specific contract that broke.

Scope is not the only axis. A second is how inputs are chosen: the examples below are all *example-based* (hand-picked inputs with known outputs). Asserting invariants over inputs rather than fixed cases is already part of the protocol as *metamorphic* tests (for example, changing whitespace or the source unit rescales the result predictably); generalising these with *property-based* testing (e.g. Hypothesis, over generated inputs) is a natural next step. Treat the scope ladder below as one well-developed dimension of the strategy, not the whole of it.

### Axes

Three independent choices describe any test, and the policy names all three so that a missing kind is visible rather than forgotten.

* **Scope**: how much code runs, from one function to the complete NOMAD pipeline. This is the ladder developed below, refined into the five parser layers of the protocol.
* **Input selection**: how the inputs are chosen. *Example-based* tests use hand-picked inputs with known outputs; *metamorphic* tests transform an input and check the relation between the two outputs; *property-based* tests generate inputs and check a law that must hold for all of them.
* **Oracle**: what decides pass or fail. An *explicit value* stated in the test or in a reviewed manifest (`case.yaml`); a *relation* between runs (metamorphic and property-based tests need no expected value); a *differential* comparison against another implementation or another revision of the same code; or a *snapshot*, a stored previous output of the code under test.

The last two oracles are the ones that most easily pass when they should fail, and the next section fixes how we use them.

Two further labels are worth recording per test, because they answer questions the three axes do not. The first is the *strength of the assertion*: a **guard** asserts only that something exists, is `None` or not, or has a type, and so passes on any refactor that keeps the attribute names; a **computed-value** test compares a value the code derived; a **round-trip** test reads back what it set; a **negative** test expects an exception or a logged failure. Guards are cheap and sometimes the right check (a recognition test returns a boolean; "the normalizer created no section" is a `None`), but a suite made of guards detects renames, not wrong values, and mutation testing on the changed lines is the way to tell the two apart. The second label is *whose behaviour is under test*: a test whose assertions only exercise the metainfo (defaults, `m_validate`, `m_from_dict(m_to_dict())`) with the parser's sections as carriers belongs in nomad-lab, and one such smoke check per package is enough. Every test in this repository can be placed on all of these axes; a test that cannot is usually one that asserts nothing.

### The role of snapshots

A snapshot records what the code produced at some point and fails when the output changes. It says nothing about whether either output was correct, so a snapshot is never the oracle for a parsing rule, a unit conversion, or a mapped value; those are checked at the lowest layer that can reproduce them, with the expected value stated in the test. Committed snapshot files of parser archives are not part of this repository.

What this repository does use is the *differential* form at end-to-end scope: the additive-output check (`tests/parsers/check_additivity.sh`, documented in [`how_to/test_parsers.md`](how_to/test_parsers.md)) runs two revisions of the parsers on the same fixtures, compares `archive.data` and `workflow2` leaf by leaf, passes on additions, warns on growth and fails on any removed or changed leaf. It is transient (the snapshots are regenerated on each run and never committed), type-sensitive, and run on manual dispatch and before a rollout, because deciding whether a change is intended is a human decision. The same differential idea, with the legacy parsers as the reference, is the parity report in the Area C distribution. Both are the end-to-end tier of the strategy: they catch what no one anticipated, and a disagreement they report is a finding to triage, not a verdict.

The `case.yaml` manifests of the protocol are the opposite of snapshots: they hold a few reviewed expectations per fixture, written by a person from the source file, and a change to a manifest is reviewed as a behaviour change. When a pull request changes parser output, the sequence is therefore: add the lowest-layer test with the new expected value, update the manifest if a public value changed, and run the additive check to see what else moved.

### Scope: different test types serve different purposes and catch different categories of bugs

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

**Example:** Testing a value-parsing function we own

````python
def test_parse_key_value_line():
    """The parser splits a raw output line into a key and a typed value."""
    key, value = parse_key_value("total energy = -1.5")

    assert key == "total energy"
    assert value == -1.5  # coerced to float, not left as the string "-1.5"
````

The assertions target the behaviour our own function implements: extracting the key and coercing the value to a float. They deliberately avoid re-asserting properties that belong to a library — for example that a `pint` quantity carries `ureg.meter`, or that `re.match` splits a string. Asserting a framework's own guarantees tests the framework, not our code; it is a common anti-pattern, listed at the end of this document.

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

The marker names above are illustrative. This repository's protocol defines a fixed set — `unit`, `integration`, `pipeline`, and `large_fixture` — mapped to the CI stages described in [`how_to/test_parsers.md`](how_to/test_parsers.md). Use those, not the generic names shown here.

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

5. **Testing the framework instead of your code**

````python
   # BAD: asserts a guarantee owned by the library, not our logic
   def test_unit_conversion():
      q = 100.0 * ureg.meter
      assert q.units == ureg.meter  # this is pint's behaviour, not ours

   # GOOD: asserts what our code computes from a raw input
   def test_energy_is_parsed_as_float():
      key, value = parse_key_value("total energy = -1.5")
      assert value == -1.5
````
