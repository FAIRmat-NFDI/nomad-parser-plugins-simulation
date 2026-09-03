# Parser testing protocol

## Purpose

Parser tests must say which contract failed. A test that only runs a complete parser
and inspects the resulting `EntryArchive` mixes several contracts:

1. mainfile and auxiliary-file recognition;
2. extraction from the source format;
3. conversion and mapping to NOMAD metainfo;
4. archive assembly and schema validation; and
5. normalization.

Such a test is useful as a smoke test, but it is too broad to be the primary way of
testing a parser. Failures are hard to locate, malformed-input behavior is difficult
to exercise, and small extraction or mapping rules require large fixtures.

This protocol makes each boundary independently testable while retaining a small
number of complete-parser tests.

## Test architecture

Each parser is tested as the following pipeline:

```text
files -> recognition -> source reader -> parsed source data -> mapper -> archive
                                                                  |
                                                                  v
                                                             normalizer
```

The objects at the two internal boundaries are test interfaces:

- **parsed source data** is the dictionary, array, or typed object returned by a
  file parser such as a `TextParser`, `FileParser`, or format-specific reader;
- **archive fragment** is the metainfo section written by an `ArchiveWriter` or
  equivalent mapping component before normalization.

New parser code must keep extraction and mapping callable independently. Existing
parsers should expose these seams when they are migrated; a rewrite to one common
implementation is not required.

### Layer 1: recognition

Test mainfile matching, auxiliary-file recognition, compression handling, and file
precedence without parsing a calculation. Use a temporary directory containing
empty or minimally named files.

Assertions should cover:

- matching and non-matching filenames or headers;
- deterministic selection when several candidates exist;
- missing, duplicate, and compressed auxiliary files; and
- paths independent of the current working directory.

These are unit tests and must not create an `EntryArchive`.

### Layer 2: source extraction

Call the lowest-level reader with a minimal source fragment and assert its parsed
source data. One fixture should demonstrate one syntax feature or source-format
variant. Prefer inline text or a fixture below 10 KiB; use a binary fixture only
when the format requires it.

Check semantic values, not merely presence or collection sizes. For example, test
the parsed SCF energies, spin channel, atom labels, and source unit. Include cases
for repeated sections, optional fields, truncated files, malformed numeric values,
and supported program-version variants.

The expected result belongs in the test when it is short. Larger results use a
reviewed YAML manifest containing only the fields relevant to that feature. Do not
store a serialized parser object or a complete archive as the oracle.

### Layer 3: mapping contract

Feed controlled parsed source data directly to the mapping/archive-writing layer.
This isolates source-to-metainfo semantics from regular expressions and third-party
file readers.

Every supported archive concept must have a mapping-contract test that checks:

- target section type and path;
- value and explicit unit conversion;
- array shape and index or spin ordering;
- references between systems, methods, outputs, and workflow tasks;
- distinction between absent data, an empty collection, zero, and `False`; and
- merge or precedence rules when multiple files provide the same concept.

Use `pytest.approx` or `numpy.testing` with a tolerance chosen for the source
precision. Convert quantities to a named unit before comparison. Avoid conditional
assertions such as `if output.electronic_dos`; if a fixture promises DOS, its
absence must fail the test.

Mapping tests create only the smallest required archive or section. They do not run
normalizers.

### Layer 4: parser integration

Run the public parser entry point on one small, representative fixture for each
calculation mode the parser supports, for example single point, geometry
optimization, molecular dynamics, or a post-processing output.

An integration test verifies the cross-component invariants rather than repeating
every extracted number:

- parsing completes without unexpected error-level logs;
- the archive can be serialized and deserialized;
- metainfo validation succeeds;
- systems, methods, outputs, and workflow references resolve;
- their step counts and array dimensions agree; and
- a short list of scientifically identifying values is correct.

The identifying values should catch a wrong file, wrong calculation, wrong unit,
or wrong ordering. They are not a snapshot of the archive.

### Layer 5: NOMAD pipeline compatibility

Keep normalization outside parser tests. A separate compatibility test runs the
parser and the required NOMAD normalizers for one fixture per parser, then checks a
small public consumer contract such as the representative system, material, and
results properties.

This layer detects integration drift between packages. It must not be used to prove
that an individual regular expression or mapping is correct. If it fails, the lower
layers should identify whether the parser contract itself is still satisfied.

### Shared simulation parser contract

Do not copy the parser-independent parts of layers 4 and 5 into every parser.
Compose the Simulation, model-system, and workflow suites with one representative,
module-scoped archive fixture and the expected program name:

```python
from tests.parsers.common import (
    SimulationParserTestSuite,
    WorkflowTestSuite,
)


class TestExampleParserSuite(
    SimulationParserTestSuite, WorkflowTestSuite
):
    archive_fixture = 'example_archive'
    expected_program_name = 'ExampleCode'
```

Use this composition in `test_integration.py`. `SimulationParserTestSuite` checks
validation, core Simulation sections, and model-system structural integrity, while
`WorkflowTestSuite` checks workflow presence and round-trip serialization. A parser
may omit a suite only when that archive concept is not part of its public contract.
Subclass `SimulationParserPipelineTestSuite` in
`test_pipeline.py` with the same fixture and program name to check normalization
compatibility. Keep recognition rules, source-reader behavior, mapping details,
supported optional sections, and scientific reference values in parser-local
tests.

## Fixtures and expected-data manifests

Organize tests by parser and layer:

```text
tests/parsers/<code>/
  test_recognition.py
  test_reader.py
  test_mapping.py
  test_integration.py
  cases/
    scf-minimal/
      case.yaml
      <source files>
```

`case.yaml` records why a fixture exists and its stable expectations:

```yaml
id: silicon-scf
purpose: spin-unpolarized SCF with two ionic steps
program_version: "6.8"
mainfile: output.out
features: [structure, scf, forces]
expect:
  program.name: ExampleCode
  model_system[0].positions:
    shape: [2, 3]
    unit: angstrom
    sample:
      "[1, 2]": 1.357
  outputs:
    length: 2
  outputs[1].total_energies[0].value:
    value: -10.42
    unit: eV
    rel: 1.0e-8
```

The manifest is a declarative aid, not a generic snapshot framework. A shared
assertion helper may implement `value`, `unit`, `shape`, `length`, and selected
`sample` checks. It must reject unknown keys and produce an error containing the
case id and archive path. References and domain invariants remain explicit Python
assertions.

Fixture rules:

- include provenance and the feature under test in `case.yaml`;
- remove unrelated output and redact user or machine paths;
- prefer the smallest source that the real reader accepts;
- never generate expected data by running the parser under test in CI;
- review manifest changes as behavior changes, not routine snapshot updates; and
- mark large or externally supplied fixtures explicitly rather than silently
  skipping when they are absent.

## Negative and metamorphic tests

Each parser must have negative tests for a truncated mainfile, a missing required
side file, and invalid content. The expected behavior—exception, warning with
partial archive, or ignored optional file—must be explicit.

Where applicable, add transformations whose result is known without another
oracle:

- changing whitespace or line endings preserves parsed values;
- parsing a compressed and uncompressed file produces equivalent source data;
- translating every atom preserves energies and shifts positions only;
- reordering auxiliary files does not change deterministic precedence; and
- changing the source unit rescales the mapped value correctly.

These tests cover entire classes of inputs more effectively than adding more
archive quantity counts.

## Test ownership and minimum coverage

For every advertised parser feature, maintain this matrix in the parser's tests or
documentation:

| Feature | Reader case | Mapping case | Integration case |
| --- | --- | --- | --- |
| Structure | required | required | representative |
| Method and numerical settings | required | required | representative |
| Energies and SCF history | required | required | representative |
| Forces and stress | when supported | when supported | one identifying value |
| DOS, bands, or eigenvalues | when supported | when supported | one shape/value |
| Workflow | source inputs if parsed | required | required by calculation mode |

A pull request that adds or fixes a parsing rule must add the lowest-layer
regression test capable of reproducing the bug. Add or change an integration
fixture only when the public parser behavior or a cross-component invariant changes.

## Pytest markers and CI

Use these markers:

- `unit`: recognition, reader, utility, and mapping-contract tests;
- `integration`: public parser entry-point tests;
- `pipeline`: parser plus NOMAD normalization;
- `large_fixture`: tests whose fixture set is too large for the default checkout.

CI runs in three stages:

1. **Pull-request fast gate:** all `unit` tests and integration tests affected by
   changed parser directories. Target: less than two minutes.
2. **Pull-request package gate:** all non-large `unit` and `integration` tests.
3. **Nightly or release gate:** `pipeline` and `large_fixture` tests using a
   versioned fixture bundle with checksums.

These stages are implemented in `.github/workflows/actions.yml`. The nightly
fixture manifest is `tests/fixtures/nightly.sha256` and is verified by
`tools/verify_fixture_bundle.py`; update its checksum in the same change as a
fixture update. The fast gate runs all unit tests, then selects integration
directories changed by the pull request. Changes to shared parser code fall back
to all integration tests.

Tests must not pass merely because a fixture is unavailable. A required fixture is
a repository or CI setup error. Optional large fixtures may be deselected by marker,
but a selected test with a missing fixture must fail.

## Additive-output check

Refactors are expected to **preserve parser output**. The `additive-output-changes`
workflow, and the `tests/parsers/check_additivity.sh` script it wraps, verify
this by snapshotting a `target` baseline ref and a `source` ref under test and
comparing their `archive.data` *leaf by leaf*. Each ref is snapshotted by running
*its own* parser code via `PYTHONPATH`, so the caller's environment is *never*
mutated. Snapshots are plain JSON written to `target_snapshots.json` and
`source_snapshots.json` at the repository root (override with the `TARGET_JSON` /
`SOURCE_JSON` environment variables), generated from scratch on each run and
**never committed** -- both paths are in `.gitignore`, so the check is a transient
comparison rather than a stored baseline.

Each parser falls into one of **three tiers**, in increasing severity:

1. **`identical`** -- no leaf added, removed, or changed: passes *silently*.
2. **`additive`** -- only additions: passes, but is surfaced as a *non-failing*
   warning annotation under GitHub Actions, so growth stays visible.
3. **removed or changed** -- a *hard failure* (exit 1), printing `old -> new`.

So additions *never* block, while any removal or value change *always* does -- an
intentional path move or value change is *meant* to fail here, and is reviewed by
running the script locally and confirming the diff. The comparison is made hard to
fool so that an `identical` verdict is trustworthy: leaf comparison is
**type-sensitive** (a value changing type is flagged, not just its value), and
every container emits a **`[type]`** leaf and each list a **`[len]`** leaf, so an
*empty or absent* container is not invisible.

The check is **manual only** (`workflow_dispatch`), *never* part of the
pull-request gates, because judging whether a change is intended is a **human
decision**. The same review holds in CI: the full comparison report is *printed
to the job log* and, together with both snapshots, **uploaded as a build
artifact** (even on a failing run), so the leaf-level diff can be inspected
without re-running locally. Reused state is guarded against staleness on two
fronts: the CI **cache key** includes both the target and source SHAs (the
autodetected parser set depends on their diff) *and* the generator's hash, so a
changed environment or generator *cannot* restore a mismatched snapshot; and each
snapshot carries a **provenance block** (interpreter, serializer, harness,
`uv.lock`, and fixture hashes), reported for context but *never diffed*, since two
refs may legitimately differ there.
