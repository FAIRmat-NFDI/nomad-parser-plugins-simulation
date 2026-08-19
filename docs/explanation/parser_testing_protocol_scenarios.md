# Draft: scenario-contract parser testing protocol

## Status

This is an alternative to the layer-oriented parser testing protocol. It is a draft
for evaluation and is not yet a contributor requirement. An executable ABINIT pilot
is available under `tests/scenarios/abinit`; it runs alongside, and does not replace,
the layer-oriented tests under `tests/parsers/abinit`.

## Design goal

Describe each supported scientific scenario once and use that description to test
all observable parser boundaries. The unit of coverage is a **capability**, such as
"spin-polarized DOS" or "variable-cell optimization", rather than a parser function
or an archive quantity.

The protocol should answer four questions for every capability:

1. Did NOMAD recognize the calculation and select the correct files?
2. Did the source reader recover the facts present in those files?
3. Did mapping preserve their meaning, units, shape, and relationships?
4. Does the public parser expose the promised user-facing result?

## Scenario contracts

Every parser owns a collection of scenario directories:

```text
tests/scenarios/abinit/
  spin-polarized-scf/
    scenario.yaml
    output.out
    calculation_o_DS2_DOS
  truncated-relaxation/
    scenario.yaml
    output.out
```

Each `scenario.yaml` is a manually reviewed contract:

```yaml
schema_version: 1
id: abinit-spin-polarized-scf
capabilities:
  - structure.periodic
  - method.dft.xc
  - electronic.scf
  - electronic.dos.spin_polarized

input:
  mainfile: output.out
  auxiliary_files:
    - calculation_o_DS2_DOS

recognition:
  matches: true
  children: []

source_facts:
  program_version: "9.10.4"
  nsppol: 2
  scf_energies:
    value: [-10.2, -10.5]
    unit: hartree
  dos:
    channels: 2
    points_per_channel: 3

archive_facts:
  data.program.name: ABINIT
  data.model_system[0].positions:
    shape: [2, 3]
    unit: angstrom
  data.outputs[0].electronic_dos:
    length: 2
  data.outputs[0].electronic_dos[0].value:
    shape: [3]
    unit: 1 / hartree

invariants:
  - particle_count_matches_positions
  - dos_axes_match_values
  - output_references_resolve

diagnostics:
  errors: []
  warnings:
    allow: []
```

Contracts contain selected scientific facts, not serialized archives. Unknown keys
must fail validation so that misspelled expectations cannot silently pass.

## Scenario runner

A shared pytest runner executes every scenario through explicit checkpoints:

```text
scenario
   |-- recognize(files) --------> recognition result
   |-- extract(files) ----------> source facts
   |-- map(source facts) -------> archive fragment
   `-- parse(files) ------------> public archive + diagnostics
```

Each parser supplies a small adapter:

```python
class AbinitScenarioAdapter(ParserScenarioAdapter):
    def recognize(self, case): ...
    def extract(self, case): ...
    def map(self, source): ...
    def parse(self, case): ...
```

The adapter exposes existing parser boundaries; it must not reimplement parsing.
If a parser cannot implement `extract` or `map` independently, that is an
architectural limitation to resolve before adding more scenarios.

The common runner owns generic checks for values, approximate values, units, array
shapes, lengths, references, and diagnostics. Parser-specific scientific invariants
are registered functions, not expressions evaluated from YAML.

## Two fixture forms

Scenarios use one of two fixture forms.

### Minimal fixtures

Minimal fixtures isolate one syntax or capability and are the default. They should
be human-readable and normally smaller than 20 KiB. They run at all four
checkpoints.

### Representative fixtures

Representative fixtures are real calculations used to verify interactions between
capabilities. They run recognition and public parsing, while selected sections may
also run extraction. Large binary inputs live in a versioned external fixture
bundle identified by URL, checksum, and license metadata.

A representative fixture must not be the only evidence for a parsing rule that can
be expressed with a minimal fixture.

## Independent expected values

Expected values must come from one of these sources, recorded in the scenario:

- a value visible in the source fixture;
- a documented analytic result;
- an independently maintained reference implementation; or
- a domain expert's reviewed calculation.

Generating a contract from the parser under test is prohibited. A developer tool
may print candidate values, but generated values enter the contract only after
manual comparison with an independent source.

Every numeric expectation specifies one of:

- exact comparison for identifiers, counts, flags, and discrete indices;
- absolute tolerance for values near zero;
- relative tolerance for measured or calculated values; or
- shape/sample comparison for large arrays.

Defaults are not applied silently. Omitting a tolerance means exact comparison.

## Capability registry

A repository-level registry defines stable capability names and required
invariants:

```yaml
electronic.dos.spin_polarized:
  requires:
    - electronic_dos
    - energy_axis
    - spin_channel_ordering
  invariants:
    - dos_axes_match_values
```

Each parser publishes the capabilities it supports. CI compares the declaration
with scenario coverage and fails if a supported capability has no scenario.

This replaces informal claims such as "the parser supports DOS" with an executable
coverage statement. A capability can be marked `experimental` with an issue and
expiry date, but it cannot be silently uncovered.

## Negative scenarios

Invalid and incomplete inputs use the same contract format. Their expected outcome
must be one of:

- `reject`: the file is not recognized as a mainfile;
- `fail`: parsing raises a named domain exception;
- `partial`: parsing succeeds with listed missing sections and diagnostics; or
- `ignore`: an invalid optional auxiliary file is ignored with a diagnostic.

At minimum, each parser provides scenarios for a truncated mainfile, a missing
required auxiliary file, an unsupported version marker, and conflicting auxiliary
files.

Broad `raises(Exception)` checks and acceptance of arbitrary warnings are not
allowed.

## Transform tests

The runner may derive additional cases from a scenario using declared
transformations:

```yaml
transforms:
  - gzip_mainfile:
      preserves: [source_facts, archive_facts]
  - crlf_line_endings:
      preserves: [source_facts, archive_facts]
  - reorder_auxiliary_files:
      preserves: [archive_facts]
```

Only transformations with a clear domain-preserving meaning are allowed. Derived
cases reuse the original oracle and do not add generated expected values.

## Diagnostics are part of the contract

Silent partial parsing is a parser behavior and must be testable. The runner captures
structured logs and compares their stable event names and severity. Contracts do
not match complete rendered messages, timestamps, or stack traces.

New diagnostics should use stable event identifiers such as
`abinit.missing_dos_file`. Until structured identifiers are available, tests may
match a short stable message fragment.

## CI policy

CI groups scenarios instead of individual test functions:

1. **Fast:** all minimal scenarios at recognition, extraction, and mapping
   checkpoints.
2. **Core:** minimal and repository-sized representative scenarios through the
   public parser.
3. **Compatibility:** representative scenarios followed by schema validation and
   NOMAD normalization.
4. **Extended:** external large fixtures, supported Python versions, and dependency
   version boundaries.

A pull request runs scenarios for changed parsers plus the shared harness tests.
The full core matrix runs before merge. Compatibility and extended matrices run on
a schedule and before release.

Failed selected fixtures are errors; tests never skip because a fixture is missing.
External-fixture jobs verify all checksums before collection.

## Change rules

- A parser bug fix adds or updates the smallest scenario that reproduces it.
- A new capability requires registry metadata and at least one positive scenario.
- A new source-format version adds a scenario even when behavior is unchanged.
- An intentional archive-contract change updates affected scenarios in the same
  pull request and explains the compatibility impact.
- Updating many unrelated expected values requires explicit reviewer approval; it
  is treated as a contract migration, not snapshot maintenance.

## Migration

1. Build and validate the scenario schema and common assertion engine.
2. Implement adapters for ABINIT and one structurally different parser such as
   VASP or GROMACS.
3. Convert existing fixtures into a capability inventory without deleting tests.
4. Create minimal scenarios for the most failure-prone extraction and mapping rules.
5. Run old and scenario suites together until capability coverage is equivalent.
6. Remove redundant assertion-heavy end-to-end tests.
7. Make capability coverage and diagnostic checks required in CI.

## Evaluation criteria

Pilot the protocol on two parsers and accept it only if:

- one scenario can expose failures independently at extraction and mapping;
- a developer can add a simple regression scenario without editing shared code;
- scenario errors name the capability, checkpoint, and failed path;
- contract review does not resemble reviewing full archive snapshots;
- fast scenarios complete within the agreed pull-request budget; and
- capability coverage reports reveal meaningful gaps missed by line coverage.

## Trade-offs

This design reduces duplicated expectations and makes supported capabilities
visible, but it introduces a scenario schema, adapter interface, and shared runner
that must remain stable. It is best suited to a repository with many parsers that
share scientific concepts. For a small parser or an early migration, explicit
layer-specific Python tests are simpler and easier to debug.
