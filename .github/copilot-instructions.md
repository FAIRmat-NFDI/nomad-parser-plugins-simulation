# Copilot instructions

## Testing

The authoritative, NOMAD-specific protocol is `docs/how_to/test_parsers.md`; follow it when
writing or reviewing parser tests. `docs/parser_testing_strategy.md` gives the background
reasoning. The rules below summarise the protocol.

Test each parser as a pipeline and test every contract boundary independently, rather than
only running the whole parser and inspecting the resulting `EntryArchive`. The layers are:

- **recognition** — mainfile/auxiliary-file matching, precedence, and compression handling; no `EntryArchive` is created.
- **source extraction** — the lowest-level reader on a minimal source fragment; assert semantic values (energies, spin channel, labels, source unit), not just presence or counts.
- **mapping contract** — feed controlled parsed source data to the mapping/archive-writing layer; assert target section path, value and explicit unit, array shape and ordering, references, and the distinction between absent, empty, zero, and `False`. Do not run normalizers.
- **parser integration** — the public entry point on one small fixture per calculation mode; assert cross-component invariants and a few scientifically identifying values, not a full archive snapshot.
- **pipeline compatibility** — parser plus the required NOMAD normalizers for one fixture, checking a small public consumer contract only.

Conventions:

- Use only the markers `unit` (recognition, reader, mapping-contract), `integration` (public entry point), `pipeline` (parser + normalization), and `large_fixture`. Do not invent others.
- Lay tests out under `tests/parsers/<code>/` as `test_recognition.py`, `test_reader.py`, `test_mapping.py`, and `test_integration.py`, with fixtures under `cases/<name>/case.yaml`.
- Compose the shared suites (`SimulationParserTestSuite`, `WorkflowTestSuite`, `SimulationParserPipelineTestSuite`) instead of copying the parser-independent checks into each parser.
- Treat `case.yaml` as a declarative aid, not a generic snapshot; keep references and domain invariants as explicit Python assertions, and review manifest changes as behaviour changes.
- Test the behaviour our own code implements, never a guarantee owned by a library (do not assert that `pint` returns `ureg.meter` or that `re.match` splits a string).
- Add negative tests (truncated mainfile, missing required side file, invalid content) and metamorphic tests (whitespace, compression, atom translation, auxiliary-file reordering, unit rescaling) where applicable.
- A pull request that adds or fixes a parsing rule must add the lowest-layer regression test that reproduces the bug.
- Refactors must preserve parser output: the additive-output check fails on any removed or changed archive leaf and only warns on additions.
