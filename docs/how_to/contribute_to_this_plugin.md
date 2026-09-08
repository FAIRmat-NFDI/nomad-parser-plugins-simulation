# Contribute to This Plugin

Before opening a pull request, run the relevant tests and include any generated
files that are affected by your changes.

## Update the parser mapping report

The parser mapping report is a committed snapshot of the file-parser quantities
and their archive mappings. Its generation is part of CI. Whenever a pull
request contains source-code changes to be merged into `develop`, regenerate
the report and commit the resulting update as part of that pull request:

```sh
uv run nomad-sim-parser mapping-report \
  --override docs/reference/parser_mapping_report_overrides.yaml
```

This updates
`docs/reference/parser_mapping_report.md`. Review the generated diff, including
any intentionally unmapped quantities, before submitting the pull request.

## Parser authoring conventions

These conventions come from recurring review feedback; follow them when writing
or changing a parser.

**Prefer stateless transformers over parser state.** Do not gate behaviour with
mutable flags on the parser instance — they have to be reset per parse and leak
across reused parser instances. Instead resolve the context a transformer needs
(for example a frame index) through the mapping annotation and pass it as an
explicit argument. For per-frame identity, the frame builder stamps a
`frame_index` on each frame and the identity transformer emits only when it is
`0`.

**Keep transform logic in the transformer.** Derive values in the annotated
transformer functions driven by the mapping, not by imperatively
post-processing `archive.data` in the `ArchiveWriter`.

**Store per-particle identity once.** In a multi-frame `model_system` sequence —
an MD trajectory or a geometry-optimization step series — the per-particle
identity (`particle_states`) is frame-independent; attach it to the first
(topology) frame only, so it is not duplicated per frame. Duplicating it scales
the archive and the Elasticsearch index document with `n_frames × n_particles`
and can push a single entry past the Elasticsearch payload limit, failing the
whole upload (see FAIRmat-NFDI/nomad-simulations#474). The convention is enforced
by `tests/parsers/common.py::assert_identity_populated_once`; wire it into a
parser's multi-frame test. Deviate only when identity genuinely changes between
frames (reactive or alchemical simulations).

**Consolidate test helpers.** Put shared test assertions in the existing suite
module (`tests/parsers/common.py`) rather than adding standalone helper files.

**Make varying behaviour configurable.** For behaviour that legitimately varies,
prefer a configuration setting consistent with the existing ones (such as the
trajectory-sampling rate) over a hard-coded policy.
