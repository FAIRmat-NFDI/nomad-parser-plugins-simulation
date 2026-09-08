# Copilot / AI agent instructions

Follow the parser authoring conventions in
`docs/how_to/contribute_to_this_plugin.md` ("Parser authoring conventions"). In
particular:

- Prefer stateless transformers; pass the context a transformer needs (for
  example a frame index) as an explicit mapping-annotation argument rather than
  mutable parser state.
- Keep transform logic in the transformer functions, not in post-conversion
  loops in the `ArchiveWriter`.
- Store per-particle identity (`particle_states`) once, on the first (topology)
  frame of a multi-frame `model_system`; never per frame. Per-frame duplication
  bloats the archive and the Elasticsearch index document and can fail an upload
  (FAIRmat-NFDI/nomad-simulations#474). Deviate only for genuinely frame-varying
  identity (reactive / alchemical simulations). Enforced by
  `tests/parsers/common.py::assert_identity_populated_once`.
- Put shared test assertions in `tests/parsers/common.py`, not in new files.
