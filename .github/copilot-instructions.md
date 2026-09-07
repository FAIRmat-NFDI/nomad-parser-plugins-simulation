# Copilot / AI agent instructions

## Store per-particle identity once, not per frame

When a parser emits a multi-frame `model_system` sequence — an MD trajectory or a
geometry-optimization / ionic-step series — the per-particle identity
(`particle_states`: `chemical_symbol`, `atomic_number`, `label`, `mass`, ...) is
frame-independent. Populate it on the **first (topology) frame only**; the
remaining frames must carry only what varies in time (positions, cell,
velocities).

Do not attach `particle_states` to every frame. Every scalar under
`archive.data` is auto-registered as an Elasticsearch `search_quantity`, so
duplicating identity scales both the archive and the index document as
`n_frames × n_particles × fields`. For long trajectories this exceeds the
Elasticsearch payload limit and fails the whole upload (see
FAIRmat-NFDI/nomad-simulations#474; the LAMMPS methane example produced a 127 MB
index document before the fix). Prefer withholding identity at the source
(e.g. attach labels only for the first frame in the parser's configuration
builder) over removing it afterwards.

Deviate — populating identity on more than one frame — only when it genuinely
changes between frames (reactive or alchemical simulations where species are
created, destroyed, or transmuted), and only for the frames that change. When in
doubt, store it once.

The invariant is enforced by
`tests/parsers/_assertions.py::assert_identity_populated_once`, which asserts that
exactly one `model_system` frame carries `particle_states`. Wire it into a
parser's multi-frame test (or rely on `SimulationParserTestSuite`, which includes
it) whenever you add or change trajectory / optimization handling.
