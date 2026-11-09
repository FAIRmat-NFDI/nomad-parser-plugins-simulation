# Bundled LAMMPS fixtures

This directory contains the small, versioned inputs used by
`tests/parsers/lammps/`. Use the **main log** below as the parser mainfile;
the other log, data, and trajectory files in its directory are dependencies.
The parser may discover files by commands in the log, so moving an individual
mainfile out of its directory changes the test.

| Directory | Main log and purpose | Expected archive-level outcome |
| --- | --- | --- |
| `1_methyl_naphthalene/` | `log.1_methyl_naphthalene`; redirects to `log.naph_298_396_20ns` and reads `data.1_methyl_naphthalene` plus `naph_298_eq.lammpstrj`. | Physical `real` units, NPT method, 4 systems of 1,134 particles. The input structure, dump, and redirected log are present. |
| `hexane_cyclohexane/` | `log.hexane_cyclohexane_nvt`; reads `data.hexane_cyclohexane`, the `.lammpstrj` dump, and the redirected `-thermo_style_multi` log. | Physical `real` units, NVT method, 201 systems of 684 particles. |
| `methane_dcd/` | `log.methane_nvt_traj_dcd_thermo_style_custom`; reads `data.64xmethane_from_restart`, a DCD dump, and the redirected log. | Physical `real` units, NVT method, 201 systems of 320 particles. |
| `methane_xyz/` | `log.methane_nvt_traj_xyz_thermo_style_custom`; reads `data.64xmethane_from_restart`, an XYZ dump, and the redirected log. | Physical `real` units, NVT method, 201 systems of 320 particles. |
| `1_xyz_files/` | `log.lammps`; one custom dump combines position, velocity, and force fields. Its local `READMe` describes how the original LAMMPS run was invoked. | Explicit **`units lj`**. Keep as a low-level reader fixture; archive writing must reject it before assigning physical units. |
| `2_xyz_files/` | `log.lammps`; separate custom position and velocity dumps. Its local `READMe` describes how the original run was invoked. | Explicit **`units lj`**. Keep as a low-level/multiple-dump reader fixture; archive writing must reject it before assigning physical units. |

The `polymer_melt/` fixture referenced by marked tests is **not bundled
here**. It lives in the external
[`nomad-simulation-parser-test-fixtures`](https://github.com/FAIRmat-NFDI/nomad-simulation-parser-test-fixtures)
repository. Configure `NOMAD_SIM_PARSERS_LARGE_FIXTURE_ROOT` to its fixture
root before running `large_fixture` tests; otherwise use
`-m 'not large_fixture'`. Do not infer a parser regression from a missing
external fixture.

The data-only and missing-file cases in `test_input_structure.py` are
generated in temporary directories by the tests; they are **not** additional
files missing from this directory. Those examples deliberately distinguish
an available LAMMPS *input* structure from an unavailable final structure.
Likewise, `log.test` in `test_reader.py` is a synthetic path used to test
filename selection, not a missing fixture.

## Central diagnostic cases are not bundled

The following public NOMAD Central entries informed the reader tests, but
their raw files are **not** package fixtures or CI inputs. They are selected
in the opt-in local cache at `distro-devs/testing_data_from_central/lammps/`;
its `README.md`, `selection.jsonl`, and `manifest.json` record provenance and
which raw files were actually downloaded. Do not assume a referenced file is
available just because it appears in a LAMMPS command.

| Case (entry ID) | Present in the local cache | Missing or not usable; expected interpretation |
| --- | --- | --- |
| BENZENE12 (`-e3gL2URfcOkRj-TBikJUkXHuAHY`) | `log.lammps`, `lmp.data` | No trajectory was written; final `data.naph_post` is not cached. The 48-particle system recovered from `read_data` is an **input structure**, not the minimized result. |
| Luca (`HJqVgYr7onG0sP3hHREtf2mN5jkO`) | Startup `log.lammps`, redirected `log.system`, `data_start.in` | Referenced `starting_structure.xyz`, `relax_traj.dcd`, and final `data.output_1` are not cached. The 44,067-particle system recovered from `read_data` is an **input structure**, not an MD frame or final minimized result. |
| Sandro (`2sU_Yal29UcX_aUIKXEd21eVXIPl`) | `lmps.log` | `read_data` names an unavailable absolute `opt.data` path, so there is no source-backed particle system. Its two one-step force evaluations are a later multi-calculation/workflow case; a one-step MD label is a known misclassification, not a passing expectation. |

The package tests use small generated analogues of these conditions; they do
not assert the Central entries themselves. See the [LAMMPS test matrix](../../parsers/lammps/README.md)
for the corresponding tests and scope.
