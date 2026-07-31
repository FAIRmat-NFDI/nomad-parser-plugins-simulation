# Exchange–correlation functionals

## Overview

Every DFT code records the exchange–correlation (XC) functional in its own way —
an input tag, an integer code, or a human-readable line in the output — and only
rarely as a single canonical name. The simulation parsers reduce that variety to
one of two forms that the `nomad-simulations` schema understands, and leave the
[LibXC](https://libxc.gitlab.io/functionals/) taxonomy (family, kind, and the
derived Jacob's ladder rung) to the schema rather than reconstructing it per
parser.

A parser produces either:

- a **canonical name** in `XCFunctional.functional_key` — used when the code
  names the functional, or when its tags reduce cleanly to a standard name
  (`PBE`, `HSE06`, `B3LYP5`, …). The schema's alias table
  (`nomad_simulations/.../libxc/aliases.json`) expands the name into LibXC
  components; or
- **free-form components** in `XCFunctional.components` — used when the code
  reports the exchange and correlation parts individually (as LibXC labels or
  ids). The schema resolves each component from the LibXC registry, filling its
  family and kind.

The choice reflects what the code actually reports: a name where one exists, the
raw parts otherwise. Either way the parser never hard-codes LibXC family/kind
data. A component that the registry cannot resolve is retained and flagged
`unidentified` instead of being dropped, so a functional's composition is never
silently misrepresented as complete.

This page explains the non-obvious conventions behind each code's mapping. It is
deliberately not an exhaustive list of supported functionals — the maps in the
parser sources are authoritative for that. It focuses on the setups where the
code's XC specification is indirect or easy to misread.

## How each code expresses XC

### VASP

VASP emits no functional name, so the canonical name is reconstructed from the
input tags in the order VASP itself resolves them: an explicit hybrid setup
wins, then `METAGGA`, then `GGA`.

The `GGA` and `METAGGA` tags are short two-letter codes (`PE` → PBE, `PS` →
PBEsol, `RP` → RPBE, `RE` → revPBE, …). Two subtleties are worth stating:

- **Hybrids** (`LHFCALC = .TRUE.`) are decided from `HFSCREEN`, `GGA` and `AEXX`
  together, not from the screening length alone. `HFSCREEN = 0.2`/`0.3` give the
  screened HSE hybrids, but the base GGA sets the variant: a PBE base gives
  HSE06/HSE03, whereas `GGA = PS` (PBEsol) gives **HSEsol**. The unscreened
  global hybrid **PBE0** additionally requires that screening is off. `AEXX = 1`
  with no DFT correlation is pure Hartree–Fock, not a DFT functional.
- `GGA = B3` and `GGA = B5` are both B3LYP but with **different correlation** —
  VWN3 (`B3` → `B3LYP`) versus VWN5 (`B5` → `B3LYP5`) — so they must not collapse
  to one key.
- When `GGA` is unset the effective functional is the **POTCAR default**, carried
  by `LEXCH` (e.g. `LEXCH = CA` is the Ceperley–Alder LDA in the Perdew–Zunger
  parametrisation, `PZ81`). PBE is assumed only when neither `GGA` nor `LEXCH`
  is present.

### ABINIT

ABINIT identifies XC by the single integer `ixc`. A non-negative `ixc` selects a
built-in preset (`1` = Teter93 LDA, `11` = PBE, `14` = revPBE, `15` = RPBE,
`23` = Wu–Cohen, …), mapped through a small table. A **negative** `ixc` packs two
LibXC functional ids positionally as `-(1000·id_x + id_c)` — for example
`ixc = -101130` is LibXC `101` (`XC_GGA_X_PBE`) plus `130` (`XC_GGA_C_PBE`),
i.e. PBE. The parser emits those ids as components and lets the schema registry
resolve them, which also covers ids the local table does not list.

### Quantum ESPRESSO

QE prints XC either as a single name (`PBE`) or as its four DFT slot codes —
exchange, correlation, gradient-exchange, gradient-correlation. `SLA PW PBX PBC`
is Slater exchange + Perdew–Wang correlation + PBE gradient terms, i.e. PBE;
`SLA PZ` is the Perdew–Zunger LDA. The parser recognises the common slot
combinations and, for anything else, keeps the raw slot string so the reported
XC is preserved rather than dropped.

### octopus

octopus prints the exchange and correlation parts as separate descriptive lines
(`Slater exchange`, `Perdew & Zunger (Modified)`), which the parser maps to LibXC
short-labels. The `XCFunctional` input variable is additionally an integer that
is the *sum* of per-piece codes; when present it is decoded as a fallback to
recover parts the descriptions did not name.

### exciting

exciting's INFO file gives an integer `xctype` code that the parser maps to LibXC
labels; the XML output gives LibXC labels directly. The integer codes track a
specific exciting/LibXC version (see *Known limitations*).

### CRYSTAL

CRYSTAL names the exchange and correlation keywords separately (`BECKE`, `PZ`,
`PBE`, `PWGGA`, …); the parser maps each keyword to its LibXC label(s) and lets
the schema assemble the functional.

### AMS / ADF

AMS reports the active functional per Jacob's-ladder rung in its density-functional
block (`LDA:`, `Gradient Corrections:`, `Meta-GGA:`). The highest present rung is
the functional in effect. A two-word gradient-correction entry names the exchange
and correlation authors together (`Becke Perdew` = BP86).

### GPAW and FHI-aims

GPAW records a standard functional name directly and the parser passes it
through unchanged. FHI-aims prints a descriptive control string per functional,
which the parser maps to the corresponding standard name.

## Known limitations

Some codes tie their XC identifiers to a specific LibXC (or code) release, so a
mapping that is correct for one version can drift in another — ABINIT's negative
`ixc` LibXC ids and exciting's integer `xctype` codes are the main examples.
Mappings should be cross-checked against the LibXC version the code is built with
when in doubt. Open verification items are tracked in the parser repository's
issues.

## References

The LibXC label taxonomy and each code's XC-input documentation:

- LibXC functional list — <https://libxc.gitlab.io/functionals/>
- VASP — [`GGA`](https://www.vasp.at/wiki/index.php/GGA),
  [`METAGGA`](https://www.vasp.at/wiki/index.php/METAGGA),
  [`LHFCALC`](https://www.vasp.at/wiki/index.php/LHFCALC) and the
  [list of hybrid functionals](https://www.vasp.at/wiki/index.php/List_of_hybrid_functionals)
- ABINIT — the [`ixc`](https://docs.abinit.org/variables/basic/#ixc) input variable
- Quantum ESPRESSO — [`input_dft`](https://www.quantum-espresso.org/Doc/INPUT_PW.html)
  (the slot codes are defined in `Modules/funct.f90` / `XClib`)
- octopus — the [`XCFunctional`](https://octopus-code.org/documentation/13/variables/hamiltonian/xc/xcfunctional/) variable
- exciting — the [`groundstate`/`xctype`](http://exciting-code.org/ref:groundstate) reference
- CRYSTAL — the `DFT`/`EXCHANGE`/`CORRELAT` keywords in the
  [CRYSTAL23 user's manual](https://www.crystal.unito.it/include/manuals/crystal23.pdf)
- AMS / ADF — [Density Functionals (XC)](https://www.scm.com/doc/ADF/Input/Density_Functional.html)
- GPAW — [exchange–correlation functionals](https://gpaw.readthedocs.io/documentation/xc/functionals.html)
- FHI-aims — the `xc` keyword ([manual §3.3](https://fhi-aims.org/uploads/manual/Ch3/S3.html))
