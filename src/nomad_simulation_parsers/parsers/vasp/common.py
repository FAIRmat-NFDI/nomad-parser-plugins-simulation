"""Shared VASP helpers used by both the OUTCAR and vasprun.xml parsers."""

from typing import Any

from nomad.utils import get_logger

LOGGER = get_logger(__name__)

VASP_TAG_TO_FUNCTIONAL = {
    # GGA tag (`--` = POTCAR default is handled by the `PE` fallback)
    'PE': 'PBE',
    'PBE': 'PBE',
    'PS': 'PBEsol',
    'RP': 'RPBE',
    'RE': 'revPBE',
    '91': 'PW91',
    'AM': 'AM05',
    'HL': 'HL',
    'WI': 'Wigner',
    'PZ': 'PZ81',
    'VW': 'VWN',
    # LEXCH-only tag: POTCAR-default LDA (Ceperley-Alder, Perdew-Zunger)
    'CA': 'PZ81',
    'B3': 'B3LYP',  # B3LYP with VWN3 correlation
    'B5': 'B3LYP5',  # B3LYP with VWN5 correlation
    'OR': 'optPBE-vdW',
    'BO': 'optB88-vdW',
    'MK': 'optB86b-vdW',
    'BF': 'BEEF-vdW',
    # METAGGA tag
    'TPSS': 'TPSS',
    'RTPSS': 'revTPSS',
    'M06L': 'M06-L',
    'MS0': 'MS0',
    'MS1': 'MS1',
    'MS2': 'MS2',
    'SCAN': 'SCAN',
    'RSCAN': 'RSCAN',
    'R2SCAN': 'R2SCAN',
    'SCANL': 'SCAN-L',
    'R2SCANL': 'R2SCAN-L',
    'MBJ': 'MBJ',
}


def _clean_tag(value: Any) -> str | None:
    """Return the XC tag as a clean string, or `None` when it is unset.

    The two sources spell "unset" differently: OUTCAR yields a non-string (bool
    `False`, from the shared parameter converter), vasprun yields an absent
    element or a `--`/`NONE` sentinel. All of them normalize to `None` here.
    """
    # TODO: OUTCAR reports an unset tag as bool `False` because `get_key_values`
    # converts every logical token; normalizing unset XC tags to `None` on the
    # OUTCAR side would let the non-string guard go, at the cost of
    # METAGGA-specific handling around that shared converter.
    if not isinstance(value, str):
        return None
    tag = value.strip().strip('"').strip()
    if not tag or tag.upper() in ('--', 'NONE'):
        return None
    return tag


def _hybrid_functional_key(parameters: dict[str, Any]) -> str | None:
    """Canonical functional name for a VASP hybrid run (LHFCALC set). Both source
    parsers type their parameters, so values arrive as bool/float already.

    The variant is decided from `HFSCREEN`, `GGA` and `AEXX` together, not the
    screening length alone: the screened HSE hybrids take their base GGA (PBE ->
    HSE, PBEsol -> HSEsol), and the unscreened global hybrid PBE0 requires that
    screening is off.
    """
    gga = _clean_tag(parameters.get('GGA'))
    aexx = parameters.get('AEXX') or 0.0
    aggac = parameters.get('AGGAC')
    aldac = parameters.get('ALDAC')
    # an absent HFSCREEN normalizes to VASP's default, 0 (unscreened)
    hfscreen = parameters.get('HFSCREEN') or 0.0
    pbe_family = gga in (None, 'PE', 'PBE')

    # B3LYP: VWN3 (GGA=B3) vs VWN5 (GGA=B5)
    if gga in ('B3', 'B5'):
        return 'B3LYP' if gga == 'B3' else 'B3LYP5'
    if aexx == 1.0 and aldac == 0.0 and aggac == 0.0:
        return None  # pure Hartree-Fock exchange, not a DFT functional
    # Screened hybrids (HSE): the base GGA sets the variant (PBE -> HSE06/03,
    # PBEsol -> HSEsol); any other base is left unresolved.
    if hfscreen == 0.2:  # noqa: PLR2004 — HSE06 screening length (1/Angstrom)
        return 'HSE06' if pbe_family else ('HSEsol' if gga == 'PS' else None)
    if hfscreen == 0.3:  # noqa: PLR2004 — HSE03 screening length (1/Angstrom)
        return 'HSE03' if pbe_family else None
    # Unscreened global hybrid PBE0: screening off (HFSCREEN 0 or unset) and a
    # PBE-family base GGA.
    if hfscreen == 0.0 and pbe_family:
        return 'PBE0'
    return None


def functional_key_from_params(parameters: dict[str, Any]) -> str | None:
    """Map VASP XC input tags to a canonical functional name.

    vasprun.xml (and the OUTCAR) carry no XC/DFT/method section and VASP emits no
    functional name, so the canonical name is reconstructed here from the input
    tags. Precedence follows VASP's own resolution: an explicit
    hybrid/Hartree-Fock setup wins, then `METAGGA`, then `GGA`. When `GGA` is
    unset the effective functional is the POTCAR default, carried by `LEXCH`
    (e.g. `CA` for an LDA POTCAR); PBE is assumed only when neither is present.
    The schema expands the returned name into LibXC components during
    normalization. The per-source `get_functional_key` transformer methods
    assemble `parameters` and delegate here.

    The precedence (METAGGA over GGA, and LHFCALC selecting a hybrid) follows the
    VASP tag documentation: https://www.vasp.at/wiki/index.php/METAGGA and
    https://www.vasp.at/wiki/index.php/LHFCALC.
    """
    if parameters.get('LHFCALC'):
        return _hybrid_functional_key(parameters)
    if metagga := _clean_tag(parameters.get('METAGGA')):
        return _mapped_functional(metagga, 'METAGGA')
    if gga := _clean_tag(parameters.get('GGA')):
        return _mapped_functional(gga, 'GGA')
    if lexch := _clean_tag(parameters.get('LEXCH')):
        return _mapped_functional(lexch, 'LEXCH')
    return _mapped_functional('PE', 'GGA')


def _mapped_functional(tag: str, source: str) -> str | None:
    """Resolve a VASP tag to a canonical name, logging a present-but-unmapped tag.

    Returns `None` for an unrecognized tag rather than the raw tag: `functional_key`
    must be a canonical name the schema can expand, and VASP's short tags are not.
    """
    name = VASP_TAG_TO_FUNCTIONAL.get(tag)
    if name is None:
        LOGGER.debug(
            'unmapped VASP %s tag %s; leaving functional_key unset', source, tag
        )
    return name
