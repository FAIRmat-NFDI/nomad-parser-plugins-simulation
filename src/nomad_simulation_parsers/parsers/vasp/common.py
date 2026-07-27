"""Shared VASP helpers used by both the OUTCAR and vasprun.xml parsers."""

from typing import Any

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
    'B3': 'B3LYP',
    'B5': 'B3LYP',
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

_HFSCREEN_HSE06, _HFSCREEN_HSE03 = 0.2, 0.3


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
    parsers type their parameters, so values arrive as bool/float already."""
    gga = _clean_tag(parameters.get('GGA'))
    aexx = parameters.get('AEXX') or 0.0
    aggac = parameters.get('AGGAC')
    aldac = parameters.get('ALDAC')
    hfscreen = parameters.get('HFSCREEN') or 0.0

    if hfscreen == _HFSCREEN_HSE06:
        return 'HSE06'
    if hfscreen == _HFSCREEN_HSE03:
        return 'HSE03'
    if gga in ('B3', 'B5'):
        return 'B3LYP'
    if aexx == 1.0 and aldac == 0.0 and aggac == 0.0:
        return None  # pure Hartree-Fock exchange, not a DFT functional
    # Default hybrid is PBE0: GGA either unset (POTCAR default PBE) or PBE-family.
    if gga in (None, 'PE', 'PBE'):
        return 'PBE0'
    return None


def functional_key_from_params(parameters: dict[str, Any]) -> str | None:
    """Map VASP XC input tags to a canonical functional name.

    vasprun.xml (and the OUTCAR) carry no XC/DFT/method section and VASP emits no
    functional name, so the canonical name is reconstructed here from the input
    tags. Precedence follows VASP's own resolution: an explicit
    hybrid/Hartree-Fock setup wins, then `METAGGA`, then `GGA` (defaulting to PBE
    when unset). The schema expands the returned name into LibXC components
    during normalization. The per-source `get_functional_key` transformer methods
    assemble `parameters` and delegate here.
    """
    if parameters.get('LHFCALC'):
        return _hybrid_functional_key(parameters)
    if metagga := _clean_tag(parameters.get('METAGGA')):
        return VASP_TAG_TO_FUNCTIONAL.get(metagga)
    gga = _clean_tag(parameters.get('GGA')) or 'PE'
    return VASP_TAG_TO_FUNCTIONAL.get(gga)
