"""Shared VASP helpers used by both the OUTCAR and vasprun.xml parsers."""

from typing import Any

XC_FUNCTIONAL_MAPPING = {
    '--': ['GGA_X_PBE', 'GGA_C_PBE'],
    'HL': ['LDA_C_HL'],
    'WI': ['LDA_C_WIGNER'],
    'PZ': ['LDA_C_PZ'],
    '91': ['GGA_X_PW91', 'GGA_C_PW91'],
    'PE': ['GGA_X_PBE', 'GGA_C_PBE'],
    'PBE': ['GGA_X_PBE', 'GGA_C_PBE'],
    'RE': ['GGA_X_PBE_R'],
    'VW': ['LDA_C_VWN'],
    'RP': ['GGA_X_RPBE', 'GGA_C_PBE'],
    'PS': ['GGA_C_PBE_SOL', 'GGA_X_PBE_SOL'],
    'AM': ['GGA_X_AM05', 'GGA_C_AM05'],
    'B3': ['HYB_GGA_XC_B3LYP3'],
    'B5': ['HYB_GGA_XC_B3LYP5'],
    'BF': ['GGA_X_BEEFVDW', 'GGA_XC_BEEFVDW'],
    'CO': [],  # TODO check if this is ever used
    'OR': ['GGA_X_OPTPBE_VDW'],
    'BO': ['GGA_X_OPTB88_VDW'],
    'MK': ['GGA_X_OPTB86B_VDW'],
    'ML': ['VDW_XC_DF2'],
    'CX': ['VDW_XC_DF_CX'],
    'TPSS': ['MGGA_X_TPSS', 'MGGA_C_TPSS'],
    'RTPSS': ['MGGA_X_RTPSS'],
    'M06L': ['MGGA_C_M06_L'],
    'MS0': ['MGGA_X_MS0'],
    'MS1': ['MGGA_X_MS1'],
    'MS2': ['MGGA_X_MS2'],
    'SCAN': ['MGGA_X_SCAN'],
    'RSCAN': ['MGGA_X_RSCAN', 'MGGA_C_RSCAN'],
    'R2SCAN': ['MGGA_X_R2SCAN', 'MGGA_C_R2SCAN'],
    'SCANL': ['MGGA_X_SCANL', 'MGGA_C_SCANL'],
    'RSCANL': [],  # not in LibXC, nor any paper, just deorbitalized SCANL
    'R2SCANL': ['MGGA_X_R2SCANL', 'MGGA_C_R2SCANL'],
    'OFR2': [],
    'MBJ': ['MGGA_X_BJ06'],
    'LBMJ': [],  # TODO ask Miguel Marquez
    'HLE17': ['MGGA_XC_HLE17'],  # TODO check if this is ever used
    'RA': ['LDA_C_PW_RPA'],  # TODO check if this is ever used
}


def get_xc_functionals(parameters: dict[str, Any]) -> list[dict[str, Any]]:
    xc_functionals = []
    if parameters.get('LHFCALC', False):
        hfscreen_06, hfscreen_03 = 0.2, 0.3
        aexx_b3, aggax_b3, aggac_b3, aldac_b3 = 0.2, 0.72, 0.81, 0.19
        xc_functional = {}
        gga = parameters.get('GGA', 'PE')
        aexx = parameters.get('AEXX', 0.0)
        aggax = parameters.get('AGGAX', 1.0)
        aggac = parameters.get('AGGAC', 1.0)
        aldac = parameters.get('ALDAC', 1.0)
        hfscreen = parameters.get('HFSCREEN', 0.0)

        if hfscreen == hfscreen_06:
            xc_functional['name'] = 'HYB_GGA_XC_HSE06'
        elif hfscreen == hfscreen_03:
            xc_functional['name'] = 'HYB_GGA_XC_HSE03'
        elif (
            gga == 'B3'
            and aexx == aexx_b3
            and aggax == aggax_b3
            and aggac == aggac_b3
            and aldac == aldac_b3
        ):
            xc_functional['name'] = 'HYB_GGA_XC_B3LYP3'
        elif aexx == 1.0 and aldac == 0.0 and aggac == 0.0:
            xc_functional['name'] = 'HF_X'
        elif gga == 'PE':
            xc_functional['name'] = 'HYB_GGA_XC_PBEH'
        else:
            xc_functional['name'] = f'HYB_GGA_XC_{gga}'
        xc_functionals.append(xc_functional)
    else:
        metagga = parameters.get('METAGGA')
        if metagga:
            functionals = XC_FUNCTIONAL_MAPPING.get(metagga, [metagga])
        else:
            # VASP defaults to PBE-like GGA if GGA is not explicitly set.
            functionals = XC_FUNCTIONAL_MAPPING.get(parameters.get('GGA', 'PE'), [])
        for functional in functionals:
            xc_functionals.append({'name': functional})
    return xc_functionals
