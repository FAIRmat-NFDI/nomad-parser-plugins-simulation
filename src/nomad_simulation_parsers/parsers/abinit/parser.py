import os
from collections.abc import Iterable
from datetime import datetime
from typing import Any

import numpy as np
from ase.data import chemical_symbols
from nomad.datamodel import EntryArchive
from nomad.parsing.parser import MatchingParser
from nomad.units import ureg
from nomad.utils import get_logger
from nomad_file_parser import ArchiveWriter, DataTextParser
from nomad_file_parser.mapping_parser import MetainfoParser, TextParser
from nomad_simulations.schema_packages.general import Program, Simulation
from nomad_simulations.schema_packages.workflow import (
    DFTGWWorkflow,
    MolecularDynamics,
    SinglePoint,
)
from nomad_simulations.schema_packages.workflow.general import (
    EnergyConvergenceTarget,
    ForceConvergenceTarget,
)
from nomad_simulations.schema_packages.workflow.geometry_optimization import (
    GeometryOptimization,
    GeometryOptimizationMethod,
)
from nomad_simulations.schema_packages.workflow.single_point import SinglePointMethod
from structlog.stdlib import BoundLogger

from nomad_simulation_parsers.parsers.utils.general import (
    calculate_band_gap_from_occupations,
)
from nomad_simulation_parsers.schema_packages import abinit

from .file_parser import AbinitOutParser

LOGGER = get_logger(__name__)


ABINIT_NATIVE_IXC = {
    0: [{}],
    1: [{'XC_functional_name': 'LDA_XC_TETER93'}],
    2: [{'XC_functional_name': 'LDA_X'}, {'XC_functional_name': 'LDA_C_PZ'}],
    # 3 - LDA, old Teter rational polynomial parametrization (4/91)
    4: [{'XC_functional_name': 'LDA_X'}, {'XC_functional_name': 'LDA_C_WIGNER'}],
    5: [{'XC_functional_name': 'LDA_X'}, {'XC_functional_name': 'LDA_C_HL'}],
    6: [{'XC_functional_name': 'LDA_X'}, {'XC_functional_name': 'LDA_C_XALPHA'}],
    7: [{'XC_functional_name': 'LDA_X'}, {'XC_functional_name': 'LDA_C_PW'}],
    # 8 - x-only part of the Perdew-Wang 92 functional
    # 9 - x- and RPA correlation part of the Perdew-Wang 92 functional
    # 10 - non-existent
    11: [{'XC_functional_name': 'GGA_X_PBE'}, {'XC_functional_name': 'GGA_C_PBE'}],
    12: [{'XC_functional_name': 'GGA_X_PBE'}],
    13: [{'XC_functional_name': 'GGA_X_LB'}, {'XC_functional_name': 'LDA_C_PW'}],
    14: [{'XC_functional_name': 'GGA_X_PBE_R'}, {'XC_functional_name': '?'}],
    15: [{'XC_functional_name': 'GGA_X_RPBE'}, {'XC_functional_name': '?'}],
    16: [{'XC_functional_name': 'GGA_XC_HCTH_93'}],
    17: [{'XC_functional_name': 'GGA_XC_HCTH_120'}],
    18: [{'XC_functional_name': 'GGA_X_B88'}, {'XC_functional_name': 'GGA_C_LYP'}],
    19: [{'XC_functional_name': 'GGA_X_B88'}, {'XC_functional_name': 'GGA_C_P86'}],
    # 20 - Fermi-Amaldi xc ( -1/N Hartree energy, where N is the number of electrons
    # per cell;
    #      G=0 is not taken into account however), for TDDFT tests.
    # 21 - same as 20, except that the xc-kernel is the LDA (ixc=1) one, for TDDFT
    # tests.
    # 22 - same as 20, except that the xc-kernel is the Burke-Petersilka-Gross hybrid,
    # for TDDFT tests.
    23: [{'XC_functional_name': 'GGA_X_WC'}, {'XC_functional_name': '?'}],
    24: [{'XC_functional_name': 'GGA_X_C09X'}, {'XC_functional_name': '?'}],
    # 25 - non-existent
    26: [{'XC_functional_name': 'GGA_XC_HCTH_147'}],
    27: [{'XC_functional_name': 'GGA_XC_HCTH_407'}],
    28: [{'XC_functional_name': 'GGA_X_OPTX'}, {'XC_functional_name': 'GGA_C_LYP'}],
    40: [{'XC_functional_name': 'HF_X'}],
    41: [{'XC_functional_name': 'HYB_GGA_XC_PBEH'}],
    42: [{'XC_functional_name': 'HYB_GGA_XC_PBE0_13'}],
}

ABINIT_LIBXC_IXC = {
    1: {'XC_functional_name': 'LDA_X'},
    2: {'XC_functional_name': 'LDA_C_WIGNER'},
    3: {'XC_functional_name': 'LDA_C_RPA'},
    4: {'XC_functional_name': 'LDA_C_HL'},
    5: {'XC_functional_name': 'LDA_C_GL'},
    6: {'XC_functional_name': 'LDA_C_XALPHA'},
    7: {'XC_functional_name': 'LDA_C_VWN'},
    8: {'XC_functional_name': 'LDA_C_VWN_RPA'},
    9: {'XC_functional_name': 'LDA_C_PZ'},
    10: {'XC_functional_name': 'LDA_C_PZ_MOD'},
    11: {'XC_functional_name': 'LDA_C_OB_PZ'},
    12: {'XC_functional_name': 'LDA_C_PW'},
    13: {'XC_functional_name': 'LDA_C_PW_MOD'},
    14: {'XC_functional_name': 'LDA_C_OB_PW'},
    15: {'XC_functional_name': 'LDA_C_2D_AMGB'},
    16: {'XC_functional_name': 'LDA_C_2D_PRM'},
    17: {'XC_functional_name': 'LDA_C_vBH'},
    18: {'XC_functional_name': 'LDA_C_1D_CSC'},
    19: {'XC_functional_name': 'LDA_X_2D'},
    20: {'XC_functional_name': 'LDA_XC_TETER93'},
    21: {'XC_functional_name': 'LDA_X_1D'},
    22: {'XC_functional_name': 'LDA_C_ML1'},
    23: {'XC_functional_name': 'LDA_C_ML2'},
    24: {'XC_functional_name': 'LDA_C_GOMBAS'},
    25: {'XC_functional_name': 'LDA_C_PW_RPA'},
    26: {'XC_functional_name': 'LDA_C_1D_LOOS'},
    27: {'XC_functional_name': 'LDA_C_RC04'},
    28: {'XC_functional_name': 'LDA_C_VWN_1'},
    29: {'XC_functional_name': 'LDA_C_VWN_2'},
    30: {'XC_functional_name': 'LDA_C_VWN_3'},
    31: {'XC_functional_name': 'LDA_C_VWN_4'},
    32: {'XC_functional_name': 'GGA_X_GAM'},
    33: {'XC_functional_name': 'GGA_C_GAM'},
    34: {'XC_functional_name': 'GGA_X_HCTH_A'},
    35: {'XC_functional_name': 'GGA_X_EV93'},
    36: {'XC_functional_name': 'HYB_MGGA_X_DLDF'},
    37: {'XC_functional_name': 'MGGA_C_DLDF'},
    38: {'XC_functional_name': 'GGA_X_BCGP'},
    39: {'XC_functional_name': 'GGA_C_BCGP'},
    40: {'XC_functional_name': 'GGA_X_LAMBDA_OC2_N'},
    41: {'XC_functional_name': 'GGA_X_B86_R'},
    42: {'XC_functional_name': 'MGGA_XC_ZLP'},
    43: {'XC_functional_name': 'LDA_XC_ZLP'},
    44: {'XC_functional_name': 'GGA_X_LAMBDA_CH_N'},
    45: {'XC_functional_name': 'GGA_X_LAMBDA_LO_N'},
    46: {'XC_functional_name': 'GGA_X_HJS_B88_V2'},
    47: {'XC_functional_name': 'GGA_C_Q2D'},
    48: {'XC_functional_name': 'GGA_X_Q2D'},
    49: {'XC_functional_name': 'GGA_X_PBE_MOL'},
    50: {'XC_functional_name': 'LDA_K_TF'},
    51: {'XC_functional_name': 'LDA_K_LP'},
    52: {'XC_functional_name': 'GGA_K_TFVW'},
    53: {'XC_functional_name': 'GGA_K_REVAPBEINT'},
    54: {'XC_functional_name': 'GGA_K_APBEINT'},
    55: {'XC_functional_name': 'GGA_K_REVAPBE'},
    56: {'XC_functional_name': 'GGA_X_AK13'},
    57: {'XC_functional_name': 'GGA_K_MEYER'},
    58: {'XC_functional_name': 'GGA_X_LV_RPW86'},
    59: {'XC_functional_name': 'GGA_X_PBE_TCA'},
    60: {'XC_functional_name': 'GGA_X_PBEINT'},
    61: {'XC_functional_name': 'GGA_C_ZPBEINT'},
    62: {'XC_functional_name': 'GGA_C_PBEINT'},
    63: {'XC_functional_name': 'GGA_C_ZPBESOL'},
    64: {'XC_functional_name': 'MGGA_XC_OTPSS_D'},
    65: {'XC_functional_name': 'GGA_XC_OPBE_D'},
    66: {'XC_functional_name': 'GGA_XC_OPWLYP_D'},
    67: {'XC_functional_name': 'GGA_XC_OBLYP_D'},
    68: {'XC_functional_name': 'GGA_X_VMT84_GE'},
    69: {'XC_functional_name': 'GGA_X_VMT84_PBE'},
    70: {'XC_functional_name': 'GGA_X_VMT_GE'},
    71: {'XC_functional_name': 'GGA_X_VMT_PBE'},
    72: {'XC_functional_name': 'MGGA_C_CS'},
    73: {'XC_functional_name': 'MGGA_C_MN12_SX'},
    74: {'XC_functional_name': 'MGGA_C_MN12_L'},
    75: {'XC_functional_name': 'MGGA_C_M11_L'},
    76: {'XC_functional_name': 'MGGA_C_M11'},
    77: {'XC_functional_name': 'MGGA_C_M08_SO'},
    78: {'XC_functional_name': 'MGGA_C_M08_HX'},
    79: {'XC_functional_name': 'GGA_C_N12_SX'},
    80: {'XC_functional_name': 'GGA_C_N12'},
    81: {'XC_functional_name': 'HYB_GGA_X_N12_SX'},
    82: {'XC_functional_name': 'GGA_X_N12'},
    83: {'XC_functional_name': 'GGA_C_REGTPSS'},
    84: {'XC_functional_name': 'GGA_C_OP_XALPHA'},
    85: {'XC_functional_name': 'GGA_C_OP_G96'},
    86: {'XC_functional_name': 'GGA_C_OP_PBE'},
    87: {'XC_functional_name': 'GGA_C_OP_B88'},
    88: {'XC_functional_name': 'GGA_C_FT97'},
    89: {'XC_functional_name': 'GGA_C_SPBE'},
    90: {'XC_functional_name': 'GGA_X_SSB_SW'},
    91: {'XC_functional_name': 'GGA_X_SSB'},
    92: {'XC_functional_name': 'GGA_X_SSB_D'},
    93: {'XC_functional_name': 'GGA_XC_HCTH_407P'},
    94: {'XC_functional_name': 'GGA_XC_HCTH_P76'},
    95: {'XC_functional_name': 'GGA_XC_HCTH_P14'},
    96: {'XC_functional_name': 'GGA_XC_B97_GGA1'},
    97: {'XC_functional_name': 'GGA_C_HCTH_A'},
    98: {'XC_functional_name': 'GGA_X_BPCCAC'},
    99: {'XC_functional_name': 'GGA_C_REVTCA'},
    100: {'XC_functional_name': 'GGA_C_TCA'},
    101: {'XC_functional_name': 'GGA_X_PBE'},
    102: {'XC_functional_name': 'GGA_X_PBE_R'},
    103: {'XC_functional_name': 'GGA_X_B86'},
    104: {'XC_functional_name': 'GGA_X_HERMAN'},
    105: {'XC_functional_name': 'GGA_X_B86_MGC'},
    106: {'XC_functional_name': 'GGA_X_B88'},
    107: {'XC_functional_name': 'GGA_X_G96'},
    108: {'XC_functional_name': 'GGA_X_PW86'},
    109: {'XC_functional_name': 'GGA_X_PW91'},
    110: {'XC_functional_name': 'GGA_X_OPTX'},
    111: {'XC_functional_name': 'GGA_X_DK87_R1'},
    112: {'XC_functional_name': 'GGA_X_DK87_R2'},
    113: {'XC_functional_name': 'GGA_X_LG93'},
    114: {'XC_functional_name': 'GGA_X_FT97_A'},
    115: {'XC_functional_name': 'GGA_X_FT97_B'},
    116: {'XC_functional_name': 'GGA_X_PBE_SOL'},
    117: {'XC_functional_name': 'GGA_X_RPBE'},
    118: {'XC_functional_name': 'GGA_X_WC'},
    119: {'XC_functional_name': 'GGA_X_MPW91'},
    120: {'XC_functional_name': 'GGA_X_AM05'},
    121: {'XC_functional_name': 'GGA_X_PBEA'},
    122: {'XC_functional_name': 'GGA_X_MPBE'},
    123: {'XC_functional_name': 'GGA_X_XPBE'},
    124: {'XC_functional_name': 'GGA_X_2D_B86_MGC'},
    125: {'XC_functional_name': 'GGA_X_BAYESIAN'},
    126: {'XC_functional_name': 'GGA_X_PBE_JSJR'},
    127: {'XC_functional_name': 'GGA_X_2D_B88'},
    128: {'XC_functional_name': 'GGA_X_2D_B86'},
    129: {'XC_functional_name': 'GGA_X_2D_PBE'},
    130: {'XC_functional_name': 'GGA_C_PBE'},
    131: {'XC_functional_name': 'GGA_C_LYP'},
    132: {'XC_functional_name': 'GGA_C_P86'},
    133: {'XC_functional_name': 'GGA_C_PBE_SOL'},
    134: {'XC_functional_name': 'GGA_C_PW91'},
    135: {'XC_functional_name': 'GGA_C_AM05'},
    136: {'XC_functional_name': 'GGA_C_XPBE'},
    137: {'XC_functional_name': 'GGA_C_LM'},
    138: {'XC_functional_name': 'GGA_C_PBE_JRGX'},
    139: {'XC_functional_name': 'GGA_X_OPTB88_VDW'},
    140: {'XC_functional_name': 'GGA_X_PBEK1_VDW'},
    141: {'XC_functional_name': 'GGA_X_OPTPBE_VDW'},
    142: {'XC_functional_name': 'GGA_X_RGE2'},
    143: {'XC_functional_name': 'GGA_C_RGE2'},
    144: {'XC_functional_name': 'GGA_X_RPW86'},
    145: {'XC_functional_name': 'GGA_X_KT1'},
    146: {'XC_functional_name': 'GGA_XC_KT2'},
    147: {'XC_functional_name': 'GGA_C_WL'},
    148: {'XC_functional_name': 'GGA_C_WI'},
    149: {'XC_functional_name': 'GGA_X_MB88'},
    150: {'XC_functional_name': 'GGA_X_SOGGA'},
    151: {'XC_functional_name': 'GGA_X_SOGGA11'},
    152: {'XC_functional_name': 'GGA_C_SOGGA11'},
    153: {'XC_functional_name': 'GGA_C_WI0'},
    154: {'XC_functional_name': 'GGA_XC_TH1'},
    155: {'XC_functional_name': 'GGA_XC_TH2'},
    156: {'XC_functional_name': 'GGA_XC_TH3'},
    157: {'XC_functional_name': 'GGA_XC_TH4'},
    158: {'XC_functional_name': 'GGA_X_C09X'},
    159: {'XC_functional_name': 'GGA_C_SOGGA11_X'},
    160: {'XC_functional_name': 'GGA_X_LB'},
    161: {'XC_functional_name': 'GGA_XC_HCTH_93'},
    162: {'XC_functional_name': 'GGA_XC_HCTH_120'},
    163: {'XC_functional_name': 'GGA_XC_HCTH_147'},
    164: {'XC_functional_name': 'GGA_XC_HCTH_407'},
    165: {'XC_functional_name': 'GGA_XC_EDF1'},
    166: {'XC_functional_name': 'GGA_XC_XLYP'},
    167: {'XC_functional_name': 'GGA_X_EB88'},
    168: {'XC_functional_name': 'GGA_C_PBE_MOL'},
    169: {'XC_functional_name': 'HYB_GGA_XC_PBE_MOL0'},
    170: {'XC_functional_name': 'GGA_XC_B97_D'},
    171: {'XC_functional_name': 'HYB_GGA_XC_PBE_SOL0'},
    172: {'XC_functional_name': 'HYB_GGA_XC_PBEB0'},
    173: {'XC_functional_name': 'GGA_XC_PBE1W'},
    174: {'XC_functional_name': 'GGA_XC_MPWLYP1W'},
    175: {'XC_functional_name': 'GGA_XC_PBELYP1W'},
    176: {'XC_functional_name': 'HYB_GGA_XC_PBE_MOLB0'},
    177: {'XC_functional_name': 'GGA_K_ABSP3'},
    178: {'XC_functional_name': 'GGA_K_ABSP4'},
    182: {'XC_functional_name': 'GGA_X_LBM'},
    183: {'XC_functional_name': 'GGA_X_OL2'},
    184: {'XC_functional_name': 'GGA_X_APBE'},
    185: {'XC_functional_name': 'GGA_K_APBE'},
    186: {'XC_functional_name': 'GGA_C_APBE'},
    187: {'XC_functional_name': 'GGA_K_TW1'},
    188: {'XC_functional_name': 'GGA_K_TW2'},
    189: {'XC_functional_name': 'GGA_K_TW3'},
    190: {'XC_functional_name': 'GGA_K_TW4'},
    191: {'XC_functional_name': 'GGA_X_HTBS'},
    192: {'XC_functional_name': 'GGA_X_AIRY'},
    193: {'XC_functional_name': 'GGA_X_LAG'},
    194: {'XC_functional_name': 'GGA_XC_MOHLYP'},
    195: {'XC_functional_name': 'GGA_XC_MOHLYP2'},
    196: {'XC_functional_name': 'GGA_XC_TH_FL'},
    197: {'XC_functional_name': 'GGA_XC_TH_FC'},
    198: {'XC_functional_name': 'GGA_XC_TH_FCFO'},
    199: {'XC_functional_name': 'GGA_XC_TH_FCO'},
    200: {'XC_functional_name': 'GGA_C_OPTC'},
    201: {'XC_functional_name': 'MGGA_X_LTA'},
    202: {'XC_functional_name': 'MGGA_X_TPSS'},
    203: {'XC_functional_name': 'MGGA_X_M06_L'},
    204: {'XC_functional_name': 'MGGA_X_GVT4'},
    205: {'XC_functional_name': 'MGGA_X_TAU_HCTH'},
    206: {'XC_functional_name': 'MGGA_X_BR89'},
    207: {'XC_functional_name': 'MGGA_X_BJ06'},
    208: {'XC_functional_name': 'MGGA_X_TB09'},
    209: {'XC_functional_name': 'MGGA_X_RPP09'},
    210: {'XC_functional_name': 'MGGA_X_2D_PRHG07'},
    211: {'XC_functional_name': 'MGGA_X_2D_PRHG07_PRP10'},
    212: {'XC_functional_name': 'MGGA_X_REVTPSS'},
    213: {'XC_functional_name': 'MGGA_X_PKZB'},
    214: {'XC_functional_name': 'MGGA_X_M05'},
    215: {'XC_functional_name': 'MGGA_X_M05_2X'},
    216: {'XC_functional_name': 'MGGA_X_M06_HF'},
    217: {'XC_functional_name': 'MGGA_X_M06'},
    218: {'XC_functional_name': 'MGGA_X_M06_2X'},
    219: {'XC_functional_name': 'MGGA_X_M08_HX'},
    220: {'XC_functional_name': 'MGGA_X_M08_SO'},
    221: {'XC_functional_name': 'MGGA_X_MS0'},
    222: {'XC_functional_name': 'MGGA_X_MS1'},
    223: {'XC_functional_name': 'MGGA_X_MS2'},
    224: {'XC_functional_name': 'HYB_MGGA_X_MS2H'},
    225: {'XC_functional_name': 'MGGA_X_M11'},
    226: {'XC_functional_name': 'MGGA_X_M11_L'},
    227: {'XC_functional_name': 'MGGA_X_MN12_L'},
    229: {'XC_functional_name': 'MGGA_C_CC06'},
    230: {'XC_functional_name': 'MGGA_X_MK00'},
    231: {'XC_functional_name': 'MGGA_C_TPSS'},
    232: {'XC_functional_name': 'MGGA_C_VSXC'},
    233: {'XC_functional_name': 'MGGA_C_M06_L'},
    234: {'XC_functional_name': 'MGGA_C_M06_HF'},
    235: {'XC_functional_name': 'MGGA_C_M06'},
    236: {'XC_functional_name': 'MGGA_C_M06_2X'},
    237: {'XC_functional_name': 'MGGA_C_M05'},
    238: {'XC_functional_name': 'MGGA_C_M05_2X'},
    239: {'XC_functional_name': 'MGGA_C_PKZB'},
    240: {'XC_functional_name': 'MGGA_C_BC95'},
    241: {'XC_functional_name': 'MGGA_C_REVTPSS'},
    242: {'XC_functional_name': 'MGGA_XC_TPSSLYP1W'},
    243: {'XC_functional_name': 'MGGA_X_MK00B'},
    244: {'XC_functional_name': 'MGGA_X_BLOC'},
    245: {'XC_functional_name': 'MGGA_X_MODTPSS'},
    246: {'XC_functional_name': 'GGA_C_PBELOC'},
    247: {'XC_functional_name': 'MGGA_C_TPSSLOC'},
    248: {'XC_functional_name': 'HYB_MGGA_X_MN12_SX'},
    249: {'XC_functional_name': 'MGGA_X_MBEEF'},
    250: {'XC_functional_name': 'MGGA_X_MBEEFVDW'},
    254: {'XC_functional_name': 'MGGA_XC_B97M_V'},
    255: {'XC_functional_name': 'GGA_XC_VV10'},
    257: {'XC_functional_name': 'MGGA_X_MVS'},
    258: {'XC_functional_name': 'GGA_C_PBEFE'},
    259: {'XC_functional_name': 'LDA_XC_KSDT'},
    260: {'XC_functional_name': 'MGGA_X_MN15_L'},
    261: {'XC_functional_name': 'MGGA_C_MN15_L'},
    262: {'XC_functional_name': 'GGA_C_OP_PW91'},
    263: {'XC_functional_name': 'MGGA_X_SCAN'},
    264: {'XC_functional_name': 'HYB_MGGA_X_SCAN0'},
    265: {'XC_functional_name': 'GGA_X_PBEFE'},
    266: {'XC_functional_name': 'HYB_GGA_XC_B97_1p'},
    267: {'XC_functional_name': 'MGGA_C_SCAN'},
    268: {'XC_functional_name': 'HYB_MGGA_X_MN15'},
    269: {'XC_functional_name': 'MGGA_C_MN15'},
    270: {'XC_functional_name': 'GGA_X_CAP'},
    401: {'XC_functional_name': 'HYB_GGA_XC_B3PW91'},
    402: {'XC_functional_name': 'HYB_GGA_XC_B3LYP'},
    403: {'XC_functional_name': 'HYB_GGA_XC_B3P86'},
    404: {'XC_functional_name': 'HYB_GGA_XC_O3LYP'},
    405: {'XC_functional_name': 'HYB_GGA_XC_mPW1K'},
    406: {'XC_functional_name': 'HYB_GGA_XC_PBEH'},
    407: {'XC_functional_name': 'HYB_GGA_XC_B97'},
    408: {'XC_functional_name': 'HYB_GGA_XC_B97_1'},
    410: {'XC_functional_name': 'HYB_GGA_XC_B97_2'},
    411: {'XC_functional_name': 'HYB_GGA_XC_X3LYP'},
    412: {'XC_functional_name': 'HYB_GGA_XC_B1WC'},
    413: {'XC_functional_name': 'HYB_GGA_XC_B97_K'},
    414: {'XC_functional_name': 'HYB_GGA_XC_B97_3'},
    415: {'XC_functional_name': 'HYB_GGA_XC_MPW3PW'},
    416: {'XC_functional_name': 'HYB_GGA_XC_B1LYP'},
    417: {'XC_functional_name': 'HYB_GGA_XC_B1PW91'},
    418: {'XC_functional_name': 'HYB_GGA_XC_mPW1PW'},
    419: {'XC_functional_name': 'HYB_GGA_XC_MPW3LYP'},
    420: {'XC_functional_name': 'HYB_GGA_XC_SB98_1a'},
    421: {'XC_functional_name': 'HYB_GGA_XC_SB98_1b'},
    422: {'XC_functional_name': 'HYB_GGA_XC_SB98_1c'},
    423: {'XC_functional_name': 'HYB_GGA_XC_SB98_2a'},
    424: {'XC_functional_name': 'HYB_GGA_XC_SB98_2b'},
    425: {'XC_functional_name': 'HYB_GGA_XC_SB98_2c'},
    426: {'XC_functional_name': 'HYB_GGA_X_SOGGA11_X'},
    427: {'XC_functional_name': 'HYB_GGA_XC_HSE03'},
    428: {'XC_functional_name': 'HYB_GGA_XC_HSE06'},
    429: {'XC_functional_name': 'HYB_GGA_XC_HJS_PBE'},
    430: {'XC_functional_name': 'HYB_GGA_XC_HJS_PBE_SOL'},
    431: {'XC_functional_name': 'HYB_GGA_XC_HJS_B88'},
    432: {'XC_functional_name': 'HYB_GGA_XC_HJS_B97X'},
    433: {'XC_functional_name': 'HYB_GGA_XC_CAM_B3LYP'},
    434: {'XC_functional_name': 'HYB_GGA_XC_TUNED_CAM_B3LYP'},
    435: {'XC_functional_name': 'HYB_GGA_XC_BHANDH'},
    436: {'XC_functional_name': 'HYB_GGA_XC_BHANDHLYP'},
    437: {'XC_functional_name': 'HYB_GGA_XC_MB3LYP_RC04'},
    438: {'XC_functional_name': 'HYB_MGGA_XC_M05'},
    439: {'XC_functional_name': 'HYB_MGGA_XC_M05_2X'},
    440: {'XC_functional_name': 'HYB_MGGA_XC_B88B95'},
    441: {'XC_functional_name': 'HYB_MGGA_XC_B86B95'},
    442: {'XC_functional_name': 'HYB_MGGA_XC_PW86B95'},
    443: {'XC_functional_name': 'HYB_MGGA_XC_BB1K'},
    444: {'XC_functional_name': 'HYB_MGGA_XC_M06_HF'},
    445: {'XC_functional_name': 'HYB_MGGA_XC_MPW1B95'},
    446: {'XC_functional_name': 'HYB_MGGA_XC_MPWB1K'},
    447: {'XC_functional_name': 'HYB_MGGA_XC_X1B95'},
    448: {'XC_functional_name': 'HYB_MGGA_XC_XB1K'},
    449: {'XC_functional_name': 'HYB_MGGA_XC_M06'},
    450: {'XC_functional_name': 'HYB_MGGA_XC_M06_2X'},
    451: {'XC_functional_name': 'HYB_MGGA_XC_PW6B95'},
    452: {'XC_functional_name': 'HYB_MGGA_XC_PWB6K'},
    453: {'XC_functional_name': 'HYB_GGA_XC_MPWLYP1M'},
    454: {'XC_functional_name': 'HYB_GGA_XC_REVB3LYP'},
    455: {'XC_functional_name': 'HYB_GGA_XC_CAMY_BLYP'},
    456: {'XC_functional_name': 'HYB_GGA_XC_PBE0_13'},
    457: {'XC_functional_name': 'HYB_MGGA_XC_TPSSH'},
    458: {'XC_functional_name': 'HYB_MGGA_XC_REVTPSSH'},
    459: {'XC_functional_name': 'HYB_GGA_XC_B3LYPs'},
    460: {'XC_functional_name': 'HYB_MGGA_XC_M08_HX'},
    461: {'XC_functional_name': 'HYB_MGGA_XC_M08_SO'},
    462: {'XC_functional_name': 'HYB_MGGA_XC_M11'},
    463: {'XC_functional_name': 'HYB_GGA_XC_WB97'},
    464: {'XC_functional_name': 'HYB_GGA_XC_WB97X'},
    465: {'XC_functional_name': 'HYB_GGA_XC_LRC_WPBEH'},
    466: {'XC_functional_name': 'HYB_GGA_XC_WB97X_V'},
    467: {'XC_functional_name': 'HYB_GGA_XC_LCY_PBE'},
    468: {'XC_functional_name': 'HYB_GGA_XC_LCY_BLYP'},
    469: {'XC_functional_name': 'HYB_GGA_XC_LC_VV10'},
    470: {'XC_functional_name': 'HYB_GGA_XC_CAMY_B3LYP'},
    471: {'XC_functional_name': 'HYB_GGA_XC_WB97X_D'},
    472: {'XC_functional_name': 'HYB_GGA_XC_HPBEINT'},
    473: {'XC_functional_name': 'HYB_GGA_XC_LRC_WPBE'},
    474: {'XC_functional_name': 'HYB_MGGA_X_MVSH'},
    475: {'XC_functional_name': 'HYB_GGA_XC_B3LYP5'},
    476: {'XC_functional_name': 'HYB_GGA_XC_EDF2'},
    477: {'XC_functional_name': 'HYB_GGA_XC_CAP0'},
    478: {'XC_functional_name': 'HYB_GGA_XC_LC_WPBE'},
    500: {'XC_functional_name': 'GGA_K_VW'},
    501: {'XC_functional_name': 'GGA_K_GE2'},
    502: {'XC_functional_name': 'GGA_K_GOLDEN'},
    503: {'XC_functional_name': 'GGA_K_YT65'},
    504: {'XC_functional_name': 'GGA_K_BALTIN'},
    505: {'XC_functional_name': 'GGA_K_LIEB'},
    506: {'XC_functional_name': 'GGA_K_ABSP1'},
    507: {'XC_functional_name': 'GGA_K_ABSP2'},
    508: {'XC_functional_name': 'GGA_K_GR'},
    509: {'XC_functional_name': 'GGA_K_LUDENA'},
    510: {'XC_functional_name': 'GGA_K_GP85'},
    511: {'XC_functional_name': 'GGA_K_PEARSON'},
    512: {'XC_functional_name': 'GGA_K_OL1'},
    513: {'XC_functional_name': 'GGA_K_OL2'},
    514: {'XC_functional_name': 'GGA_K_FR_B88'},
    515: {'XC_functional_name': 'GGA_K_FR_PW86'},
    516: {'XC_functional_name': 'GGA_K_DK'},
    517: {'XC_functional_name': 'GGA_K_PERDEW'},
    518: {'XC_functional_name': 'GGA_K_VSK'},
    519: {'XC_functional_name': 'GGA_K_VJKS'},
    520: {'XC_functional_name': 'GGA_K_ERNZERHOF'},
    521: {'XC_functional_name': 'GGA_K_LC94'},
    522: {'XC_functional_name': 'GGA_K_LLP'},
    523: {'XC_functional_name': 'GGA_K_THAKKAR'},
    524: {'XC_functional_name': 'GGA_X_WPBEH'},
    525: {'XC_functional_name': 'GGA_X_HJS_PBE'},
    526: {'XC_functional_name': 'GGA_X_HJS_PBE_SOL'},
    527: {'XC_functional_name': 'GGA_X_HJS_B88'},
    528: {'XC_functional_name': 'GGA_X_HJS_B97X'},
    529: {'XC_functional_name': 'GGA_X_ITYH'},
    530: {'XC_functional_name': 'GGA_X_SFAT'},
    531: {'XC_functional_name': 'HYB_MGGA_XC_WB97M_V'},
}


# TODO temporary fix for structlog unable to propagate logger
class AbinitMetainfoParser(MetainfoParser):
    @property
    def logger(self):
        return LOGGER


class MainfileParser(TextParser):
    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER

    def __init__(self):
        super().__init__()
        self.text_parser = AbinitOutParser()

    def get_workflow_method(self) -> str:
        ionmov = self.get_input_var('ionmov', 1, 0, scalar=True)
        return {
            1: 'viscous_damped_md',
            2: 'bfgs',
            3: 'bfgs',
            4: 'conjugate_gradient',
            5: 'steepest_descent',
            7: 'quenched_md',
            10: 'dic_bfgs',
            11: 'dic_bfgs',
            20: 'diis',
        }.get(ionmov)

    def get_input_var(self, name, n_dataset=1, default=None, scalar=False) -> Any:
        val = self.data_object.input_vars.get(name)
        if val is None or n_dataset > len(val) or val[n_dataset - 1] is None:
            val = [default] * n_dataset
        val = val[n_dataset - 1]
        if scalar and isinstance(val, np.ndarray | list):
            return val[-1]
        return val

    def get_datetime(self, date: str, time: str) -> datetime:
        return datetime.strptime(f'{date} {time}', '%a %d %b %Y %Hh%M')

    def get_systems(self) -> list[dict[str, Any]]:
        def get_positions() -> np.ndarray:
            natom = self.get_input_var('natom', scalar=True)
            xcart = self.get_input_var('xcart')
            if xcart is not None:
                return np.reshape(xcart, (natom, 3))
            xred = self.get_input_var('xred')
            rprim = self.get_input_var('rprim')
            if xred is not None and rprim is not None:
                xred = np.reshape(xred, (natom, 3))
                rprim = np.reshape(rprim, (natom, 3))
                return np.dot(xred, rprim.transpose())
            if natom == 1:  # handling exception
                return np.array([[0.0, 0.0, 0.0]])

        # initial system
        systems = [dict(cartesian_coordinates=get_positions())]
        for dataset in self.data_object.get('dataset', []):
            # relaxation steps
            systems.extend(dataset.get('relaxation', []))
        return systems

    def get_energy_contributions(self, source: dict[str, Any]) -> list[dict[str, Any]]:
        return [
            dict(value=val, name=key)
            for key, val in source.items()
            if key.startswith('energy_') and key != 'energy_total'
        ]

    def get_outputs(self) -> list[dict[str, Any]]:
        outputs = []
        for dataset in self.data_object.get('dataset', []):
            results = dict(dataset.get('results') or {})
            scf_steps = self.get_scf_steps(dataset)
            if scf_steps:
                results['scf_steps'] = scf_steps
            outputs.append(results)
            # relaxation steps
            for step in dataset.get('relaxation', []):
                step_output = dict(step)
                scf_steps = self.get_scf_steps(step)
                if scf_steps:
                    step_output['scf_steps'] = scf_steps
                outputs.append(step_output)
        return outputs

    def get_atoms(self) -> list[dict[str, Any]]:
        znucl = self.get_input_var('znucl', default=1)
        typat = self.get_input_var('typat', default=1)
        if znucl is None or typat is None:
            return []
        return [dict(label=chemical_symbols[int(znucl[n_at - 1])]) for n_at in typat]

    def get_periodic_boundary_conditions(
        self, lattice_vectors: Any = None
    ) -> list[bool] | None:
        if lattice_vectors is None:
            return None
        return [vec is not None for vec in lattice_vectors]

    def get_xc_functionals(self) -> list[dict[str, Any]]:
        ixc = self.get_input_var('ixc', 1, 1, scalar=True)
        if ixc >= 0:
            return [
                {'unidentified': True}
                if functional.get('XC_functional_name') == '?'
                else functional
                for functional in ABINIT_NATIVE_IXC.get(ixc, [])
            ]
        # LibXC path: the negative value packs two LibXC ids positionally. Id 0
        # means the slot carries no functional; a non-zero id absent from the
        # table is handed to the schema as a raw LibXC id to resolve or mark.
        functional1 = -ixc // 1000
        functional2 = -ixc - functional1 * 1000
        components = []
        for functional_id in (functional1, functional2):
            if functional_id == 0:
                continue
            mapped = ABINIT_LIBXC_IXC.get(functional_id)
            components.append(mapped if mapped else {'libxc_id': functional_id})
        return components

    def get_bandstructures(
        self, eigenvalues: np.ndarray, occupations: np.ndarray
    ) -> list[dict[str, Any]]:
        n_spin_channels = 2
        nsppol = self.get_input_var('nsppol', 2, 1, scalar=True)
        eigs = np.reshape(
            eigenvalues,
            (
                nsppol,
                len(eigenvalues) // nsppol,
                np.size(eigenvalues) // len(eigenvalues),
            ),
        )

        kpts = eigs.T[3:6].T[0]
        # if len(kpts) == 1:  # no bs for one kpoint (atoms or molecules)
        #     return []

        nband = int(eigs.T[1].T[0][0])
        eigs = eigs.T[6 : 6 + nband].T
        is_spin_polarized = nsppol == n_spin_channels
        bandstructures = []
        for n, eig in enumerate(eigs):
            entry = dict(energies=eig, k_points=kpts)
            if is_spin_polarized:
                entry['spin_channel'] = n
            bandstructures.append(entry)

        if occupations is not None:
            occs = np.reshape(
                occupations,
                (
                    nsppol,
                    len(occupations) // nsppol,
                    np.size(occupations) // len(occupations),
                ),
            )
            if np.shape(eigs) != np.shape(occs):
                self.logger.error('Inconsistent shape of eigenvalues and occupations')
            for n, occ in enumerate(occs):
                bandstructures[n]['occupations'] = occ

        return bandstructures

    def get_geometry_convergence(self):
        # TODO consider to do this for each dataset and not only the first one
        tolmxde = self.get_input_var('tolmxde', n_dataset=1, default=0.0, scalar=True)
        tolmxf = self.get_input_var(
            name='tolmxf', n_dataset=1, default=0.0, scalar=True
        )

        return [
            EnergyConvergenceTarget(
                threshold=tolmxde * ureg.hartree,
                threshold_type='absolute',
            ),
            ForceConvergenceTarget(
                threshold=tolmxf * ureg.hartree / ureg.bohr,
                threshold_type='maximum',
            ),
        ]

    def get_single_point_convergence(self) -> list[EnergyConvergenceTarget]:
        toldfe = self.get_input_var('toldfe', n_dataset=1, default=None, scalar=True)
        if toldfe is None:
            return []
        return [
            EnergyConvergenceTarget(
                threshold=toldfe,
                threshold_type='absolute',
            )
        ]

    def get_scf_steps(self, source: dict[str, Any]) -> dict[str, Any]:
        min_columns_for_delta = 2
        self_consistent = source.get('self_consistent', {})
        iterations = self_consistent.get('energy_total_scf_iteration', [])
        if not iterations:
            return {}

        energies_total = []
        delta_energies_total = []
        extra_columns = {}

        for step in iterations:
            row = np.array(step, dtype=float).reshape(-1)
            if row.size == 0:
                continue
            energies_total.append(row[0] * ureg.hartree)
            if row.size >= min_columns_for_delta:
                delta_energies_total.append(abs(row[1]) * ureg.hartree)
            for col_idx, val in enumerate(row[2:], start=3):
                extra_columns.setdefault(f'column_{col_idx}', []).append(float(val))

        scf_steps = {'energies_total': energies_total}
        if delta_energies_total:
            scf_steps['delta_energies_total'] = delta_energies_total
        if extra_columns:
            scf_steps['code_specific_quantities'] = extra_columns
        return scf_steps

    def get_band_gaps(
        self, eigenvalues: np.ndarray, occupations: np.ndarray
    ) -> list[dict[str, Any]]:
        """Calculate band gaps from eigenvalues and occupations using common utility."""
        if eigenvalues is None or occupations is None:
            return []

        bandstructures = self.get_bandstructures(eigenvalues, occupations)
        gaps = []
        for bandstructure in bandstructures:
            eigs = bandstructure.get('energies')
            occs = bandstructure.get('occupations')
            spin_channel = bandstructure.get('spin_channel')

            # Use common utility for band gap calculation
            gap_result = calculate_band_gap_from_occupations(
                eigs, occs, spin_channel=spin_channel
            )
            if gap_result is not None:
                gaps.append(gap_result)

        return gaps


class DosParser(TextParser):
    # TODO temporary fix for structlog unable to propagate logger
    @property
    def logger(self):
        return LOGGER

    def __init__(self):
        super().__init__()
        self.text_parser = DataTextParser()

    def get_dos(self, source: np.ndarray) -> list[dict[str, Any]]:
        nsp = self.data.get('nspinpol')
        dos = []
        for dos_sp in np.reshape(
            source, (nsp, len(source) // nsp, np.size(source) // len(source))
        ):
            dos_sp_t = dos_sp.T
            dos.append(dict(energies=dos_sp_t[0], value=dos_sp_t[1]))
        return dos


class AbinitArchiveWriter(ArchiveWriter):
    code_name = 'ABINIT'
    annotation_key = abinit.OUT_KEY

    def __init__(self):
        super().__init__()
        self.mainfile_parser = MainfileParser()
        self.metainfo_parser = AbinitMetainfoParser()
        self.dos_parser = DosParser()

    def parse_workflow(self):
        ionmov = self.mainfile_parser.get_input_var('ionmov', 1, [0])[0]
        vis = self.mainfile_parser.get_input_var('vis', 1, [100.0])[0]
        convergence = None
        if ionmov in [2, 3, 4, 5, 7, 10, 11, 20] or (ionmov == 1 and vis > 0.0):
            workflow = GeometryOptimization()
            workflow.method = GeometryOptimizationMethod()

            # Set optimization type based on optcell parameter
            optcell = self.mainfile_parser.get_input_var('optcell', 1, [0])[0]
            if optcell == 0:
                workflow.method.optimization_type = 'atomic'
            elif optcell == 1:
                workflow.method.optimization_type = 'cell_volume'
            else:
                workflow.method.optimization_type = 'cell_shape'

            convergence = self.mainfile_parser.get_geometry_convergence()
        elif ionmov in [6, 8, 9, 12, 13, 14, 23] or (ionmov == 1 and vis == 0.0):
            workflow = MolecularDynamics()
        else:
            workflow = SinglePoint()
            workflow.method = SinglePointMethod()
            convergence = self.mainfile_parser.get_single_point_convergence()
        self.archive.workflow2 = workflow
        self.metainfo_parser.annotation_key = self.annotation_key
        self.metainfo_parser.data_object = self.archive.workflow2
        self.mainfile_parser.convert(self.metainfo_parser)
        # Assign convergence targets only after convert() to preserve the
        # polymorphic EnergyConvergenceTarget/ForceConvergenceTarget subclasses;
        # see the `add_mapping_annotation` docstring for why the ordering matters.
        if convergence:
            self.archive.workflow2.method.convergence_targets = convergence

    def write_to_archive(self):
        self.archive.data = Simulation(program=Program(name=self.code_name))
        self.metainfo_parser.annotation_key = self.annotation_key
        self.metainfo_parser.data_object = self.archive.data

        self.mainfile_parser.filepath = self.mainfile
        self.mainfile_parser.convert(self.metainfo_parser)

        # parse dos from dos file
        self.metainfo_parser.annotation_key = abinit.DOS_KEY
        # DS2_DOS files
        file_root = self.mainfile_parser.data_object.get('x_abinit_output_files_root')
        if file_root is None:
            file_root = f'{os.path.basename(self.mainfile).rstrip(".out")}_o'
        self.dos_parser.filepath = os.path.join(
            os.path.dirname(self.mainfile), f'{file_root}_DS2_DOS'
        )
        # read nspin from mainfile parser
        self.dos_parser.data.setdefault(
            'nspinpol', self.mainfile_parser.get_input_var('nsppol', 2, 1, scalar=True)
        )
        self.dos_parser.convert(self.metainfo_parser, update_mode='merge@-1')

        self.parse_workflow()

        gw_archive = self.child_archives.get('GW')
        if gw_archive is not None:
            gw_archive.data = Simulation(program=Program(name=self.code_name))

            writer = AbinitArchiveWriter()
            writer.annotation_key = 'gw_out'
            writer.write(self.mainfile, gw_archive, self.logger)

            workflow_archive = self.child_archives['GW_workflow']
            workflow_archive.workflow2 = DFTGWWorkflow(
                tasks=[self.archive.workflow2, gw_archive.workflow2]
            )

        self.metainfo_parser.close()
        self.mainfile_parser.close()
        self.dos_parser.close()


class AbinitParser(MatchingParser):
    """
    Main parser interface to NOMAD.
    """

    def is_mainfile(
        self,
        filename: str,
        mime,
        buffer: bytes,
        decoded_buffer: str,
        compression: str = None,
    ) -> bool | Iterable:
        is_mainfile = super().is_mainfile(
            filename, mime, buffer, decoded_buffer, compression
        )
        if is_mainfile:
            out_parser = AbinitOutParser()
            out_parser.findall = False
            out_parser.mainfile = filename
            ds_numbers = out_parser.dataset_numbers
            optdriver = out_parser.input_vars.get('optdriver', [])
            out_parser.findall = True
            n_gw = [4, 66]
            if n_gw[0] in ds_numbers and (1 and 2 and 3) not in ds_numbers:
                return True
            if len(optdriver) == n_gw[0] and (optdriver[-1] in n_gw):
                self.creates_children = True
                return ['GW', 'GW_workflow']
            return True

    def parse(
        self,
        mainfile: str,
        archive: EntryArchive,
        logger: BoundLogger = None,
        child_archives: dict[str, EntryArchive] = {},
    ):
        archive_writer = AbinitArchiveWriter()
        archive_writer.write(mainfile, archive, logger, child_archives)
     
        # run the old parser
        # TODO remove
        from electronicparsers.abinit.parser import AbinitParser  # noqa

        AbinitParser().parse(mainfile, archive, logger)
