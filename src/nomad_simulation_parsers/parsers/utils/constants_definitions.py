import numpy as np
from ase import data as asedata

RE_FLOAT = r'[-+]?\d+\.*\d*(?:[Ee][-+]\d+)?'
RE_N = r'[\n\r]'
MOLE = 6.022140857e23  # Avogadro number
REFERENCE_MASSES = np.array(asedata.atomic_masses)
CHEMICAL_SYMBOLS = asedata.chemical_symbols
