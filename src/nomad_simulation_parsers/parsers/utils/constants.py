import numpy as np
from ase import data as asedata

RE_FLOAT = r'[-+]?\d+\.*\d*(?:[Ee][-+]\d+)?'
RE_N = r'[\n\r]'
# TODO: replace with (1 * ureg.avogadro_number).to_base_units()
MOLE = 6.022140857e23  # Avogadro number was updated to 6.02214076e23 in 2019!
REFERENCE_MASSES = np.array(asedata.atomic_masses)
CHEMICAL_SYMBOLS = asedata.chemical_symbols
