import re

import numpy as np
from nomad_file_parser import Quantity, TextParser

from ..common import general_quantities, header_quantities


class GIPAWFileParser(TextParser):
    def __init__(self):
        super().__init__(None)

    def init_quantities(self):
        re_float = r' *[-+]?\d+\.\d*(?:[Ee][-+]\d+)? *'

        def str_to_ms_data_list(val_in):
            pattern = re.compile(
                r'Atom\s+(\d+)\s+(\w+).*?\n'
                r'\s*([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\n'
                r'\s*([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\n'
                r'\s*([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)',
                re.DOTALL,
            )

            data = []
            for match in pattern.findall(val_in):
                atom_num = int(match[0])
                atom_type = match[1]
                values = [float(x) for x in match[2:]]
                data.append([atom_type, atom_num] + values)
            return data

        def str_to_chi_tensor(val_in):
            lines = val_in.strip().splitlines()
            tensor = []
            for line in lines:
                if line.strip():
                    row = [float(x) for x in line.strip().split()]
                    tensor.append(row)

            # Convert from units of 10^{-6} cm^3/mol to m^3/mol
            FACTOR = 1e-12
            res = np.array(tensor) * FACTOR

            return res

        def parse_tensor_block(val_in: str):
            min_length = 3
            lines = [
                line.strip() for line in val_in.strip().splitlines() if line.strip()
            ]
            result = []
            for i in range(0, len(lines), 3):
                block = lines[i : i + 3]
                if len(block) < min_length:
                    continue

                values = []
                atom_type = None
                atom_index = None

                for row in block:
                    parts = row.split()
                    if atom_type is None:
                        atom_type = parts[0]
                        atom_index = int(parts[1])
                    values.extend([float(p) for p in parts[2:]])

                result.append([atom_type, atom_index] + values)
            return result

        def parse_scalar_block(val):
            min_length = 6
            lines = [line.strip() for line in val.strip().splitlines() if line.strip()]
            results = []
            for line in lines:
                parts = line.split()
                if len(parts) != min_length:
                    continue
                results.append([parts[0], int(parts[1]), float(parts[-1])])
            return results

        def str_to_gtensor(val_in):
            lines = val_in.strip().splitlines()
            tensor = []
            for line in lines:
                if line.strip():
                    row = [float(x) for x in line.strip().split()]
                    tensor.append(row)
            return tensor

        self._quantities = [
            Quantity(
                'header',
                r'([Pp]rogram GIPAW[\s\S]+?)GIPAW job',
                repeats=False,
                sub_parser=TextParser(
                    quantities=header_quantities + general_quantities
                ),
            ),
            Quantity(
                'ms_list',
                r'Total NMR chemical shifts in ppm:\s*((?:.*?\n)*?)\s*Initialization:',
                str_operation=str_to_ms_data_list,
                convert=False,
            ),
            Quantity(
                'chi_bare_pGv',
                rf'chi_bare\s+pGv\s+\(\w+\)\s+in\s+10\^{{-6}}\s+cm\^3/mol:\s*\n'
                rf'((?:\s*{re_float}\s+{re_float}\s+{re_float}\s*\n?){{1,}})',
                repeats=False,
                str_operation=str_to_chi_tensor,
                convert=False,
            ),
            Quantity(
                'chi_bare_vGv',
                rf'chi_bare\s+vGv\s+\(\w+\)\s+in\s+10\^{{-6}}\s+cm\^3/mol:\s*\n'
                rf'((?:\s*{re_float}\s+{re_float}\s+{re_float}\s*\n?){{1,}})',
                repeats=False,
                str_operation=str_to_chi_tensor,
                convert=False,
            ),
            Quantity(
                'efg',
                r'----- total EFG \(symmetrized\) -----\n((?:.*?\n)*?)\s+'
                r'NQR/NMR SPECTROSCOPIC PARAMETERS:',
                str_operation=parse_tensor_block,
                convert=False,
            ),
            Quantity(
                'hyperfine_dipolar',
                r'----- total dipolar -----\n((?:.*?\n)*?)\s+----- '
                r'total dipolar \(symmetrized\) -----',
                str_operation=parse_tensor_block,
                convert=False,
            ),
            Quantity(
                'hyperfine_fermi_contact',
                r'----- Fermi contact in G -----\n((?:.*?\n)*?)\s+Initialization:',
                str_operation=parse_scalar_block,
                convert=False,
            ),
            Quantity(
                'delta_g_total_paratec',
                rf'Delta_g total \(SOO a la Paratec\):\s*-+\s*\n'
                rf'((?:\s*{re_float}\s+{re_float}\s+{re_float}\s*\n){{3}})',
                repeats=False,
                str_operation=str_to_gtensor,
                convert=False,
            ),
            Quantity(
                'delta_g_total',
                rf'Delta_g total \(SOO as in Eq\.\(7\)\):\s*-+\s*\n'
                rf'((?:\s*{re_float}\s+{re_float}\s+{re_float}\s*\n){{3}})',
                repeats=False,
                str_operation=str_to_gtensor,
                convert=False,
            ),
        ]
