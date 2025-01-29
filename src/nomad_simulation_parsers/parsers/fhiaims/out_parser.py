from datetime import datetime
from typing import Any, Optional, Union

import numpy as np
import pint
from nomad.parsing.file_parser import Quantity, TextParser
from nomad.units import ureg
from nomad.utils import get_logger

from nomad_simulation_parsers.parsers.utils.general import log

RE_FLOAT = r'[-+]?\d+\.\d*(?:[Ee][-+]\d+)?'
RE_N = r'[\n\r]'
RE_GW_FLAG = rf'{RE_N}\s*(?:qpe_calc|sc_self_energy)\s*([\w]+)'

LOGGER = get_logger(__name__)

units_mapping = {'Ha': ureg.hartree, 'eV': ureg.eV}


# TODO move this to general utils
def str_to_unit(val_in: str) -> Optional[pint.Unit]:
    """
    Convert string to pint.Unit
    """
    val = val_in.strip().lower()
    unit = None
    if val.startswith('a'):
        unit = 1 / ureg.angstrom
    elif val.startswith('b'):
        unit = 1 / ureg.bohr
    return unit


# TODO merge with str_to_scf_convergence
def str_to_energy_components(val_in: str) -> dict[str, Union[float, pint.Quantity]]:
    """
    Parse key, value pairs of energy components from string
    """
    val = [v.strip() for v in val_in.strip().splitlines()]
    res = dict()
    min_n = 2
    for val_n in val:
        v = val_n.lstrip(' |').strip().split(':')
        if len(v) < min_n or not v[1]:
            continue
        vi = v[1].split()
        if not vi[0][-1].isdecimal() or len(vi) < min_n:
            continue
        unit = units_mapping.get(vi[1], None)
        res[v[0].strip()] = float(vi[0]) * unit if unit is not None else float(vi[0])
    return res


def str_to_scf_convergence(val_in: str) -> dict[str, Union[float, pint.Quantity]]:
    """
    Parse key, value pairs of scf convergence parameters.
    """
    res = dict()
    min_n = 2
    for val_n in val_in.strip().splitlines():
        v = val_n.lstrip(' |').split(':')
        if len(v) != min_n:
            break
        vs = v[1].split()
        unit = None
        if len(vs) > 1:
            unit = units_mapping.get(vs[1], None)
        res[v[0].strip()] = float(vs[0]) * unit if unit is not None else float(vs[0])
    return res


def str_to_scf_convergence2(val_in: str) -> dict[str, pint.Quantity]:
    """
    Parse scf energy.
    """
    n = 7
    val = val_in.split('|')
    if len(val) != n:
        return {}
    energy = float(val[3]) * ureg.eV
    return {'Change of total energy': energy}


def str_to_atomic_forces(val_in: str) -> pint.Quantity:
    """
    Parse atomic forces data from string
    """
    val = [v.lstrip(' |').split() for v in val_in.strip().splitlines()]
    n = 4
    forces = np.array([v[1:4] for v in val if len(v) == n], dtype=float)
    return forces * ureg.eV / ureg.angstrom


def str_to_dos_files(val_in: str) -> tuple[str, list[str]]:
    """
    Parse dos files and species data from string.
    """
    val = [v.strip() for v in val_in.strip().splitlines()]
    files = []
    species = []
    for v in val[1:]:
        if v.startswith('| writing') and 'raw data' in v:
            files.append(v.split('to file')[1].strip(' .'))
            if 'for species' in v:
                species.append(v.split('for species')[1].split()[0])
        elif not v.startswith('|'):
            break
    return files, list(set(species))


def str_to_array_size_parameters(val_in: str) -> dict[str, int]:
    """
    Parse array parameters from string.
    """
    n = 2
    val = [v.lstrip(' |').split(':') for v in val_in.strip().splitlines()]
    return {v[0].strip(): int(v[1]) for v in val if len(v) == n}


# TODO merge with str_to_species
def str_to_species_in(val_in: str) -> list[dict[str, Any]]:
    """
    Parse the input species data from string.
    """
    val = [v.strip() for v in val_in.splitlines()]
    data = []
    min_n = 2
    species = dict()
    for i in range(len(val)):
        if val[i].startswith('Reading configuration options for species'):
            if species:
                data.append(species)
            species = dict(species=val[i].split('species')[1].split()[0])
        elif not val[i].startswith('| Found'):
            continue
        val[i] = val[i].split(':')
        if len(val[i]) == 1:
            val[i] = val[i][0].split('treatment for')
        if len(val[i]) < min_n:
            continue
        k = val[i][0].split('Found')[1].strip()
        v = val[i][1].replace(',', '').split()
        if 'Gaussian basis function' in k and 'elementary' in v:
            n_gaussians = int(v[v.index('elementary') - 1])
            for j in range(n_gaussians):
                v.extend(val[i + j + 1].lstrip('|').split())
        v = v[0] if len(v) == 1 else v
        if val[i][0] in species:
            species[k].extend([v])
        else:
            species[k] = [v]
    data.append(species)
    return data


def str_to_species(val_in: str) -> dict[str, list[str]]:
    """
    Parse species data from string.
    """
    data = dict()
    min_n = 2
    val = [v.strip() for v in val_in.splitlines()]
    for i in range(len(val)):
        if val[i].startswith('species'):
            data['species'] = val[i].split()[1]
        elif not val[i].startswith('| Found'):
            continue
        val[i] = val[i].split(':')
        if len(val[i]) == 1:
            val[i] = val[i][0].split('treatment for')
        if len(val[i]) < min_n:
            continue
        k = val[i][0].split('Found')[1].strip()
        v = val[i][1].replace(',', '').split()
        if 'Gaussian basis function' in k and 'elementary' in v:
            n_gaussians = int(v[v.index('elementary') - 1])
            for j in range(n_gaussians):
                v.extend(val[i + j + 1].lstrip('|').split())
        v = v[0] if len(v) == 1 else v
        if k in data:
            data[k].extend([v])
        else:
            data[k] = [v]
    return data


def str_to_hirshfeld(val_in: str) -> dict[str, list[float]]:
    """
    Parse Hirshfeld parameters
    """
    val = [v.strip() for v in val_in.strip().splitlines()]
    data = dict(atom=val[0])
    for val_n in val[1:]:
        if val_n.startswith('|'):
            v = val_n.strip(' |').split(':')
            if v[0][0].isalpha():
                key = v[0].strip()
                data[key] = []
            data[key].extend([float(vi) for vi in v[-1].split()])
    return data


def str_to_frequency(val_in: str) -> tuple[int, float]:
    """
    Parse frequency value
    """
    val = val_in.strip().split()
    return int(val[0]), float(val[1])


def str_to_gw_eigs(val_in: str) -> dict[str, np.ndarray]:
    """
    Parse GW eigenvalues
    """
    val = [v.split() for v in val_in.splitlines()]
    keys = val[0]
    data = []
    for v in val[1:]:
        if len(keys) == len(v) and v[0].isdecimal():
            data.append(v)
    data = np.array(data, dtype=float)
    data = np.transpose(data)
    res = {keys[i]: data[i] for i in range(len(data))}
    return res


def str_to_gw_scf(val_in: str) -> dict[str, Union[float, pint.Quantity]]:
    """
    Parse GW results
    """
    val = [v.split(':') for v in val_in.splitlines()]
    data = {}
    n = 2
    for v in val:
        if len(v) == n:
            data[v[0].strip(' |')] = float(v[1].split()[0]) * ureg.eV
        if 'Fit accuracy for G' in v[0]:
            data['Fit accuracy for G(w)'] = float(v[0].split()[-1])
    return data


def str_to_md_calculation_info(val_in: str) -> dict[str, Union[float, pint.Quantity]]:
    """
    Parse molecular dynamics parameters
    """
    val = [v.strip() for v in val_in.strip().splitlines()]
    min_n = 2
    res = dict()
    for val_n in val:
        v = val_n.lstrip(' |').strip().split(':')
        if len(v) < min_n or not v[1]:
            continue
        vi = v[1].split()
        if not vi[0][-1].isdecimal():
            continue
        elif len(vi) < min_n:
            res[v[0].strip()] = float(vi[0])
        else:
            unit = units_mapping.get(vi[1], None)
            res[v[0].strip()] = (
                float(vi[0]) * unit if unit is not None else float(vi[0])
            )
    return res


def str_to_quantity(val_in: str) -> Union[float, pint.Quantity, None]:
    val = val_in.split()
    n = 2
    if len(val) == 1:
        return float(val[0])
    elif len(val) == n:
        return float(val[0]) * ureg(val[1])
    else:
        return None


@log
def str_to_ureg(val_in, logger=LOGGER):
    return ureg(val_in.replace('^', '**'))


def str_to_md_control_in(val_in):
    val = val_in.split()
    return {val[0]: ' '.join(val[1:])}


class FHIAimsOutFileParser(TextParser):
    def init_quantities(self):
        structure_quantities = [
            Quantity(
                'labels',
                rf'(?:Species\s*([A-Z][a-z]*)|([A-Z][a-z]*)\w*{RE_N})',
                repeats=True,
            ),
            Quantity(
                'positions',
                rf'({RE_FLOAT})\s+({RE_FLOAT})\s+({RE_FLOAT}) *{RE_N}',
                dtype=np.dtype(np.float64),
                repeats=True,
            ),
            Quantity(
                'positions',
                rf'atom +({RE_FLOAT})\s+({RE_FLOAT})\s+({RE_FLOAT})',
                dtype=np.dtype(np.float64),
                repeats=True,
            ),
            Quantity(
                'velocities',
                rf'velocity\s+({RE_FLOAT})\s+({RE_FLOAT})\s+({RE_FLOAT})',
                dtype=np.dtype(np.float64),
                repeats=True,
            ),
        ]

        eigenvalues = Quantity(
            'eigenvalues',
            r'Writing Kohn\-Sham eigenvalues\.([\s\S]+?'
            rf'State[\s\S]+?)(?:{RE_N}{RE_N} +[A-RT-Z])',
            repeats=True,
            sub_parser=TextParser(
                quantities=[
                    Quantity(
                        'kpoints',
                        rf'{RE_N} *K-point:\s*\d+ at\s*'
                        rf'({RE_FLOAT})\s*({RE_FLOAT})\s*({RE_FLOAT})',
                        dtype=float,
                        repeats=True,
                    ),
                    Quantity(
                        'occupation_eigenvalue',
                        rf'{RE_N} *\d+\s*({RE_FLOAT})\s*({RE_FLOAT})\s*{RE_FLOAT}',
                        repeats=True,
                    ),
                ]
            ),
        )

        date_time = Quantity(
            'date_time',
            r'Date\s*:\s*(\d+), Time\s*:\s*([\d\.]+)\s*',
            repeats=False,
            convert=False,
            str_operation=lambda x: datetime.strptime(
                x, '%Y%m%d %H%M%S.%f'
            ).timestamp(),
        )

        scf_quantities = [
            # TODO add section_eigenvalues to scf_iteration
            date_time,
            eigenvalues,
            Quantity(
                'energy_components',
                rf'{RE_N} *Total energy components:([\s\S]+?)((?:{RE_N}{RE_N}|\| '
                rf'Electronic free energy per atom\s*:\s*[Ee\d\.\-]+ eV))',
                repeats=False,
                str_operation=str_to_energy_components,
                convert=False,
            ),
            Quantity(
                'forces',
                rf'{RE_N} *Total forces\([\s\d]+\)\s*:([\s\d\.\-\+Ee]+){RE_N}',
                repeats=True,
            ),
            Quantity(
                'stress_tensor',
                rf'{RE_N} *Sum of all contributions\s*:\s*([\d\.\-\+Ee ]+{RE_N})',
                repeats=False,
            ),
            Quantity('pressure', r' *\|\s*Pressure:\s*([\d\.\-\+Ee ]+)', repeats=False),
            Quantity(
                'scf_convergence',
                rf'{RE_N} *Self-consistency convergence accuracy:([\s\S]+?)(\| '
                rf'Change of total energy\s*:\s*[\d\.\-\+Ee V]+)',
                repeats=False,
                str_operation=str_to_scf_convergence,
                convert=False,
            ),
            Quantity(
                'humo',
                r'Highest occupied state \(VBM\) '
                r'at\s*([\d\.\-\+Ee ]+) (?P<__unit>\w+)',
                repeats=False,
                dtype=float,
            ),
            Quantity(
                'lumo',
                r'Lowest unoccupied state \(CBM\) '
                r'at\s*([\d\.\-\+Ee ]+) (?P<__unit>\w+)',
                repeats=False,
                dtype=float,
            ),
            Quantity(
                'fermi_level',  # older version
                rf'{RE_N} *\| Chemical potential \(Fermi level\) in '
                rf'(\w+)\s*:([\d\.\-\+Ee ]+)',
                str_operation=lambda x: float(x.split()[1])
                * units_mapping.get(x.split()[0]),
            ),
            Quantity(
                'fermi_level',  # newer version
                rf'{RE_N} *\| Chemical potential \(Fermi level\)\:'
                rf'\s*([\-\d\.]+)\s*(\w+)',
                str_operation=lambda x: float(x.split()[0])
                * units_mapping.get(x.split()[1], 1),
            ),
            Quantity(
                'time_calculation',
                r'Time for this iteration +: +[\d\.]+ s +([\d\.]+) s',
                dtype=float,
            ),
        ]

        calculation_quantities = [
            Quantity(
                'self_consistency',
                r'Begin self\-consistency iteration #\s*\d+([\s\S]+?'
                'Total energy evaluation[s:\d\. ]+)',
                repeats=True,
                sub_parser=TextParser(quantities=scf_quantities),
            ),
            # different format for scf loop
            Quantity(
                'self_consistency',
                rf'{RE_N} *SCF\s*\d+\s*:([ \|\-\+Ee\d\.s]+)',
                repeats=True,
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'scf_convergence',
                            r'([\s\S]+)',
                            str_operation=str_to_scf_convergence2,
                            repeats=False,
                            convert=False,
                        )
                    ]
                ),
            ),
            Quantity(
                'structure',
                rf'Atomic structure(.|\n)*\| *Atom *x \[A\] *y \[A\] *z \[A\]([\s\S]+?'
                rf'Species[\s\S]+?(?:{RE_N} *{RE_N}| 1\: ))',
                repeats=False,
                convert=False,
                sub_parser=TextParser(quantities=structure_quantities),
            ),
            Quantity(
                'structure',
                rf'{RE_N} *(atom +{RE_FLOAT}[\s\S]+?(?:{RE_N} *{RE_N}|\-\-\-))',
                repeats=False,
                convert=False,
                sub_parser=TextParser(quantities=structure_quantities),
            ),
            Quantity(  # This quantity is double defined in self._quantities
                'lattice_vectors',
                rf'{RE_N} *lattice_vector([\d\.\- ]+){RE_N} *lattice_vector([\d\.\- ]+)'
                rf'{RE_N} *lattice_vector([\d\.\- ]+)',
                unit='angstrom',
                repeats=False,
                shape=(3, 3),
                dtype=float,
            ),
            Quantity(
                'energy',
                rf'{RE_N} *Energy and forces in a compact form:([\s\S]+?(?:{RE_N}{RE_N}'
                '|Electronic free energy\s*:\s*[\d\.\-Ee]+ eV))',
                str_operation=str_to_energy_components,
                repeats=False,
                convert=False,
            ),
            # in some cases, the energy components are also printed for after a
            # calculation same format as in scf iteration, they are printed also in
            # initialization so we should get last occurence
            Quantity(
                'energy_components',
                rf'{RE_N} *Total energy components:([\s\S]+?)((?:{RE_N}{RE_N}|\| '
                rf'Electronic free energy per atom\s*:\s*[\d\.\-Ee]+ eV))',
                repeats=True,
                str_operation=str_to_energy_components,
                convert=False,
            ),
            Quantity(
                'energy_xc',
                rf'{RE_N} *Start decomposition of the XC Energy([\s\S]+?)'
                rf'End decomposition of the XC Energy',
                str_operation=str_to_energy_components,
                repeats=False,
                convert=False,
            ),
            eigenvalues,
            Quantity(
                'forces',
                rf'{RE_N} *Total atomic forces.*?\[eV/Ang\]:\s*([\d\.Ee\-\+\s\|]+)',
                str_operation=str_to_atomic_forces,
                repeats=False,
                convert=False,
            ),
            # TODO no metainfo for scf forces but old parser put it in
            # atom_forces_free_raw
            Quantity(
                'forces_raw',
                rf'{RE_N} *Total forces\([\s\d]+\)\s*:([\s\d\.\-\+Ee]+){RE_N}',
                repeats=True,
                dtype=float,
            ),
            Quantity(
                'force_maximum',
                rf'Maximum force component is +({RE_FLOAT}) *eV/A',
                dtype=float,
                unit='eV/angstrom',
            ),
            Quantity(
                'time_calculation',
                rf'{RE_N} *\| Time for this force evaluation\s*:'
                rf'\s*[\d\.]+ s\s*([\d\.]+) s',
                repeats=False,
                dtype=float,
            ),
            Quantity(
                'total_dos_files',
                r'Calculating total density of states([\s\S]+?)\-{5}',
                str_operation=str_to_dos_files,
                repeats=False,
                convert=False,
            ),
            Quantity(
                'atom_projected_dos_files',
                r'Calculating atom\-projected density of states([\s\S]+?)\-{5}',
                str_operation=str_to_dos_files,
                repeats=False,
                convert=False,
            ),
            Quantity(
                'species_projected_dos_files',
                r'Calculating angular momentum projected density of states'
                r'([\s\S]+?)\-{5}',
                str_operation=str_to_dos_files,
                repeats=False,
                convert=False,
            ),
            Quantity(
                'vdW_TS',
                rf'(Evaluating non\-empirical van der Waals correction[\s\S]+?)'
                rf'(?:\|\s*Converged\.|\-{5}{RE_N}{RE_N})',
                repeats=False,
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'kind',
                            r'Evaluating non\-empirical van der Waals correction '
                            '\(([\w /]+)\)',
                            repeats=False,
                            convert=False,
                            flatten=False,
                        ),
                        Quantity(
                            'atom_hirshfeld',
                            r'\| Atom\s*\d+:([\s\S]+?)\-{5}',
                            str_operation=str_to_hirshfeld,
                            repeats=True,
                            convert=False,
                        ),
                    ]
                ),
            ),
            Quantity(
                'converged',
                r'Self\-consistency cycle (converged)\.',
                repeats=False,
                dtype=str,
            ),
            date_time,
        ]

        molecular_dynamics_quantities = [
            Quantity(
                'md_run',
                r' *Running\s*Born-Oppenheimer\s*molecular\s*dynamics\s*in\s*([A-Z]{3})'
                r'\s*ensemble*\D*\s*with\s*([A-Za-z\-]*)\s*thermostat',
                repeats=False,
                convert=False,
            ),
            Quantity(
                'md_timestep',
                rf'{RE_N} *Molecular dynamics time step\s*=\s*'
                rf'({RE_FLOAT} [A-Za-z]*)\s*{RE_N}',
                str_operation=str_to_quantity,
                repeats=False,
                convert=False,
            ),
            Quantity(
                'md_simulation_time',
                rf'{RE_N} *\| *simulation time\s*=\s*({RE_FLOAT} [A-Za-z]*)\s*{RE_N}',
                str_operation=str_to_quantity,
                repeats=False,
                convert=False,
            ),
            Quantity(
                'md_temperature',
                rf'{RE_N} *\| *at temperature\s*=\s*({RE_FLOAT} [A-Za-z]*)\s*{RE_N}',
                str_operation=str_to_quantity,
                repeats=False,
                convert=False,
            ),
            Quantity(
                'md_thermostat_mass',
                rf'{RE_N} *\| *thermostat effective mass\s*=\s*({RE_FLOAT})\s*{RE_N}',
                str_operation=str_to_quantity,
                repeats=False,
                convert=False,
            ),
            Quantity(
                'md_thermostat_units',
                r'Thermostat\s*units\s*for\s*molecular\s*dynamics\s*:\s*([A-Za-z\^\-0-9]*)',
                str_operation=str_to_ureg,
                repeats=False,
                convert=False,
            ),
            Quantity(
                'md_calculation_info',
                rf'{RE_N} *Advancing structure using '
                rf'Born-Oppenheimer Molecular Dynamics:\s*{RE_N}'
                rf' *Complete information for previous time-step:'
                rf'([\s\S]+?)((?:{RE_N}{RE_N}|\| Nose-Hoover Hamiltonian\s*:'
                rf'\s*[Ee\d\.\-\+]+ eV))',
                str_operation=str_to_md_calculation_info,
                repeats=False,
                convert=False,
            ),
            Quantity(
                'md_system_info',
                rf'Atomic structure.*as used in the preceding time step:\s*{RE_N}'
                rf'([\s\S]+?)((?:{RE_N}{RE_N}|\s*Begin self-consistency loop))',
                repeats=False,
                convert=False,
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'positions',
                            rf'atom +({RE_FLOAT})\s+({RE_FLOAT})\s+({RE_FLOAT})',
                            dtype=np.dtype(np.float64),
                            repeats=True,
                        ),
                        Quantity(
                            'velocities',
                            rf'velocity\s+({RE_FLOAT})\s+({RE_FLOAT})\s+({RE_FLOAT})',
                            dtype=np.dtype(np.float64),
                            repeats=True,
                        ),
                    ]
                ),
            ),
        ]

        tail = '|'.join(
            [
                r'Time for this force evaluation\s*:\s*[s \d\.]+',
                r'Final output of selected total energy values',
                r'No geometry change',
                r'Leaving FHI\-aims',
                r'\Z',
            ]
        )

        self._quantities = [
            Quantity(
                'version',
                r'(?:Version|FHI\-aims version)\s*\:*\s*([\d\.]+)\s*',
                repeats=False,
            ),
            Quantity(
                'program_compilation_date',
                r'Compiled on ([\d\/]+)',
                repeats=False,
            ),
            Quantity(
                'program_compilation_time',
                r'at (\d+\:\d+\:\d+)',
                repeats=False,
            ),
            Quantity('compilation_host', r'on host ([\w\.\-]+)', repeats=False),
            date_time,
            Quantity(
                'cpu1_start',
                r'Time zero on CPU 1\s*:\s*([0-9\-E\.]+)\s*(?P<__unit>\w+)\.',
                repeats=False,
                dtype=float,
            ),
            Quantity(
                'wall_start',
                r'Internal wall clock time zero\s*:'
                r'\s*([0-9\-E\.]+)\s*(?P<__unit>\w+)\.',
                repeats=False,
                dtype=float,
            ),
            Quantity('raw_id', r'aims_uuid\s*:\s*([\w\-]+)', repeats=False),
            Quantity(
                'x_fhi_aims_number_of_tasks',
                r'Using\s*(\d+)\s*parallel tasks',
                repeats=False,
            ),
            Quantity(
                'parallel_task_nr',
                r'Task\s*(\d+)\s*on host',
                repeats=True,
            ),
            Quantity(
                'parallel_task_host',
                r'Task\s*\d+\s*on host\s*([\s\S]+?)reporting',
                repeats=True,
                flatten=False,
            ),
            Quantity(
                'fhi_aims_files',
                r'(?:FHI\-aims file:|Parsing)\s*([\w\/\.]+)',
                repeats=True,
            ),
            Quantity(
                'array_size_parameters',
                r'Basic array size parameters:\s*([\|:\s\w\.\/]+:\s*\d+)',
                repeats=False,
                str_operation=str_to_array_size_parameters,
                convert=False,
            ),
            Quantity(
                'controlInOut_hse_unit',
                r'hse_unit: Unit for the HSE06 hybrid functional '
                r'screening parameter set to\s*(\w)',
                str_operation=str_to_unit,
                repeats=False,
            ),
            Quantity(
                'controlInOut_hybrid_xc_coeff',
                r'hybrid_xc_coeff: Mixing coefficient for hybrid-functional '
                r'exact exchange modified to\s*([\d\.]+)',
                repeats=False,
                dtype=float,
            ),
            Quantity(
                'k_grid', rf'{RE_N} *Found k-point grid:\s*([\d ]+)', repeats=False
            ),  # taken from tests/data/fhi_aims
            Quantity(
                'controlInOut_MD_time_step',
                rf'{RE_N} *Molecular dynamics time step\s*=\s*([\d\.]+)\s*'
                rf'(?P<__unit>[\w]+)',
                repeats=False,
            ),
            Quantity(
                'controlInOut_relativistic',
                rf'{RE_N} *Scalar relativistic treatment of kinetic energy:'
                rf'\s*([\w\- ]+)',
                repeats=False,
            ),
            Quantity(
                'controlInOut_relativistic',
                rf'{RE_N} *(Non-relativistic) treatment of kinetic energy',
                repeats=False,
            ),
            Quantity(
                'controlInOut_relativistic_threshold',
                rf'{RE_N} *Threshold value for ZORA:\s*([\d\.Ee\-\+])',
                repeats=False,
            ),
            Quantity(
                'convergence_energy',
                rf'Convergence accuracy of total energy: *({RE_FLOAT})',
                dtype=float,
            ),
            Quantity(
                'convergence_forces',
                rf'Convergence accuracy for geometry relaxation: '
                rf'Maximum force < *({RE_FLOAT}) *eV/A',
                dtype=float,
                unit='eV/angstrom',
            ),
            Quantity(
                'geometry_relaxation_method',
                r'Geometry relaxation: ([\w ]+)',
                flatten=False,
            ),
            Quantity(
                'controlInOut_xc',
                rf'{RE_N} *XC:\s*(?:Using)*\s*([\w\- ]+) with OMEGA ='
                rf'\s*([\d\.Ee\-\+]+)',
                repeats=False,
                dtype=None,
                flatten=False,
            ),
            Quantity(
                'petukhov',
                rf'{RE_N} *Fixing petukhov mixing factor to\s+(\d?\.[\d]+)',
                repeats=False,
                dtype=np.dtype(np.float64),
            ),
            Quantity(
                'controlInOut_xc',
                r'XC: (?:Running|Using) ([\-\w \(\) ]+)',
                repeats=False,
                flatten=False,
            ),
            Quantity(
                'controlInOut_xc',
                rf'{RE_N} *(Hartree-Fock) calculation starts \.\.\.',
                repeats=False,
                flatten=False,
            ),
            Quantity(
                'band_segment_points',
                r'Plot band\s*\d+\s*\|\s*begin[ \d\.\-]+\s*\|\s*end[ \d\.\-]+\s*\|\s*'
                r'number of points:\s*(\d+)',
                repeats=True,
            ),
            Quantity(
                'species',
                rf'(Reading configuration options for species [\s\S]+?)(?:{RE_N} *'
                rf'Finished|{RE_N} *{RE_N})',
                str_operation=str_to_species_in,
                repeats=False,
            ),
            Quantity(
                'control_inout',
                rf'{RE_N} *Reading file control\.in\.\s*\-*\s*([\s\S]+?)'
                rf'(?:Finished reading input file \'control\.in\'|'
                rf'Input file control\.in ends\.)',
                repeats=False,
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'species',
                            r'Reading configuration options for (species[\s\S]+?)grid'
                            r' points\.',
                            repeats=True,
                            str_operation=str_to_species,
                        ),
                        *molecular_dynamics_quantities,
                    ]
                ),
            ),
            Quantity(
                'control_in_verbatim',
                rf'{RE_N}  *Parsing control\.in([\S\s]*)Completed first pass over'
                rf' input file control\.in',
                repeats=False,
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            'md_controlin',
                            rf' *([\_a-zA-Z\d\-]*MD[\_a-zA-Z\d\-]*)\s+'
                            rf'([a-zA-Z\d\.\-\_\^]+.*){RE_N}',
                            str_operation=str_to_md_control_in,
                            repeats=True,
                            convert=False,
                        )
                    ]
                ),
            ),
            # GW input quantities
            Quantity('gw_flag', RE_GW_FLAG, repeats=False),
            Quantity('anacon_type', rf'{RE_N}\s*anacon_type\s*(\d+)', repeats=False),
            Quantity(
                'gw_analytical_continuation',
                rf'{RE_N} (?:Using)*\s*([\w\-\s]+) for analytical continuation',
                repeats=False,
                flatten=True,
                str_operation=lambda x: [
                    y.lower() for v in x.split(' ') for y in v.split('-')
                ],
            ),
            Quantity('k_grid', rf'{RE_N} *k\_grid\s*([\d ]+)', repeats=False),
            Quantity(
                'freq_grid_type',
                rf'{RE_N}\s*Initialising([\w\-\s]+)time and frequency grids',
                repeats=False,
            ),
            Quantity(
                'n_freq',
                rf'{RE_N}\s*frequency_points\s*(\d+)',
                repeats=False,
                dtype=int,
            ),
            Quantity(
                'frequency_data',
                r'\s*\|*\s*i_freq\s*([\d*\s*.+eE\-\+]+)',
                repeats=True,
                str_operation=str_to_frequency,
            ),
            Quantity(
                'frozen_core',
                rf'{RE_N}\s*frozen_core_scf\s*(\d+)',
                repeats=False,
                dtype=int,
            ),
            Quantity(
                'n_states_gw',
                r'\|\s*Number of Kohn-Sham states \(occupied \+ empty\)\s*\:\s*(\d+)',
                repeats=False,
            ),
            Quantity(
                'gw_self_consistency',
                r'GW Total Energy Calculation([\s\S]+?)\-{5}',
                repeats=True,
                str_operation=str_to_gw_scf,
                convert=False,
            ),
            Quantity(
                'gw_eigenvalues',
                r'(state\s*occ_num\s*e_gs[\s\S]+?)\s*\| Total time',
                str_operation=str_to_gw_eigs,
                repeats=False,
                convert=False,
            ),
            # assign the initial geometry to full scf as no change in structure is done
            # during the initial scf step
            Quantity(
                'lattice_vectors',
                r'Input geometry:\s*\|\s*Unit cell:\s*'
                r'\s*\|\s*([\d\.\-\+eE\s]+)\s*\|\s*([\d\.\-\+eE\s]+)\s*\|\s*([\d\.\-\+eE\s]+)',
                repeats=False,
                unit='angstrom',
                shape=(3, 3),
                dtype=float,
            ),
            Quantity(
                'structure',
                rf'Atomic structure(.|\n)*\| *Atom *x \[A\] *y \[A\] *z \[A\]([\s\S]+?'
                rf'Species[\s\S]+?(?:{RE_N} *{RE_N}| 1\: ))',
                repeats=False,
                convert=False,
                sub_parser=TextParser(quantities=structure_quantities),
            ),
            Quantity(
                'lattice_vectors_reciprocal',
                r'Quantities derived from the lattice vectors:\s*'
                r'\s*\|\s*Reciprocal lattice vector \d:([\d\.\-\+eE\s]+)\s*\|\s*'
                r'Reciprocal lattice vector \d:([\d\.\-\+eE\s]+)\s*\|\s*Reciprocal '
                r'lattice vector \d:([\d\.\-\+eE\s]+)',
                repeats=False,
                unit='1/angstrom',
                shape=(3, 3),
                dtype=float,
            ),
            Quantity(
                'full_scf',
                r'Begin self-consistency loop: Initialization' rf'([\s\S]+?(?:{tail}))',
                repeats=True,
                sub_parser=TextParser(quantities=calculation_quantities),
            ),
            Quantity(
                'geometry_optimization',
                rf'{RE_N} *Geometry optimization: Attempting to predict improved'
                rf' coordinates\.([\s\S]+?(?:{tail}))',
                repeats=True,
                sub_parser=TextParser(quantities=calculation_quantities),
            ),
            Quantity(
                'molecular_dynamics',
                rf'{RE_N} *Molecular dynamics: Attempting to update all '
                rf'nuclear coordinates\.([\s\S]+?(?:{tail}))',
                repeats=True,
                sub_parser=TextParser(
                    quantities=[*calculation_quantities, *molecular_dynamics_quantities]
                ),
            ),
            Quantity(
                'timing',
                r'(Date.+\s+Computational steps[\s\S]+?\Z)',
                sub_parser=TextParser(
                    quantities=[
                        date_time,
                        Quantity(
                            'total_time',
                            r'\| Total time +: +[\d\.]+ s +([\d\.]+) s',
                            dtype=np.float64,
                        ),
                    ]
                ),
            ),
        ]
        # TODO add SOC perturbed eigs, dielectric function
