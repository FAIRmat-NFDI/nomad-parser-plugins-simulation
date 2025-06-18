import re

import numpy as np
from nomad.parsing.file_parser.text_parser import Quantity, TextParser

RE_FLOAT = r'[-+]?\d+\.\d*(?:[Ee][-+]\d+)?'
RE_N = r'[\n\r]'
TOL = 0.01


def str_to_profiling(
    val_in: str,
) -> tuple[list[str], list[str], list[str], list[float], list[float], list[int]]:
    sections = val_in.strip().split('\n\n')
    caller, category, function, cpu_time, wall_time, calls = (
        [],
        [],
        [],
        [],
        [],
        [],
    )
    section_pattern = re.compile(r'(?:Called by ([\S\s]+?):|([\S\s]+?) routines)')
    function_pattern = re.compile(r'(\S+)\s*:')
    cpu_time_pattern = re.compile(
        r'(?:(\d+)h)?(?:(\d+)m)?(?:([\d\.]+)s)? (?:CPU|CPU time)'
    )
    wall_time_pattern = re.compile(
        r'(?:(\d+)h)?(?:(\d+)m)?(?:([\d\.]+)s)? (?:WALL|WALL time|wall|wall time)'
    )
    calls_pattern = re.compile(r'(\d+)\s*calls')
    for section in sections:
        sub_sections = section.split('\n')
        name = section_pattern.findall(sub_sections[0])
        if name:
            caller_name, category_name = name[0][:2]
            sub_sections = sub_sections[1:]
        else:
            caller_name, category_name = '', ''

        for sub_section in sub_sections:
            function_name = function_pattern.search(sub_section)
            if function_name is None:
                continue
            caller.append(caller_name.strip())
            category.append(category_name.strip())
            function.append(function_name.group(1).strip())
            res = cpu_time_pattern.findall(sub_section)
            if len(res) == 0:
                return caller, category, function, cpu_time, wall_time, calls
            time = sum(
                [
                    float(res[0][i]) * 60 ** (2 - i) if res[0][i] else 0
                    for i in range(len(res[0]))
                ]
            )
            cpu_time.append(0.0 if res is None else time)
            res = wall_time_pattern.findall(sub_section)
            wall_time.append(
                0.0
                if not res
                else sum(
                    [
                        float(res[0][i]) * 60 ** (2 - i) if res[0][i] else 0
                        for i in range(len(res[0]))
                    ]
                )
            )
            res = calls_pattern.search(sub_section)
            calls.append(0 if res is None else int(res.group(1)))
    return caller, category, function, cpu_time, wall_time, calls


def str_to_dispersion(val_in: str) -> dict[dict[str, float]]:
    val = [v.split() for v in val_in]
    return {v[0].strip(): {'vdw_radius': float(v[1]), 'C_6': float(v[2])} for v in val}


def str_to_sticks(val_in: str) -> dict[str, list[int]]:
    val = [v.split() for v in val_in.strip().split('\n')]
    return {v[0]: [int(i) for i in v[1:]] for v in val}


def str_to_atom_data(val_in: str) -> list[list[float | str]]:
    val = [
        v.replace('(', ' ').replace(')', ' ').split()
        for v in val_in.strip().split('\n')
    ]
    return [[float(vi) if not vi[0].isalpha() else vi for vi in v] for v in val]


def str_to_arrays(val_in: str) -> tuple[list[str], list[float], list[str]]:
    val = [v for v in val_in.strip().split('\n')]
    pattern = re.compile(r'([\S\s]+?)(\d+\.\d+)\s*Mb\s*\(\s*([\d, ]+)\)')
    names = []
    sizes = []
    dimensions = []
    for v in val:
        res = pattern.search(v)
        if res is None:
            continue
        g = res.groups()
        names.append(g[0].strip())
        sizes.append(float(g[1]))
        dimensions.append(g[2].strip())
    return names, sizes, dimensions


header_quantities = [
    Quantity(
        'program_name_version',
        r'Program\s*(\w+\s*v\.\S+\s*(?:\(svn rev\.\s*\d+\))*)',
    ),
    Quantity(
        'start_date_time',
        r'starts (?:on|\.\.\.\s*Today is)\s*(\w+)\s*at\s*([\d: ]+)',
        flatten=False,
    ),
    Quantity(
        'compile_parallel_version',
        r'(Serial multi\-threaded|Serial|Parallel)\s*version\s*(\(MPI\))*',
        flatten=False,
    ),
    Quantity('nthreads', r', running on\s*(\d+)\s*processor cores', dtype=int),
    Quantity(
        'nproc',
        r'(?:Number of processors in use:|, running on)' r'\s*(\d+)\s*(?:processors)*',
        dtype=int,
    ),
    Quantity('npool', r'npool[\w/]*\s*=\s*(\d+)', dtype=int),
    Quantity('input_filename', r'Reading input from (.+)', flatten=False),
    Quantity('save_directory', r'Reading xml data from directory:\s+?(\S+)'),
    Quantity('ntypx', r'ntypx.*?=\s*(\d+)', dtype=int),
    Quantity('npk', r'npk.*?=\s*(\d+)', dtype=int),
    Quantity('lmaxx', r'lmaxx.*?=\s*(\d+)', dtype=int),
    Quantity('nchix', r'nchix.*?=\s*(\d+)', dtype=int),
    Quantity('ndmx', r'ndmx.*?=\s*(\d+)', dtype=int),
    Quantity('nbrx', r'nbrx.*?=\s*(\d+)', dtype=int),
    Quantity(
        'pseudopotential_report',
        r'(pseudopotential report for atomic species:[\s\S]+?)={10}',
        sub_parser=TextParser(
            quantities=[
                Quantity('species', r'atomic species:\s*(\d+)', dtype=int),
                Quantity(
                    'version',
                    r'pseudo potential version\s*([\d ]+)',
                    flatten=False,
                ),
                Quantity('contents', r'\-\s*([\s\S]+)', flatten=False),
            ]
        ),
        repeats=True,
    ),
    Quantity(
        'gamma_algorithms',
        r'(gamma\-point specific algorithms are used)',
        str_operation=lambda x: True,
    ),
    Quantity(
        'diagonalization_algorithm',
        r'of the eigenvalue problem:\s*a (serial) algorithm will be used',
    ),
    Quantity(
        'atom_radii',
        r'new r_m :\s*[\d\.]+\s*\(alat units\)\s*([\d\.]+) \(a\.u\.\) '
        r'for type\s*(\d+)',
        repeats=True,
    ),
    Quantity(
        'input_positions_cell_dirname',
        r'Atomic positions and unit cell read from directory:\s*(\S+)',
    ),
    Quantity(
        'supercell',
        rf'(?:operation: I \+ \(|translation:)\s*'
        rf'({RE_FLOAT})\s*({RE_FLOAT})\s*({RE_FLOAT})',
        dtype=float,
        repeats=True,
    ),
    Quantity(
        'renormalized_wavefunction',
        r'file (\S+): wavefunction\(s\)([\w ]+)renormalized',
        repeats=True,
    ),
    Quantity('exchange_correlation', r'Exchange\-correlation *= *(.+)', flatten=False),
]

energy_quantities = [
    Quantity(
        'energy_total', rf'total energy\s*=\s*({RE_FLOAT})', dtype=float, unit='rydberg'
    ),
    Quantity(
        'energy_total_harris_foulkes_estimate',
        rf'Harris-Foulkes estimate\s*=\s*({RE_FLOAT})',
        dtype=float,
        unit='rydberg',
    ),
    Quantity(
        'energy_total_accuracy_estimate',
        rf'estimated scf accuracy\s*<\s*({RE_FLOAT})',
        dtype=float,
        unit='rydberg',
    ),
    Quantity(
        'energy_total_paw_all_electron',
        rf'total all-electron energy\s*=\s*({RE_FLOAT})',
        dtype=float,
        unit='rydberg',
    ),
]

calculation_quantities = [
    Quantity(
        'charge_negative_spin',
        rf'negative rho \(up, down\):\s*({RE_FLOAT})\s*({RE_FLOAT})',
        dtype=float,
    ),
    Quantity(
        'magnetic_moments',
        rf'atom:\s*(\d+)\s*charge:\s*({RE_FLOAT})\s*magn:'
        rf'\s*({RE_FLOAT})\s*constr:\s*({RE_FLOAT})',
        repeats=True,
        dtype=float,
    ),
    Quantity(
        'energies',
        r'(total energy[\s\S]+?)\n\s*\n',
        sub_parser=TextParser(quantities=energy_quantities),
    ),
    Quantity(
        'magnetization_total',
        rf'total magnetization\s*=\s*({RE_FLOAT})\s*Bohr mag/cell',
        dtype=float,
        unit='bohr_magneton',
    ),
    Quantity(
        'magnetization_absolute',
        rf'absolute magnetization\s*=\s*({RE_FLOAT})' rf'\s*Bohr mag/cell',
        dtype=float,
        unit='bohr_magneton',
    ),
    Quantity(
        'fermi_energy_shift',
        rf'Fermi energy shift \(Ry\) *= *({RE_FLOAT}) +({RE_FLOAT})',
        dtype=np.float64,
    ),
]

general_quantities = [
    Quantity(
        'dispersion',
        r'atom\s*VdW radius\s*C_6\s*([\s\S]+?)\n\s*\n',
        str_operation=str_to_dispersion,
        convert=False,
    ),
    Quantity(
        'g_vector_sticks',
        r'G-vecs:\s*dense\s*smooth\s*PW\s*([\s\S]+?)\n *\n',
        str_operation=str_to_sticks,
        convert=False,
    ),
    Quantity('ibrav', r'bravais\-lattice index\s*=\s*(\d+)', dtype=int),
    Quantity(
        'alat',
        rf'lattice parameter \((?:alat|a_0)\)\s*=\s*({RE_FLOAT})',
        unit='bohr',
        dtype=float,
    ),
    Quantity(
        'cell_volume',
        rf'unit-cell volume\s*=\s*({RE_FLOAT})',
        unit='bohr**3',
        dtype=float,
    ),
    Quantity('number_of_atoms', r'number of atoms/cell\s*=\s*(\d+)', dtype=int),
    Quantity(
        'number_of_species',
        r'number of atomic types\s*=\s*(\d+)',
        dtype=int,
    ),
    Quantity(
        'number_of_electrons',
        r'(number of electrons\s*=[^\n]*)',
        sub_parser=TextParser(
            quantities=[
                Quantity(
                    'total',
                    rf'number of electrons\s*=\s*({RE_FLOAT})',
                    dtype=float,
                ),
                Quantity('up', rf'up:\s*({RE_FLOAT})', dtype=float),
                Quantity('down', rf'down:\s*({RE_FLOAT})', dtype=float),
            ],
        ),
    ),
    Quantity(
        'number_of_states',
        r'number of Kohn\-Sham states\s*=\s*(\d+)',
        dtype=int,
    ),
    Quantity(
        'wavefunction_cutoff',
        rf'kinetic\-energy cutoff\s*=\s*({RE_FLOAT})',
        dtype=float,
        unit='rydberg',
    ),
    Quantity(
        'density_cutoff',
        rf'charge density cutoff\s*=\s*({RE_FLOAT})',
        dtype=float,
        unit='rydberg',
    ),
    Quantity(
        'fock_cutoff',
        rf'cutoff for Fock operator\s*=\s*({RE_FLOAT})',
        dtype=float,
        unit='rydberg',
    ),
    Quantity(
        'scf_threshold_energy_change',
        rf'convergence threshold\s*=\s*({RE_FLOAT})',
        dtype=float,
        unit='rydberg',
    ),
    Quantity(
        'potential_mixing_beta',
        rf'mixing beta\s*=\s*({RE_FLOAT})',
        dtype=float,
    ),
    Quantity('mixing_scheme', r'number of iterations used\s*=\s*(\d+)\s*(.*)mixing'),
    Quantity(
        'xc_functional',
        r'Exchange\-correlation\s*=\s*(.+)\s*(\([\d ]+\))',
        convert=False,
        flatten=False,
    ),
    Quantity(
        'exact_exchange_fraction',
        rf'EXX\-fraction\s*=\s*({RE_FLOAT})',
        dtype=float,
    ),
    Quantity('md_max_steps', r'\s*nstep\s*=\s*(\d+)', dtype=int),
    Quantity(
        'spin_orbit_mode',
        r'\s*(.*?)\s*calculation\s*(with(?:out)?)\s*spin-orbit',
    ),
    Quantity(
        'berry_efield',
        r'(Using Berry phase electric field[\s\S]+?'
        r'Number of iterative cycles\s*:\s*\d+)',
        sub_parser=TextParser(
            quantities=[
                Quantity('direction', r'Direction\s*:\s*(\d+)', dtype=int),
                Quantity(
                    'intensity',
                    rf'Intensity \((?:Ry\s*)?a.u.\)\s*:\s*({RE_FLOAT})',
                    dtype=float,
                ),
                Quantity('strings', r'Strings composed by\s*:\s*(\d+)', dtype=int),
                Quantity(
                    'niter',
                    r'Number of iterative cycles\s*:\s*(\d+)',
                    dtype=int,
                ),
            ]
        ),
    ),
    Quantity('assume_isolated', r'Assuming isolated system,\s*(.*?)\s*method'),
    Quantity(
        'celldm',
        rf'celldm\(1\)=\s*({RE_FLOAT})\s*celldm\(2\)=\s*({RE_FLOAT})\s*celldm\(3\)=\s*({RE_FLOAT})\s*'
        rf'celldm\(4\)=\s*({RE_FLOAT})\s*celldm\(5\)=\s*({RE_FLOAT})\s*celldm\(6\)=\s*({RE_FLOAT})\s*',
        dtype=float,
        unit='bohr',
    ),
    Quantity('units', r'crystal axes: \(cart\. coord\. in units of ([\w ]+)\)\s*'),
    Quantity(
        'simulation_cell',
        r'a\(1\) = \(([\-\d\. ]+)\)\s*a\(2\) = \(([\-\d\. ]+)\)\s*a\(3\) = '
        r'\(([\-\d\. ]+)\)\s*',
        dtype=float,
        shape=(3, 3),
    ),
    Quantity(
        'reciprocal_cell_units',
        r'reciprocal axes: \(cart\. coord\. in units of ([\w ]+)\)\s*',
    ),
    Quantity(
        'reciprocal_cell',
        r'b\(1\) = \(([\-\d\. ]+)\)\s*b\(2\) = \(([\-\d\. ]+)\)\s*b\(3\) = '
        r'\(([\-\d\. ]+)\)\s*',
        dtype=float,
        shape=(3, 3),
    ),
    Quantity(
        'pseudopotential',
        r'(PseudoPot\. #[\s\S]+?\n\s*\n)',
        repeats=True,
        sub_parser=TextParser(
            quantities=[
                Quantity('idx', r'PseudoPot\. # (\d+)'),
                Quantity('label', r'for (\w+)'),
                Quantity('filename', r'read from file:?\s*\s*(\S+)'),
                Quantity('md5sum', r'MD5 check sum:\s*(\S+)'),
                Quantity('type', r'Pseudo is\s*(.*?),', flatten=False),
                Quantity('valence', rf'Zval\s*=\s*({RE_FLOAT})'),
                Quantity('comment', r'\s*(.+?)\s*Using radial', flatten=False),
                Quantity(
                    'n_radial_grid_points',
                    r'Using radial grid of *(\d+) points',
                    dtype=int,
                ),
                Quantity(
                    'integral_ndirections',
                    r'Setup to integrate on\s*(\d+)\s+directions:',
                ),
                Quantity('integral_lmax_exact', r'integral exact up to l =\s*(\d+)'),
                Quantity(
                    'augmentation_shape',
                    r'Shape of augmentation charge:\s*(.*?)\s*',
                ),
                Quantity('ndmx', r'grid of\s*(\d+) points', dtype=int),
                Quantity('nbeta', r',\s*(\d+) beta functions', dtype=int),
                Quantity('beta', r'l\((\d+)\)\s*=\s*(\d+)', repeats=True, dtype=int),
                Quantity(
                    'ncoefficients',
                    r'Q\(r\) pseudized with\s*(\d+)\s*coefficients',
                    dtype=int,
                ),
                Quantity(
                    'rinner',
                    r'rinner\s*=\s*([\-\d\.\s]+)',
                    str_operation=lambda x: ' '.join(x.split()),
                    convert=False,
                ),
            ]
        ),
    ),
    Quantity('point_group', r', *(.+?) +point group', dtype=str, flatten=False),
    Quantity(
        'atom_species_pp',
        r'atomic species\s*valence\s*mass\s*pseudopotential([\s\S]+?)\n\s*\n',
        str_operation=str_to_atom_data,
        convert=False,
    ),
    Quantity(
        'starting_magnetization',
        r'Starting magnetic structure\s*atomic species\s*magnetization'
        r'([\s\S]+?)\n\s*\n',
        str_operation=str_to_atom_data,
        convert=False,
    ),
    Quantity(
        'md_cell_mass',
        rf'cell mass\s*=\s*({RE_FLOAT})\s*AMU/\(a\.u\.\)\^2',
        dtype=float,
    ),
    Quantity(
        'symmetry',
        r'(\d+\s*Sym\.\s*Ops\.[\s\S]+?)\n\s*\n',
        sub_parser=TextParser(
            quantities=[
                Quantity('nsymm', r'(\d+) Sym\.', dtype=int),
                Quantity(
                    'symm_inversion',
                    r'\((\S+)\s*inversion',
                    str_operation=lambda x: 'with' in x,
                    convert=False,
                ),
                Quantity(
                    'nsymm_with_fractional_translation',
                    r'(\d+)\s*have fractional translation',
                    dtype=int,
                ),
                Quantity(
                    'nsymm_ignored',
                    r'(\d+)\s*additional sym\.ops\. were found but ignored',
                    dtype=int,
                ),
            ]
        ),
    ),
    Quantity(
        'labels_positions',
        r'(Cartesian axes\s*site n\.\s*atom.+?positions[\s\S]+?)\n\s*\n',
        repeatas=False,
        sub_parser=TextParser(
            quantities=[
                Quantity('axes', r'(\w+)\s*axes'),
                Quantity('units', r'site n\.\s*atom.+?positions\s*\(\s*(\S+)'),
                Quantity('labels', r' ([A-Z][a-z]?)\S* ', repeats=True),
                Quantity(
                    'positions',
                    rf'=\s*\(\s*({RE_FLOAT}\s*{RE_FLOAT}\s*{RE_FLOAT})\s*\)',
                    repeats=True,
                    dtype=float,
                ),
            ]
        ),
    ),
    Quantity(
        'k_points',
        r'(number of k points=[\s\S]+?)\n\s*\n',
        repeats=False,
        sub_parser=TextParser(
            quantities=[
                Quantity('nk', r'number of k points=\s*(\d+)', dtype=int),
                Quantity(
                    'gaussian_broadening',
                    rf'gaussian broad\. \(Ry\)= +({RE_FLOAT})',
                    dtype=float,
                    unit='rydberg',
                ),
                Quantity('n_gauss', r'ngauss *= *(\d+)', dtype=np.int32),
                Quantity(
                    'smearing',
                    r'([\w\-]+)\s*(?:broad|smearing|method),?',
                    dtype=str,
                    flatten=False,
                ),
                Quantity(
                    'width',
                    rf'width\s*\(Ry\)=\s*({RE_FLOAT})',
                    dtype=float,
                    unit='rydberg',
                ),
                Quantity('units', r'cart\. coord\. in units (2pi/alat)'),
                Quantity('ik', r'k\(\s*(\d+)\s*\)', repeats=True),
                Quantity(
                    'points',
                    rf'=\s*\(\s*({RE_FLOAT}\s*{RE_FLOAT}\s*{RE_FLOAT})\)',
                    repeats=True,
                    dtype=float,
                ),
                Quantity('wk', rf',\s*wk\s*=\s*({RE_FLOAT})', repeats=True),
                Quantity('warning', r'(Number of k-points >= 100: set verbosity)'),
            ]
        ),
    ),
    Quantity(
        'dense_grid',
        rf'(?:G\s+cutoff\s*=\s*({RE_FLOAT})\s*\(|Dense\s*grid:)\s*(\d+)\s*'
        r'G\-vectors\)*\s*FFT\s+(?:dimensions|grid):\s*\(\s*([\d ,]+)\)',
        str_operation=lambda x: x.replace(',', ' ').split(),
    ),
    Quantity(
        'smooth_grid',
        rf'(?:G\s+cutoff\s*=\s*({RE_FLOAT})\s*\(|Smooth\s*grid:)\s*(\d+)\s*'
        r'G\-vectors\s*(?:smooth grid|FFT dimensions)\s*:\s*\(\s*([\d ,]+)\)',
        str_operation=lambda x: x.replace(',', ' ').split(),
    ),
    Quantity(
        'alpha_ewald',
        rf' Alpha used in Ewald sum = *({RE_FLOAT})',
        dtype=np.float64,
    ),
    Quantity(
        'core_charge_realspace',
        r'(Real space treatment of Q\(r\))',
        str_operation=lambda x: True,
    ),
    Quantity(
        'input_occupation',
        r'Occupations\s*read\s*from\s*input\s*'
        r'(?:Spin-up)?([\d\.Ee\s]+)(?:Spin-down)?([\d\.Ee\s]+)',
        str_operation=lambda x: np.array(x.strip().split(), dtype=float),
        convert=False,
    ),
    Quantity(
        'allocated_arrays',
        r'allocated arrays\s*est\. size \(Mb\)\s*dimensions([\s\S]+?)Largest',
        str_operation=str_to_arrays,
        convert=False,
    ),
    Quantity(
        'temporary_arrays',
        r'temporary arrays\s*est\. size \(Mb\)\s*dimensions([\s\S]+?)\n\s*\n',
        str_operation=str_to_arrays,
        convert=False,
    ),
    Quantity(
        'martyna_tuckerman_parameters',
        rf'alpha, beta MT =\s*({RE_FLOAT})\s*({RE_FLOAT})',
        dtype=float,
    ),
    Quantity(
        'core_charge_check',
        rf'Check: negative/imaginary core charge\s*=\s*'
        rf'({RE_FLOAT})\s*({RE_FLOAT})',
        dtype=float,
    ),
    Quantity(
        'input_potential_recalculated_file',
        r'The potential is recalculated from file\s*:\s*(\S+)',
    ),
    Quantity(
        'starting_density_file',
        r'The initial density is read from file\s*:\s*(\S+)',
    ),
    Quantity(
        'starting_potential',
        r'Initial potential from\s*(.+)',
        flatten=False,
    ),
    Quantity(
        'starting_charge_negative',
        rf'Check: negative starting charge\s*=\s*({RE_FLOAT})',
        dtype=float,
    ),
    Quantity(
        'initial_charge',
        rf'starting charge\s*({RE_FLOAT})\s*, renormalised to\s*({RE_FLOAT})',
        dtype=float,
    ),
    Quantity('starting_wfc', r'Starting wfc\s*(.+)', flatten=False),
    Quantity(
        'time_setup_cpu1_end',
        rf'total cpu time spent up to now is\s*({RE_FLOAT})',
    ),
    Quantity(
        'per_process_mem',
        r'per\-process dynamical memory:\s*([\d\.]+)\s*Mb',
        dtype=float,
        unit='mebibyte',
    ),
    Quantity(
        'profiling',
        r'(.+?CPU.+?WALL[\s\S]+?)(?:This run|\Z)',
        repeats=False,
        str_operation=str_to_profiling,
        convert=False,
    ),
    Quantity(
        'memory',
        r'per\-process dynamical memory:\s*([\d\.]+)\s*Mb',
        dtype=float,
        unit='mebibyte',
    ),
    Quantity('output_datafile', r'Writing output data file ([\w\.]+)', dtype=str),
]

scf_iteration_quantities = [
    Quantity('number', r'n\s*#\s*(\d+)'),
    Quantity('ecutwfc', r'ecut=\s*([\d\.]+)', unit='rydberg'),
    Quantity('beta', r'beta=\s*([\d\.]+)'),
    Quantity(
        'total_time',
        r'total cpu time spent up to now is\s*([\d\.]+)',
        unit='s',
    ),
    Quantity(
        'time',
        rf'av\.it\.: *({RE_FLOAT})',
        dtype=float,
    ),
    Quantity(
        'threshold',
        rf'thresh *= *({RE_FLOAT})',
        dtype=float,
    ),
    Quantity(
        'alpha_mix',
        rf'alpha\_mix *= *({RE_FLOAT})',
        dtype=float,
    ),
    Quantity(
        'ddv_scf',
        rf'\|ddv_scf\|\^2 *= *({RE_FLOAT})',
        dtype=float,
    ),
]

tail_quantities = [
    Quantity(
        'end_date_time',
        r'was terminated on:\s*([\d: ]+)\s*(\w+)',
        flatten=False,
    ),
    Quantity('job_done', r'(JOB DONE)'),
]


# origin: espresso-5.4.0/Modules/funct.f90
# update:
# . New exchange-correlation functionals exist in
# .     espresso-6.5.0/Modules/funct.f90
#   short comments mark the corresponding new metainfo
_exchange_map = [
    None,
    {
        'xc_terms': [
            {
                'XC_functional_name': 'LDA_X',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_iexch_name': 'sla',
            'x_qe_xc_iexch_comment': 'Slater (alpha=2/3)',
            'x_qe_xc_iexch': 1,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'LDA_X',
                'XC_functional_parameters': {'alpha': 1.0},
            }
        ],
        'xc_section_method': {
            'x_qe_xc_iexch_name': 'sl1',
            'x_qe_xc_iexch_comment': 'Slater (alpha=1.0)',
            'x_qe_xc_iexch': 2,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'x_qe_LDA_X_RELATIVISTIC',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_iexch_name': 'rxc',
            'x_qe_xc_iexch_comment': 'Relativistic Slater',
            'x_qe_xc_iexch': 3,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'OEP_EXX',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_iexch_name': 'oep',
            'x_qe_xc_iexch_comment': 'Optimized Effective Potential',
            'x_qe_xc_iexch': 4,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'HF_X',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_iexch_name': 'hf',
            'x_qe_xc_iexch_comment': 'Hartree-Fock',
            'x_qe_xc_iexch': 5,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'HF_X',
                'exx_compute_weight': lambda exx: exx,
                'XC_functional_weight': 0.25,
            },
            {
                'XC_functional_name': 'LDA_X',
                'exx_compute_weight': lambda exx: (1.0 - exx),
                'XC_functional_weight': 0.75,
            },
        ],
        'xc_section_method': {
            'x_qe_xc_iexch_name': 'pb0x',
            'x_qe_xc_iexch_comment': 'PBE0 (Slater*0.75+HF*0.25)',
            'x_qe_xc_iexch': 6,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'HF_X',
                'exx_compute_weight': lambda exx: exx,
                'XC_functional_weight': 0.20,
            },
            {
                'XC_functional_name': 'LDA_X',
                'exx_compute_weight': lambda exx: (1.0 - exx),
                'XC_functional_weight': 0.8,
            },
        ],
        'xc_section_method': {
            'x_qe_xc_iexch_name': 'b3lp',
            'x_qe_xc_iexch_comment': 'B3LYP(Slater*0.80+HF*0.20)',
            'x_qe_xc_iexch': 7,
        },
    },
    # LDA_X_KZK is not part of libXC. Look up it at
    # 'https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-info/wikis/metainfo/XC-functional'
    {
        'xc_terms': [
            {
                'XC_functional_name': 'LDA_X_KZK',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_iexch_name': 'kzk',
            'x_qe_xc_iexch_comment': 'Finite-size corrections',
            'x_qe_xc_iexch': 8,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'HF_X',
                'exx_compute_weight': lambda exx: exx,
                'XC_functional_weight': 0.218,
            },
            {
                'XC_functional_name': 'LDA_X',
                'exx_compute_weight': lambda exx: (1.0 - exx),
                'XC_functional_weight': 0.782,
            },
        ],
        'xc_section_method': {
            'x_qe_xc_iexch_name': 'x3lp',
            'x_qe_xc_iexch_comment': 'X3LYP(Slater*0.782+HF*0.218)',
            'x_qe_xc_iexch': 9,
        },
    },
    # update for espresso-6.5.0: KLI
    {
        'xc_terms': [
            {
                'XC_functional_name': 'LDA_X_KLI',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_iexch_name': 'kli',
            'x_qe_xc_iexch_comment': 'KLI aproximation for exx',
            'x_qe_xc_iexch': 10,
        },
    },
]

# Correlation functionals UNchanged between espresso v5.4 & v6.5
_correlation_map = [
    None,
    {
        'xc_terms': [
            {
                'XC_functional_name': 'LDA_C_PZ',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_icorr_name': 'pz',
            'x_qe_xc_icorr_comment': 'Perdew-Zunger',
            'x_qe_xc_icorr': 1,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'LDA_C_VWN',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_icorr_name': 'vwn',
            'x_qe_xc_icorr_comment': 'Vosko-Wilk-Nusair',
            'x_qe_xc_icorr': 2,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'LDA_C_LYP',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_icorr_name': 'lyp',
            'x_qe_xc_icorr_comment': 'Lee-Yang-Parr',
            'x_qe_xc_icorr': 3,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'LDA_C_PW',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_icorr_name': 'pw',
            'x_qe_xc_icorr_comment': 'Perdew-Wang',
            'x_qe_xc_icorr': 4,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'LDA_C_WIGNER',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_icorr_name': 'wig',
            'x_qe_xc_icorr_comment': 'Wigner',
            'x_qe_xc_icorr': 5,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'LDA_C_HL',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_icorr_name': 'hl',
            'x_qe_xc_icorr_comment': 'Hedin-Lunqvist',
            'x_qe_xc_icorr': 6,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'LDA_C_OB_PZ',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_icorr_name': 'obz',
            'x_qe_xc_icorr_comment': 'Ortiz-Ballone form for PZ',
            'x_qe_xc_icorr': 7,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'LDA_C_OB_PW',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_icorr_name': 'obw',
            'x_qe_xc_icorr_comment': 'Ortiz-Ballone form for PW',
            'x_qe_xc_icorr': 8,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'LDA_C_GL',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_icorr_name': 'gl',
            'x_qe_xc_icorr_comment': 'Gunnarson-Lunqvist',
            'x_qe_xc_icorr': 9,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'LDA_C_KZK',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_icorr_name': 'kzk',
            'x_qe_xc_icorr_comment': 'Finite-size corrections',
            'x_qe_xc_icorr': 10,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'LDA_C_VWN_RPA',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_icorr_name': 'vwn-rpa',
            'x_qe_xc_icorr_comment': 'Vosko-Wilk-Nusair, alt param',
            'x_qe_xc_icorr': 11,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'LDA_C_VWN',
                'XC_functional_weight': 0.19,
            },
            {
                'XC_functional_name': 'LDA_C_LYP',
                'XC_functional_weight': 0.81,
            },
        ],
        'xc_section_method': {
            'x_qe_xc_icorr_name': 'b3lp',
            'x_qe_xc_icorr_comment': 'B3LYP (0.19*vwn+0.81*lyp)',
            'x_qe_xc_icorr': 12,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'LDA_C_VWN_RPA',
                'XC_functional_weight': 0.19,
            },
            {
                'XC_functional_name': 'LDA_C_LYP',
                'XC_functional_weight': 0.81,
            },
        ],
        'xc_section_method': {
            'x_qe_xc_icorr_name': 'b3lpv1r',
            'x_qe_xc_icorr_comment': 'B3LYP-VWN-1-RPA (0.19*vwn_rpa+0.81*lyp)',
            'x_qe_xc_icorr': 13,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'LDA_C_VWN_RPA',
                'XC_functional_weight': 0.129,
            },
            {
                'XC_functional_name': 'LDA_C_LYP',
                'XC_functional_weight': 0.871,
            },
        ],
        'xc_section_method': {
            'x_qe_xc_icorr_name': 'x3lp',
            'x_qe_xc_icorr_comment': 'X3LYP (0.129*vwn_rpa+0.871*lyp)',
            'x_qe_xc_icorr': 14,
        },
    },
]

# New 'exchange_gradient_correction' functionals for q-espresso (qe) v6.5
#    igcx=[1..28] unchanged between qe-v5.4 & v6.5
# New additions: igcx=[29..42]

_exchange_gradient_correction_map = [
    None,
    {
        'xc_terms': [
            {
                'XC_functional_name': 'GGA_X_B88',
            }
        ],
        'xc_terms_remove': [
            {
                'XC_functional_name': 'LDA_X',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_igcx_name': 'b88',
            'x_qe_xc_igcx_comment': 'Becke88 (beta=0.0042)',
            'x_qe_xc_igcx': 1,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'GGA_X_PW91',
            }
        ],
        'xc_terms_remove': [
            {
                'XC_functional_name': 'LDA_X',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_igcx_name': 'ggx',
            'x_qe_xc_igcx_comment': 'Perdew-Wang 91',
            'x_qe_xc_igcx': 2,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'GGA_X_PBE',
            }
        ],
        'xc_terms_remove': [
            {
                'XC_functional_name': 'LDA_X',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_igcx_name': 'pbx',
            'x_qe_xc_igcx_comment': 'Perdew-Burke-Ernzenhof exch',
            'x_qe_xc_igcx': 3,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'GGA_X_PBE_R',
            }
        ],
        'xc_terms_remove': [
            {
                'XC_functional_name': 'LDA_X',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_igcx_name': 'rpb',
            'x_qe_xc_igcx_comment': 'revised PBE by Zhang-Yang',
            'x_qe_xc_igcx': 4,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'GGA_XC_HCTH_120',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_igcx_name': 'hcth',
            'x_qe_xc_igcx_comment': 'Cambridge exch, Handy et al, HCTH/120',
            'x_qe_xc_igcx': 5,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'GGA_X_OPTX',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_igcx_name': 'optx',
            'x_qe_xc_igcx_comment': "Handy's exchange functional",
            'x_qe_xc_igcx': 6,
        },
    },
    {
        # igcx=7 is not defined in 5.4's funct.f90
        #        definition taken from 5.0, which did not have separate imeta
        'xc_terms': [
            {
                'XC_functional_name': 'MGGA_X_TPSS',
            }
        ],
        'xc_terms_remove': [
            {
                'XC_functional_name': 'LDA_X',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_igcx_name': 'tpss',
            'x_qe_xc_igcx_comment': 'TPSS Meta-GGA (Espresso-version < 5.1)',
            'x_qe_xc_igcx': 7,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'GGA_X_PBE',
                'XC_functional_weight': 0.75,
                'exx_compute_weight': lambda exx: (1.0 - exx),
            }
        ],
        'xc_terms_remove': [
            {
                'XC_functional_name': 'LDA_X',
                'XC_functional_weight': 0.75,
                'exx_compute_weight': lambda exx: (1.0 - exx),
            }
        ],
        'xc_section_method': {
            'x_qe_xc_igcx_name': 'pb0x',
            'x_qe_xc_igcx_comment': 'PBE0 (PBE exchange*0.75)',
            'x_qe_xc_igcx': 8,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'GGA_X_B88',
                'XC_functional_weight': 0.72,
                'exx_compute_weight': lambda exx: 0.72 if abs(exx) > TOL else 1.0,
            }
        ],
        'xc_terms_remove': [
            {
                'XC_functional_name': 'LDA_X',
                'XC_functional_weight': 0.8,
                'exx_compute_weight': lambda exx: (1.0 - exx),
            }
        ],
        'xc_section_method': {
            'x_qe_xc_igcx_name': 'b3lp',
            'x_qe_xc_igcx_comment': 'B3LYP (Becke88*0.72)',
            'x_qe_xc_igcx': 9,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'GGA_X_PBE_SOL',
            }
        ],
        'xc_terms_remove': [
            {
                'XC_functional_name': 'LDA_X',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_igcx_name': 'psx',
            'x_qe_xc_igcx_comment': 'PBEsol exchange',
            'x_qe_xc_igcx': 10,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'GGA_X_WC',
            }
        ],
        'xc_terms_remove': [
            {
                'XC_functional_name': 'LDA_X',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_igcx_name': 'wcx',
            'x_qe_xc_igcx_comment': 'Wu-Cohen',
            'x_qe_xc_igcx': 11,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'HYB_GGA_XC_HSE06',
                'exx_compute_weight': lambda exx: 1.0 if (abs(exx) > TOL) else 0.0,
            },
            {
                'XC_functional_name': 'GGA_X_PBE',
                'exx_compute_weight': lambda exx: 0.0 if (abs(exx) > TOL) else 1.0,
            },
        ],
        'xc_terms_remove': [
            {
                'XC_functional_name': 'LDA_X',
            },
            {
                'XC_functional_name': 'GGA_C_PBE',
                'exx_compute_weight': lambda exx: 1.0 if (abs(exx) > TOL) else 0.0,
            },
        ],
        'xc_section_method': {
            'x_qe_xc_igcx_name': 'hse',
            'x_qe_xc_igcx_comment': 'HSE screened exchange',
            'x_qe_xc_igcx': 12,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'GGA_X_RPW86',
            }
        ],
        'xc_terms_remove': [
            {
                'XC_functional_name': 'LDA_X',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_igcx_name': 'rw86',
            'x_qe_xc_igcx_comment': 'revised PW86',
            'x_qe_xc_igcx': 13,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'GGA_X_PBE',
            }
        ],
        'xc_terms_remove': [
            {
                'XC_functional_name': 'LDA_X',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_igcx_name': 'pbe',
            'x_qe_xc_igcx_comment': 'same as PBX, back-comp.',
            'x_qe_xc_igcx': 14,
        },
    },
    {
        # igcx=15 is not defined in 5.4's funct.f90
        #        definition taken from 5.0, which did not have separate imeta
        'xc_terms': [
            {
                'XC_functional_name': 'MGGA_X_TB09',
            }
        ],
        'xc_terms_remove': [
            {
                'XC_functional_name': 'LDA_X',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_igcx_name': 'tb09',
            'x_qe_xc_igcx_comment': 'TB09 Meta-GGA (Espresso-version < 5.1)',
            'x_qe_xc_igcx': 15,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'GGA_X_C09X',
            }
        ],
        'xc_terms_remove': [
            {
                'XC_functional_name': 'LDA_X',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_igcx_name': 'c09x',
            'x_qe_xc_igcx_comment': 'Cooper 09',
            'x_qe_xc_igcx': 16,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'GGA_X_SOGGA',
            }
        ],
        'xc_terms_remove': [
            {
                'XC_functional_name': 'LDA_X',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_igcx_name': 'sox',
            'x_qe_xc_igcx_comment': 'sogga',
            'x_qe_xc_igcx': 17,
        },
    },
    {
        # igcx=18 is not defined in 5.4's funct.f90
        #        definition taken from 5.0, which did not have separate imeta
        'xc_terms': [
            {
                'XC_functional_name': 'MGGA_X_M06_L',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_igcx_name': 'm6lx',
            'x_qe_xc_igcx_comment': 'M06L Meta-GGA (Espresso-version < 5.1)',
            'x_qe_xc_igcx': 18,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'GGA_X_Q2D',
            }
        ],
        'xc_terms_remove': [
            {
                'XC_functional_name': 'LDA_X',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_igcx_name': 'q2dx',
            'x_qe_xc_igcx_comment': 'Q2D exchange grad corr',
            'x_qe_xc_igcx': 19,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'HYB_GGA_XC_GAU_PBE',
                'exx_compute_weight': lambda exx: 1.0 if (abs(exx) > TOL) else 0.0,
            },
            {
                'XC_functional_name': 'GGA_X_PBE',
                'exx_compute_weight': lambda exx: 0.0 if (abs(exx) > TOL) else 1.0,
            },
        ],
        'xc_terms_remove': [
            {
                'XC_functional_name': 'LDA_X',
            },
            {
                'XC_functional_name': 'GGA_C_PBE',
                'exx_compute_weight': lambda exx: 1.0 if (abs(exx) > TOL) else 0.0,
            },
        ],
        'xc_section_method': {
            'x_qe_xc_igcx_name': 'gaup',
            'x_qe_xc_igcx_comment': 'Gau-PBE hybrid exchange',
            'x_qe_xc_igcx': 20,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'GGA_X_PW86',
            }
        ],
        'xc_terms_remove': [
            {
                'XC_functional_name': 'LDA_X',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_igcx_name': 'pw86',
            'x_qe_xc_igcx_comment': 'Perdew-Wang (1986) exchange',
            'x_qe_xc_igcx': 21,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'GGA_X_B86_MGC',
            }
        ],
        'xc_terms_remove': [
            {
                'XC_functional_name': 'LDA_X',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_igcx_name': 'b86b',
            'x_qe_xc_igcx_comment': 'Becke (1986) exchange',
            'x_qe_xc_igcx': 22,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'GGA_X_OPTB88_VDW',
            }
        ],
        'xc_terms_remove': [
            {
                'XC_functional_name': 'LDA_X',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_igcx_name': 'obk8',
            'x_qe_xc_igcx_comment': 'optB88 exchange',
            'x_qe_xc_igcx': 23,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'GGA_X_OPTB86B_VDW',
            }
        ],
        'xc_terms_remove': [
            {
                'XC_functional_name': 'LDA_X',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_igcx_name': 'ob86',
            'x_qe_xc_igcx_comment': 'optB86b exchange',
            'x_qe_xc_igcx': 24,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'GGA_X_EV93',
            }
        ],
        'xc_terms_remove': [
            {
                'XC_functional_name': 'LDA_X',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_igcx_name': 'evx',
            'x_qe_xc_igcx_comment': 'Engel-Vosko exchange',
            'x_qe_xc_igcx': 25,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'GGA_X_B86_R',
            }
        ],
        'xc_terms_remove': [
            {
                'XC_functional_name': 'LDA_X',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_igcx_name': 'b86r',
            'x_qe_xc_igcx_comment': 'revised Becke (b86b)',
            'x_qe_xc_igcx': 26,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'GGA_X_LV_RPW86',
            }
        ],
        'xc_terms_remove': [
            {
                'XC_functional_name': 'LDA_X',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_igcx_name': 'cx13',
            'x_qe_xc_igcx_comment': 'consistent exchange',
            'x_qe_xc_igcx': 27,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'GGA_X_B88',
                'XC_functional_weight': 0.542,
                'exx_compute_weight': lambda exx: 0.542 if (abs(exx) > TOL) else 1.0,
            },
            {
                'XC_functional_name': 'GGA_X_PW91',
                'XC_functional_weight': 0.167,
                'exx_compute_weight': lambda exx: 0.167 if (abs(exx) > TOL) else 0.0,
            },
        ],
        'xc_terms_remove': [
            {
                'XC_functional_name': 'LDA_X',
                'exx_compute_weight': lambda exx: 0.709 if (abs(exx) > TOL) else 1.0,
            }
        ],
        'xc_section_method': {
            'x_qe_xc_igcx_name': 'x3lp',
            'x_qe_xc_igcx_comment': 'X3LYP (Becke88*0.542  + Perdew-Wang91*0.167)',
            'x_qe_xc_igcx': 28,
        },
    },
    # New additions for qe-v6.5.0: igcx=[29..42]
    # - - - - - -
    # igcx: 29. The ingredient 'vdW-DF-cx' is documented in the nomad-meta-info, where
    # it has the name 'vdw_c_df_cx'
    # 'https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-info/-/wikis/metainfo/XC-functional':
    # 'vdW-DF-cx' implies igcx=27 => 'cx13', hence LDA_X is implicit. Full weight.
    {
        'xc_terms': [
            {
                'XC_functional_name': 'vdw_c_df_cx',
            },
            {
                'XC_functional_name': 'HF_X',
                'exx_compute_weight': lambda exx: exx,
                'XC_functional_weight': 0.25,
            },
        ],
        'xc_terms_remove': [
            {
                'XC_functional_name': 'LDA_X',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_igcx_name': 'cx0',
            'x_qe_xc_igcx_comment': 'vdW-DF-cx+HF/4 (cx13-0)',
            'x_qe_xc_igcx': 29,
        },
    },
    # - - - - - -
    # igcx:30. Needs full LDA_X removal, due to 'GGA_X_RPW86' (see igcx:27)
    {
        'xc_terms': [
            {
                'XC_functional_name': 'GGA_X_RPW86',
            },
            {
                'XC_functional_name': 'HF_X',
                'exx_compute_weight': lambda exx: exx,
                'XC_functional_weight': 0.25,
            },
        ],
        'xc_terms_remove': [
            {
                'XC_functional_name': 'LDA_X',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_igcx_name': 'r860',
            'x_qe_xc_igcx_comment': 'rPW86+HF/4 (rw86-0); (for DF0)',
            'x_qe_xc_igcx': 30,
        },
    },
    # - - - - - -
    # igcx:31. Similar comments as in 'igcx:29'
    {
        'xc_terms': [
            {
                'XC_functional_name': 'vdw_c_df_cx',
            },
            {
                'XC_functional_name': 'HF_X',
                'exx_compute_weight': lambda exx: exx,
                'XC_functional_weight': 0.20,
            },
        ],
        'xc_terms_remove': [
            {
                'XC_functional_name': 'LDA_X',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_igcx_name': 'cx0p',
            'x_qe_xc_igcx_comment': 'vdW-DF-cx+HF/5 (cx13-0p)',
            'x_qe_xc_igcx': 31,
        },
    },
    # - - - - - -
    {
        'xc_terms': [
            {
                'XC_functional_name': 'GGA_X_RESERVED',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_igcx_name': 'ahcx',
            'x_qe_xc_igcx_comment': 'vdW-DF-cx based; not yet in use (reserved PH)',
            'x_qe_xc_igcx': 32,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'GGA_X_RESERVED',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_igcx_name': 'ahf2',
            'x_qe_xc_igcx_comment': 'vdW-DF2 based; not yet in use (reserved PH)',
            'x_qe_xc_igcx': 33,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'GGA_X_RESERVED',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_igcx_name': 'ahpb',
            'x_qe_xc_igcx_comment': 'PBE based; not yet in use (reserved PH)',
            'x_qe_xc_igcx': 34,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'GGA_X_RESERVED',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_igcx_name': 'ahps',
            'x_qe_xc_igcx_comment': 'PBE-sol based; not in use (reserved PH)',
            'x_qe_xc_igcx': 35,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'GGA_X_RESERVED',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_igcx_name': 'cx14',
            'x_qe_xc_igcx_comment': 'Exporations (typo?: explorations), (reserved PH)',
            'x_qe_xc_igcx': 36,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'GGA_X_RESERVED',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_igcx_name': 'cx15',
            'x_qe_xc_igcx_comment': 'Exporations (typo? explorations?)(reserved PH)',
            'x_qe_xc_igcx': 37,
        },
    },
    #
    # igcx': 38. Ingredients:
    #   'b86r' -> 'igcx:26' -> 'GGA_X_B86_R'
    #   'vdW-DF2' -> 'vdw_c_df2' . See nomad's gitlab:
    #   'https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-info/-/wikis/metainfo/XC-functional':
    {
        'xc_terms': [
            {
                'XC_functional_name': 'vdw_c_df2',
            },
            {
                'XC_functional_name': 'GGA_X_B86_R',
            },
            {
                'XC_functional_name': 'HF_X',
                'exx_compute_weight': lambda exx: exx,
                'XC_functional_weight': 0.25,
            },
        ],
        'xc_terms_remove': [
            {
                'XC_functional_name': 'LDA_X',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_igcx_name': 'br0',
            'x_qe_xc_igcx_comment': 'vdW-DF2-b86r+HF/4 (b86r-0)',
            'x_qe_xc_igcx': 38,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'GGA_X_RESERVED',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_igcx_name': 'cx16',
            'x_qe_xc_igcx_comment': 'Exporations (typo?, explorations?)(reserved PH)',
            'x_qe_xc_igcx': 39,
        },
    },
    # - - - -
    {
        'xc_terms': [
            {
                'XC_functional_name': 'vdw_c_df1',
            },
            {
                'XC_functional_name': 'GGA_X_C09X',
            },
            {
                'XC_functional_name': 'HF_X',
                'exx_compute_weight': lambda exx: exx,
                'XC_functional_weight': 0.25,
            },
        ],
        'xc_terms_remove': [
            {
                'XC_functional_name': 'LDA_X',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_igcx_name': 'c090',
            'x_qe_xc_igcx_comment': 'vdW-DF-c09+HF/4 (c09-0)',
            'x_qe_xc_igcx': 40,
        },
    },
    # - - - - - - -
    # 'igcx:41' Note: 'B86b' is defined in 'igcx:22'
    {
        'xc_terms': [
            {
                'XC_functional_name': 'GGA_X_B86_MGC',
            }
        ],
        'xc_terms_remove': [
            {
                'XC_functional_name': 'LDA_X',
                'XC_functional_weight': 0.75,
                'exx_compute_weight': lambda exx: (1.0 - exx),
            }
        ],
        'xc_section_method': {
            'x_qe_xc_igcx_name': 'b86x',
            'x_qe_xc_igcx_comment': 'B86b exchange * 0.75',
            'x_qe_xc_igcx': 41,
        },
    },
    # - - - - - - -
    # 'B88' is defined in 'igcx:1'
    {
        'xc_terms': [
            {
                'XC_functional_name': 'GGA_X_B88',
            }
        ],
        'xc_terms_remove': [
            {
                'XC_functional_name': 'LDA_X',
                'XC_functional_weight': 0.50,
                'exx_compute_weight': lambda exx: (1.0 - exx),
            }
        ],
        'xc_section_method': {
            'x_qe_xc_igcx_name': 'b88x',
            'x_qe_xc_igcx_comment': 'B88 exchange * 0.50',
            'x_qe_xc_igcx': 42,
        },
    },
]

# UNchanged between espresso v5.4 & v6.5
_correlation_gradient_correction_map = [
    None,
    {
        'xc_terms': [
            {
                'XC_functional_name': 'GGA_C_P86',
            }
        ],
        'xc_terms_remove': [
            {
                'XC_functional_name': 'LDA_C_PW',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_igcc_name': 'p86',
            'x_qe_xc_igcc_comment': 'Perdew86',
            'x_qe_xc_igcc': 1,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'GGA_C_PW91',
            }
        ],
        'xc_terms_remove': [
            {
                'XC_functional_name': 'LDA_C_PW',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_igcc_name': 'ggc',
            'x_qe_xc_igcc_comment': 'Perdew-Wang 91 corr.',
            'x_qe_xc_igcc': 2,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'GGA_C_LYP',
            }
        ],
        'xc_terms_remove': [
            {
                'XC_functional_name': 'LDA_C_LYP',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_igcc_name': 'blyp',
            'x_qe_xc_igcc_comment': 'Lee-Yang-Parr',
            'x_qe_xc_igcc': 3,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'GGA_C_PBE',
            }
        ],
        'xc_terms_remove': [
            {
                'XC_functional_name': 'LDA_C_PW',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_igcc_name': 'pbc',
            'x_qe_xc_igcc_comment': 'Perdew-Burke-Ernzenhof corr',
            'x_qe_xc_igcc': 4,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'GGA_XC_HCTH_120',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_igcc_name': 'hcth',
            'x_qe_xc_igcc_comment': 'Cambridge exch, Handy et al, HCTH/120',
            'x_qe_xc_igcc': 5,
        },
    },
    {
        # igcc=6 is not defined in 5.4's funct.f90
        #        definition taken from 5.0, which did not have separate imeta
        'xc_terms': [
            {
                'XC_functional_name': 'MGGA_C_TPSS',
            }
        ],
        'xc_terms_remove': [
            {
                'XC_functional_name': 'LDA_C_PW',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_igcc_name': 'tpss',
            'x_qe_xc_igcc_comment': 'TPSS Meta-GGA (Espresso-version < 5.1)',
            'x_qe_xc_igcc': 6,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'GGA_C_LYP',
                'XC_functional_weight': 0.81,
            }
        ],
        'xc_terms_remove': [
            {
                'XC_functional_name': 'LDA_C_LYP',
                'XC_functional_weight': 0.81,
            }
        ],
        'xc_section_method': {
            'x_qe_xc_igcc_name': 'b3lp',
            'x_qe_xc_igcc_comment': 'B3LYP (Lee-Yang-Parr*0.81)',
            'x_qe_xc_igcc': 7,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'GGA_C_PBE_SOL',
            }
        ],
        'xc_terms_remove': [
            {
                'XC_functional_name': 'LDA_C_PW',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_igcc_name': 'psc',
            'x_qe_xc_igcc_comment': 'PBEsol corr',
            'x_qe_xc_igcc': 8,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'GGA_C_PBE',
            }
        ],
        'xc_terms_remove': [
            {
                'XC_functional_name': 'LDA_C_PW',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_igcc_name': 'pbe',
            'x_qe_xc_igcc_comment': 'same as PBX, back-comp.',
            'x_qe_xc_igcc': 9,
        },
    },
    {
        # igcc=10 is not defined in 5.4's funct.f90
        #        definition taken from 5.0, which did not have separate imeta
        #        functionals.f90 tells that correlation is taken from tpss
        'xc_terms': [
            {
                'XC_functional_name': 'MGGA_C_TPSS',
            }
        ],
        'xc_terms_remove': [
            {
                'XC_functional_name': 'LDA_C_PW',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_igcc_name': 'tb09',
            'x_qe_xc_igcc_comment': 'TB09 Meta-GGA (Espresso-version < 5.1)',
            'x_qe_xc_igcc': 10,
        },
    },
    {
        # igcc=11 is not defined in 5.4's funct.f90
        #        definition taken from 5.0, which did not have separate imeta
        'xc_terms': [
            {
                'XC_functional_name': 'MGGA_C_M06_L',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_igcc_name': 'm6lc',
            'x_qe_xc_igcc_comment': 'M06L Meta-GGA (Espresso-version < 5.1)',
            'x_qe_xc_igcc': 11,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'GGA_C_Q2D',
            }
        ],
        'xc_terms_remove': [
            {
                'XC_functional_name': 'LDA_C_PW',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_igcc_name': 'q2dc',
            'x_qe_xc_igcc_comment': 'Q2D correlation grad corr',
            'x_qe_xc_igcc': 12,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'GGA_C_LYP',
                'XC_functional_weight': 0.871,
            }
        ],
        'xc_terms_remove': [
            {
                'XC_functional_name': 'LDA_C_LYP',
                'XC_functional_weight': 0.871,
            }
        ],
        'xc_section_method': {
            'x_qe_xc_igcc_name': 'x3lp',
            'x_qe_xc_igcc_comment': 'X3LYP (Lee-Yang-Parr*0.871)',
            'x_qe_xc_igcc': 13,
        },
    },
    {
        #  igcc=14 is not defined in NEITHER of v5.1, v6.1, v6.4's Modules/funct.f90
        # 'BEEF-vdW, a GGA with vdW-DF2 type nonlocal correlation'
        'xc_terms': [
            {
                'XC_functional_name': 'GGA_C_BEEF-vdW',
            }
        ],
        'xc_terms_remove': [
            {
                'XC_functional_name': 'LDA_C_PW',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_igcc_name': 'BEEF-vdW',
            'x_qe_xc_igcc_comment': 'libbeef V0.1.1 library',
            'x_qe_xc_igcc': 14,
        },
    },
]

# New additions for espresso-6.5.0: imeta=[4, 5, 6]
_meta_gga_map = [
    None,
    {
        'xc_terms': [
            {
                'XC_functional_name': 'MGGA_X_TPSS',
            },
            {
                'XC_functional_name': 'MGGA_C_TPSS',
            },
        ],
        'xc_terms_remove': [
            {
                'XC_functional_name': 'LDA_X',
            },
            {
                'XC_functional_name': 'LDA_C_PW',
            },
        ],
        'xc_section_method': {
            'x_qe_xc_imeta_name': 'tpss',
            'x_qe_xc_imeta_comment': 'TPSS Meta-GGA',
            'x_qe_xc_imeta': 1,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'MGGA_X_M06_L',
            },
            {
                'XC_functional_name': 'MGGA_C_M06_L',
            },
        ],
        'xc_section_method': {
            'x_qe_xc_imeta_name': 'm6lx',
            'x_qe_xc_imeta_comment': 'M06L Meta-GGA',
            'x_qe_xc_imeta': 2,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'MGGA_X_TB09',
            },
            {
                # confirmed by looking into functionals.f90
                'XC_functional_name': 'MGGA_C_TPSS',
            },
        ],
        'xc_terms_remove': [
            {
                'XC_functional_name': 'LDA_X',
            },
            {
                'XC_functional_name': 'LDA_C_PW',
            },
        ],
        'xc_section_method': {
            'x_qe_xc_imeta_name': 'tb09',
            'x_qe_xc_imeta_comment': 'TB09 Meta-GGA',
            'x_qe_xc_imeta': 3,
        },
    },
    # imeta = [4,5,6] are new espresso-6.5.0/Modules/funct.f90
    {
        'xc_terms': [
            {
                'XC_functional_name': 'MGGA_X_TPSS',
            },
            {
                'XC_functional_name': 'MGGA_C_TPSS',
            },
        ],
        'xc_terms_remove': [
            {
                'XC_functional_name': 'LDA_X',
            },
            {
                'XC_functional_name': 'LDA_C_PW',
            },
        ],
        'xc_section_method': {
            'x_qe_xc_imeta_name': '+meta',
            'x_qe_xc_imeta_comment': 'activate MGGA even without MGGA-XC',
            'x_qe_xc_imeta': 4,
        },
    },
    # - - - - - - - -
    {
        'xc_terms': [
            {
                'XC_functional_name': 'MGGA_X_SCAN',
            },
            {
                'XC_functional_name': 'MGGA_C_SCAN',
            },
        ],
        'xc_terms_remove': [
            {
                'XC_functional_name': 'LDA_X',
            },
            {
                'XC_functional_name': 'LDA_C_PW',
            },
        ],
        'xc_section_method': {
            'x_qe_xc_imeta_name': 'scan',
            'x_qe_xc_imeta_comment': 'SCAN Meta-GGA ',
            'x_qe_xc_imeta': 5,
        },
    },
    # - - - - - - - -
    {
        'xc_terms': [
            {
                'XC_functional_name': 'HYB_MGGA_X_SCAN0',
            },
            {
                'XC_functional_name': 'MGGA_C_SCAN',
            },
        ],
        'xc_terms_remove': [
            {
                'XC_functional_name': 'LDA_X',
            },
            {
                'XC_functional_name': 'LDA_C_PW',
            },
        ],
        'xc_section_method': {
            'x_qe_xc_imeta_name': 'sca0',
            'x_qe_xc_imeta_comment': 'SCAN0  Meta-GGA',
            'x_qe_xc_imeta': 6,
        },
    },
]

_van_der_waals_map = [
    None,
    {
        'xc_terms': [
            {
                'XC_functional_name': 'VDW_XC_DF1',
            }
        ],
        'xc_terms_remove': [
            {
                'XC_functional_name': 'LDA_X',
            },
            {
                'XC_functional_name': 'LDA_C_PW',
            },
        ],
        'xc_section_method': {
            'x_qe_xc_inlc_name': 'vdw1',
            'x_qe_xc_inlc_comment': 'vdW-DF1',
            'x_qe_xc_inlc': 1,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'VDW_XC_DF2',
            }
        ],
        'xc_terms_remove': [
            {
                'XC_functional_name': 'LDA_X',
            },
            {
                'XC_functional_name': 'LDA_C_PW',
            },
        ],
        'xc_section_method': {
            'x_qe_xc_inlc_name': 'vdw2',
            'x_qe_xc_inlc_comment': 'vdW-DF2',
            'x_qe_xc_inlc': 2,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'VDW_C_RVV10',
            }
        ],
        'xc_terms_remove': [
            {
                'XC_functional_name': 'GGA_C_PBE',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_inlc_name': 'vv10',
            'x_qe_xc_inlc_comment': 'rVV10',
            'x_qe_xc_inlc': 3,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'VDW_DFX_x_qe',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_inlc_name': 'vdwx',
            'x_qe_xc_inlc_comment': 'vdW-DF-x (reserved Thonhauser, not implemented)',
            'x_qe_xc_inlc': 4,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'VDW_DFY_x_qe',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_inlc_name': 'vdwy',
            'x_qe_xc_inlc_comment': 'vdW-DF-y (reserved Thonhauser, not implemented)',
            'x_qe_xc_inlc': 5,
        },
    },
    {
        'xc_terms': [
            {
                'XC_functional_name': 'VDW_DFZ_x_qe',
            }
        ],
        'xc_section_method': {
            'x_qe_xc_inlc_name': 'vdwz',
            'x_qe_xc_inlc_comment': 'vdW-DF-z (reserved Thonhauser, not implemented)',
            'x_qe_xc_inlc': 6,
        },
    },
]

libxc_shortcut = {
    '0.810*GGA_C_LYP+0.720*GGA_X_B88+0.200*HF_X+0.190*LDA_C_VWN': {
        'xc_terms': [
            {
                'XC_functional_name': 'HYB_GGA_XC_B3LYP',
            }
        ]
    },
    'GGA_C_PBE+0.750*GGA_X_PBE+0.250*HF_X': {
        'xc_terms': [
            {
                'XC_functional_name': 'HYB_GGA_XC_PBEH',
            }
        ]
    },
}

xc_functional_map = [
    _exchange_map,
    _correlation_map,
    _exchange_gradient_correction_map,
    _correlation_gradient_correction_map,
    _van_der_waals_map,
    _meta_gga_map,
]
