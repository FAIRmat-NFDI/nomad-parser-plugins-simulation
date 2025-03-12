import re

import numpy as np
from nomad.parsing.file_parser.text_parser import Quantity, TextParser

RE_FLOAT = r'[-+]?\d+\.\d*(?:[Ee][-+]\d+)?'
RE_N = r'[\n\r]'


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
        r'(?:Number of processors in use:|, running on)'
        r'\s*(\d+)\s*(?:processors)*',
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
        rf'absolute magnetization\s*=\s*({RE_FLOAT})'
        rf'\s*Bohr mag/cell',
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
        r'(Cartesian axes\s*site n\.\s*atom\s*positions[\s\S]+?)\n\s*\n',
        repeatas=False,
        sub_parser=TextParser(
            quantities=[
                Quantity('axes', r'(\w+)\s*axes'),
                Quantity('units', r'site n\.\s*atom\s*positions\s*\(\s*(\S+)'),
                Quantity('labels', r'(\w+)\s*tau', repeats=True),
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
