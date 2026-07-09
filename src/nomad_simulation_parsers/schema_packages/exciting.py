from nomad.metainfo import SchemaPackage
from nomad_simulations.schema_packages import (
    atoms_state,
    general,
    model_method,
    model_system,
    numerical_settings,
    outputs,
    properties,
    variables,
    workflow,
)

from nomad_simulation_parsers.schema_packages.utils import add_mapping_annotation

m_package = SchemaPackage()

INFO_KEY = 'info'
INPUT_XML_KEY = 'input_xml'
EIGVAL_KEY = 'eigval'
BANDSTRUCTURE_XML_KEY = 'bandstructure_xml'
DOS_XML_KEY = 'dos_xml'
GEO_OPT_KEY = 'geo_opt'


def add_mapping_annotations(*args):
    for mapping_tuple in args:
        add_mapping_annotation(*mapping_tuple)


# TODO Use this structure?:
# add_mapping_annotations(
#     model_method.DFT.m_def,
#     [
#         (INFO_KEY, '.initialization.xc_functional'),
#         (INPUT_XML_KEY, '.input.groundstate'),
#         (BANDSTRUCTURE_XML_KEY, '.@')
#     ],
# )

# simulation

add_mapping_annotations(
    (general.Simulation.m_def, INFO_KEY, '@'),
    (general.Simulation.m_def, INPUT_XML_KEY, '@'),
    (general.Simulation.m_def, EIGVAL_KEY, '@'),
    (general.Simulation.m_def, BANDSTRUCTURE_XML_KEY, '@'),
    (general.Simulation.m_def, DOS_XML_KEY, '@'),
    (general.Simulation.m_def, GEO_OPT_KEY, '@'),
)

add_mapping_annotation(
    general.Simulation.model_system,
    INFO_KEY,
    ('get_configurations', ['.@']),
    cache=True,
)
add_mapping_annotations(
    (general.Simulation.program, INFO_KEY, '.@'),
    (general.Simulation.outputs, INFO_KEY, '.@'),
    (general.Simulation.outputs, EIGVAL_KEY, '.@'),
    (general.Simulation.outputs, BANDSTRUCTURE_XML_KEY, '.@'),
    (general.Simulation.outputs, DOS_XML_KEY, '.@'),
)

# program

add_mapping_annotation(general.Program.version, INFO_KEY, '.program_version')


# model method
add_mapping_annotations(
    (model_method.DFT.m_def, INFO_KEY, '.initialization.xc_functional'),
    (model_method.DFT.m_def, INPUT_XML_KEY, '.input.groundstate'),
    (model_method.DFT.m_def, BANDSTRUCTURE_XML_KEY, '.@'),
    (model_method.DFT.xc, INFO_KEY, '.@'),
    (model_method.DFT.xc, INPUT_XML_KEY, '.@'),
    (model_method.XCFunctional.components, INFO_KEY, ('get_xc_functionals', ['.type'])),
    (
        model_method.XCFunctional.components,
        INPUT_XML_KEY,
        ('get_xc_functionals', ['.libxc']),
    ),
    (model_method.XCComponent.canonical_label, INFO_KEY, '.libxc'),
    (model_method.XCComponent.canonical_label, INPUT_XML_KEY, '.libxc'),
)

# numerical_settings
add_mapping_annotations(
    (numerical_settings.KSpace.m_def, BANDSTRUCTURE_XML_KEY, '.@'),
    (numerical_settings.KSpace.k_line_path, BANDSTRUCTURE_XML_KEY, '.@'),
    (
        numerical_settings.KLinePath.high_symmetry_path_names,
        BANDSTRUCTURE_XML_KEY,
        r'bandstructure.vertex[*]."@label"',
    ),
    (
        numerical_settings.KLinePath.high_symmetry_path_values,
        BANDSTRUCTURE_XML_KEY,
        ('reshape_coords', [r'bandstructure.vertex[*]."@coord"']),
    ),
)

# model system
add_mapping_annotations(
    (model_system.Representation.m_def, INFO_KEY, '.@'),
    (model_system.ModelSystem.positions, INFO_KEY, '.positions'),
    (model_system.AtomsState.m_def, INFO_KEY, '.atoms'),
    (model_system.Representation.lattice_vectors, INFO_KEY, '.lattice_vectors'),
    (
        model_system.Representation.periodic_boundary_conditions,
        INFO_KEY,
        '.periodic_boundary_conditions',
    ),
)

# atoms state
add_mapping_annotation(atoms_state.AtomsState.chemical_symbol, INFO_KEY, '.symbol')

# properties
add_mapping_annotations(
    (properties.TotalEnergy.value, INFO_KEY, '.final.energy_total || energy_total'),
    (properties.forces.TotalForce.value, INFO_KEY, '.forces'),
)

# outputs
add_mapping_annotations(
    (outputs.Outputs.total_energies, INFO_KEY, ('get_energies', ['.@'])),
    (outputs.Outputs.total_forces, INFO_KEY, ('get_forces', ['.@'])),
    (outputs.Outputs.electronic_eigenvalues, EIGVAL_KEY, ('get_eigenvalues', ['.@'])),
    (outputs.Outputs.electronic_band_gaps, EIGVAL_KEY, ('get_band_gaps', ['.@'])),
    (
        outputs.Outputs.electronic_band_structures,
        BANDSTRUCTURE_XML_KEY,
        ('get_bandstructures', ['.@']),
    ),
    (outputs.Outputs.electronic_dos, DOS_XML_KEY, 'dos.totaldos.diagram'),
)

# eigenvalues
add_mapping_annotations(
    # TODO: check whether this section is k-dependent
    (outputs.ElectronicEigenvalues.n_levels, EIGVAL_KEY, '.n_states'),
    (outputs.ElectronicEigenvalues.value, EIGVAL_KEY, '.eigenvalues'),
    (outputs.ElectronicEigenvalues.occupation, EIGVAL_KEY, '.occupancies'),
)

# bandstructure
add_mapping_annotations(
    (outputs.ElectronicBandStructure.n_levels, BANDSTRUCTURE_XML_KEY, '.n_states'),
    (outputs.ElectronicBandStructure.value, BANDSTRUCTURE_XML_KEY, '.energies'),
    (outputs.ElectronicBandStructure.k_path, BANDSTRUCTURE_XML_KEY, '.k_path'),
    (
        outputs.ElectronicDensityOfStates.projected_dos,
        DOS_XML_KEY,
        'dos.partialdos.diagram',
    ),
)


add_mapping_annotations(
    (outputs.ElectronicBandGap.spin_channel, EIGVAL_KEY, '.spin_channel'),
)
add_mapping_annotation(
    outputs.ElectronicBandGap.value,
    EIGVAL_KEY,
    '.value',
    unit='eV',
)

add_mapping_annotation(
    variables.Energy2.points,
    DOS_XML_KEY,
    ('to_float', [r'.point[*]."@e"']),
    unit='hartree',
)

add_mapping_annotation(variables.Energy2.m_def, DOS_XML_KEY, '.@')

###### TODO read unit from axis
add_mapping_annotation(
    outputs.ElectronicDensityOfStates.value,
    DOS_XML_KEY,
    ('to_float', [r'.point[*]."@dos"']),
    unit='1/hartree',
)

# workflow

add_mapping_annotations(
    (workflow.GeometryOptimization.m_def, GEO_OPT_KEY, '@'),
    (
        workflow.geometry_optimization.GeometryOptimizationMethod.m_def,
        GEO_OPT_KEY,
        '.@',
    ),
    # TODO: Mapping annotations don't work for convergence targets because
    # parser methods return fully-formed metainfo objects, not dictionaries.
    # The mapper expects dict data. Convergence targets are now populated
    # manually in the parser. Consider refactoring the mapping annotation
    # system to support object instantiation or keep manual approach.
    # (workflow.geometry_optimization.GeometryOptimizationMethod.convergence_targets,
    #  GEO_OPT_KEY, ('get_geometry_convergence', ['.@'])),
    # (workflow.geometry_optimization
    #  .GeometryOptimizationMethod.single_point_convergence_targets,
    #  GEO_OPT_KEY, ('get_single_point_convergence', ['.@'])),
    (workflow.single_point.SinglePointMethod.m_def, INFO_KEY, '.@'),
    (workflow.SinglePoint.m_def, INFO_KEY, '.@'),
    # TODO: Same issue as above - manually populated in parser
    # (workflow.single_point.SinglePointMethod.convergence_targets,
    #  INFO_KEY, ('get_single_point_convergence', ['.@']))
)

# scf steps
add_mapping_annotations(
    (outputs.Outputs.scf_steps, INFO_KEY, ('get_scf_steps', ['@'])),
    (outputs.SCFSteps.durations, INFO_KEY, '.durations'),
    (outputs.SCFSteps.energies_total, INFO_KEY, '.energies_total'),
    (outputs.SCFSteps.delta_energies_total, INFO_KEY, '.delta_energies_total'),
    (outputs.SCFSteps.delta_potential_rms, INFO_KEY, '.delta_potential_rms'),
    (outputs.SCFSteps.delta_density_rms, INFO_KEY, '.delta_density_rms'),
    (outputs.SCFSteps.delta_force_abs, INFO_KEY, '.delta_force_abs'),
)

try:
    m_package.__init_metainfo__()
except Exception:
    pass
