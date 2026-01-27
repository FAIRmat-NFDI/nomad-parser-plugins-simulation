from nomad.metainfo import SchemaPackage
from nomad_simulations.schema_packages import (
    atoms_state,
    general,
    model_method,
    model_system,
    numerical_settings,
    outputs,
    properties,
    workflow,
)

from nomad_simulation_parsers.schema_packages.utils import add_mapping_annotation

m_package = SchemaPackage()

INFO_KEY = 'info'
INPUT_XML_KEY = 'input_xml'
EIGVAL_KEY = 'eigval'
BANDSTRUCTURE_XML_KEY = 'bandstructure_xml'
DOS_XML_KEY = 'dos_xml'

# simulation
add_mapping_annotation(general.Simulation.m_def, INFO_KEY, '@')
add_mapping_annotation(general.Simulation.m_def, INPUT_XML_KEY, '@')
add_mapping_annotation(general.Simulation.m_def, EIGVAL_KEY, '@')
add_mapping_annotation(general.Simulation.m_def, BANDSTRUCTURE_XML_KEY, '@')
add_mapping_annotation(general.Simulation.m_def, DOS_XML_KEY, '@')

# geometry optimization

class Simulation(general.Simulation):
    add_mapping_annotation(general.Simulation.program, INFO_KEY, '.@')
    # DFT method
    add_mapping_annotation(
        model_method.DFT.m_def, INFO_KEY, '.initialization.xc_functional'
    )
    add_mapping_annotation(model_method.DFT.m_def, INPUT_XML_KEY, '.input.groundstate')
    add_mapping_annotation(model_method.DFT.m_def, BANDSTRUCTURE_XML_KEY, '.@')
    add_mapping_annotation(
        general.Simulation.model_system,
        INFO_KEY,
        ('get_configurations', ['.@']),
        cache=True,
    )
    add_mapping_annotation(
        general.Simulation.outputs, INFO_KEY, ('get_configurations', ['.@'])
    )
    add_mapping_annotation(general.Simulation.outputs, EIGVAL_KEY, '.@')
    add_mapping_annotation(general.Simulation.outputs, BANDSTRUCTURE_XML_KEY, '.@')
    add_mapping_annotation(general.Simulation.outputs, DOS_XML_KEY, '.@')


class Program(general.Program):
    add_mapping_annotation(general.Program.version, INFO_KEY, '.program_version')

workflow.single_point.SinglePointMethod.m_def.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(dict(info=Mapper(mapper='.@')))

workflow.SinglePoint.m_def.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(dict(info=Mapper(mapper='.@')))

workflow.single_point.SinglePointMethod.convergence.m_annotations.setdefault( 
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(info=Mapper(mapper=('get_single_point_convergence', ['.@']))))


# general workflow convergence mapping

workflow.general.WorkflowConvergenceTarget.convergence_parameter_name.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(info=Mapper(mapper='.convergence_parameter_name'),
                  geo_opt=Mapper(mapper='.convergence_parameter_name')))

workflow.general.WorkflowConvergenceTarget.convergence_threshold.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(info=Mapper(mapper='.convergence_threshold'),
                  geo_opt=Mapper(mapper='.convergence_threshold')))

workflow.general.WorkflowConvergenceTarget.threshold_type.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(info=Mapper(mapper='.threshold_type'),
                  geo_opt=Mapper(mapper='.threshold_type')))

workflow.general.WorkflowConvergenceTarget.threshold_unit.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(info=Mapper(mapper='.convergence_threshold_unit'),
                  geo_opt=Mapper(mapper='.convergence_threshold_unit')))

# outputs
general.Simulation.outputs.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    info=Mapper(mapper='.@'),
    eigval=Mapper(mapper='.@'),
    bandstructure_xml=Mapper(mapper='.@'),
    dos_xml=Mapper(mapper='.@'),

)
# outputs.Outputs.total_energies.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
#     info=Mapper(mapper='.@')
# )
# outputs.Outputs.m_def.m_annotations[MAPPING_ANNOTATION_KEY]=dict(
#     info=Mapper(mapper='.@')
# )
outputs.Outputs.total_forces.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    info=Mapper(mapper=('get_forces', ['@']))
)
outputs.Outputs.electronic_eigenvalues.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    eigval=Mapper(mapper=('get_eigenvalues', ['.@']))
)
# outputs.Outputs.electronic_band_structures.m_annotations[MAPPING_ANNOTATION_KEY]=dict(
#     bandstructure_xml=Mapper(mapper=('get_bandstructures', ['.@']))
# )
# outputs.Outputs.electronic_dos.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
#     dos_xml=Mapper(mapper=('dos.totaldos.diagram'))
# )
outputs.Outputs.scf_steps.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    info=Mapper(mapper=('get_scf_steps', ['@']))
)

# output quantities
# outputs.SCFSteps.m_def.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
#     info=Mapper(mapper=('get_scf_steps', ['.@']))
# )
outputs.SCFSteps.durations.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    info=Mapper(mapper=('.durations'))
)
outputs.SCFSteps.energies_total.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    info=Mapper(mapper=('.energies_total'))
)
outputs.SCFSteps.delta_energies_total.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    info=Mapper(mapper=('.delta_energies_total'))
)
outputs.SCFSteps.delta_potential_rms.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    info=Mapper(mapper=('.delta_potential_rms'))
)
outputs.SCFSteps.delta_density_rms.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    info=Mapper(mapper=('.delta_density_rms'))
)
outputs.SCFSteps.delta_force_abs.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    info=Mapper(mapper=('.delta_force_abs'))
)


# class Simulation(general.Simulation):
#     general.Simulation.program.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
#         info=Mapper(mapper='.@')
#     )
#     # DFT method
#     model_method.DFT.m_def.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
#         info=Mapper(mapper='.initialization.xc_functional'),
#         input_xml=Mapper(mapper='.input.groundstate'),
#         bandstructure_xml=Mapper(mapper='.@'),
#     )
#     general.Simulation.model_system.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
#         info=Mapper(mapper=('get_configurations', ['.@']), cache=True)
#     )
#     general.Simulation.outputs.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
#         #info=Mapper(mapper=('get_configurations', ['.@'])),
#         eigval=Mapper(mapper='.@'),
#         bandstructure_xml=Mapper(mapper='.@'),
#         dos_xml=Mapper(mapper='.@'),
#     )


general.Program.m_def.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    info=Mapper(mapper='.@')
)

general.Program.version.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    info=Mapper(mapper='.program_version')
)

# properties

properties.forces.TotalForce.value.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    info=Mapper(mapper='.forces')
)

class ModelMethod(model_method.ModelMethod):
    add_mapping_annotation(numerical_settings.KSpace.m_def, BANDSTRUCTURE_XML_KEY, '.@')


class KSpace(numerical_settings.KSpace):
    add_mapping_annotation(
        numerical_settings.KSpace.k_line_path, BANDSTRUCTURE_XML_KEY, '.@'
    )


class KLinePath(numerical_settings.KLinePath):
    add_mapping_annotation(
        numerical_settings.KLinePath.high_symmetry_path_names,
        BANDSTRUCTURE_XML_KEY,
        r'bandstructure.vertex[*]."@label"',
    )
    add_mapping_annotation(
        numerical_settings.KLinePath.high_symmetry_path_values,
        BANDSTRUCTURE_XML_KEY,
        ('reshape_coords', [r'bandstructure.vertex[*]."@coord"']),
    )


class DFT(model_method.DFT):
    add_mapping_annotation(model_method.DFT.xc, INFO_KEY, '.@')
    add_mapping_annotation(model_method.DFT.xc, INPUT_XML_KEY, '.@')


class XCFunctional(model_method.XCFunctional):
    add_mapping_annotation(
        model_method.XCFunctional.components,
        INFO_KEY,
        ('get_xc_functionals', ['.type']),
    )
    add_mapping_annotation(
        model_method.XCFunctional.components,
        INPUT_XML_KEY,
        ('get_xc_functionals', ['.libxc']),
    )


class XCComponent(model_method.XCComponent):
    add_mapping_annotation(model_method.XCComponent.canonical_label, INFO_KEY, '.libxc')
    add_mapping_annotation(
        model_method.XCComponent.canonical_label, INPUT_XML_KEY, '.libxc'
    )


class ModelSystem(model_system.ModelSystem):
    add_mapping_annotation(model_system.Representation.m_def, INFO_KEY, '.@')
    add_mapping_annotation(model_system.ModelSystem.positions, INFO_KEY, '.positions')
    add_mapping_annotation(model_system.AtomsState.m_def, INFO_KEY, '.atoms')


class Representation(model_system.Representation):
    add_mapping_annotation(
        model_system.Representation.lattice_vectors, INFO_KEY, '.lattice_vectors'
    )


class AtomsState(atoms_state.AtomsState):
    add_mapping_annotation(atoms_state.AtomsState.chemical_symbol, INFO_KEY, '.symbol')


class Outputs(outputs.Outputs):
    add_mapping_annotation(outputs.Outputs.total_energies, INFO_KEY, '.@')
    add_mapping_annotation(
        outputs.Outputs.total_forces, INFO_KEY, ('get_forces', ['.@'])
    )
    add_mapping_annotation(
        outputs.Outputs.electronic_eigenvalues, EIGVAL_KEY, ('get_eigenvalues', ['.@'])
    )
    add_mapping_annotation(
        outputs.Outputs.electronic_band_structures,
        BANDSTRUCTURE_XML_KEY,
        ('get_bandstructures', ['.@']),
    )
    add_mapping_annotation(
        outputs.Outputs.electronic_dos, DOS_XML_KEY, 'dos.totaldos.diagram'
    )


class TotalEnergy(properties.TotalEnergy):
    add_mapping_annotation(
        properties.TotalEnergy.value, INFO_KEY, '.final.energy_total || energy_total'
    )


class TotalForce(properties.forces.TotalForce):
    add_mapping_annotation(properties.forces.TotalForce.value, INFO_KEY, '.forces')


# TODO: check whether this section is k-dependent
class ElectronicEigenvalues(outputs.ElectronicEigenvalues):
    add_mapping_annotation(
        outputs.ElectronicEigenvalues.n_levels, EIGVAL_KEY, '.n_states'
    )
    add_mapping_annotation(
        outputs.ElectronicEigenvalues.value, EIGVAL_KEY, '.eigenvalues'
    )
    add_mapping_annotation(
        outputs.ElectronicEigenvalues.occupation, EIGVAL_KEY, '.occupancies'
    )


class ElectronicBandStructure(outputs.ElectronicBandStructure):
    add_mapping_annotation(
        outputs.ElectronicBandStructure.n_levels, BANDSTRUCTURE_XML_KEY, '.n_states'
    )
    add_mapping_annotation(
        outputs.ElectronicBandStructure.value, BANDSTRUCTURE_XML_KEY, '.energies'
    )


class ElectronicDensityOfStates(outputs.ElectronicDensityOfStates):
    ###### TODO read unit from axis
    add_mapping_annotation(
        outputs.ElectronicDensityOfStates.value,
        DOS_XML_KEY,
        ('to_float', [r'.point[*]."@dos"']),
        unit='1/hartree',
    )
    add_mapping_annotation(
        outputs.ElectronicDensityOfStates.projected_dos,
        DOS_XML_KEY,
        'dos.partialdos.diagram',
    )


try:
    m_package.__init_metainfo__()
except Exception:
    pass
