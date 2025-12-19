from nomad.datamodel.metainfo.annotations import Mapper
from nomad.metainfo import SchemaPackage
from nomad.parsing.file_parser.mapping_parser import MAPPING_ANNOTATION_KEY
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

m_package = SchemaPackage()

# simulation
general.Simulation.m_def.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    info=Mapper(mapper='@'),
    geo_opt=Mapper(mapper='@'),
    input_xml=Mapper(mapper='@'),
    eigval=Mapper(mapper='@'),
    bandstructure_xml=Mapper(mapper='@'),
    dos_xml=Mapper(mapper='@'),
)

# geometry optimization

workflow.GeometryOptimization.m_def.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(dict(geo_opt=Mapper(mapper='.@')))

workflow.geometry_optimization.GeometryOptimizationMethod.m_def.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(dict(geo_opt=Mapper(mapper='.@')))

"""
workflow.geometry_optimization.GeometryOptimizationMethod.optimization_method.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(out=Mapper(mapper=('get_workflow_method', []))))
"""

workflow.geometry_optimization.GeometryOptimizationMethod.convergence.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(geo_opt=Mapper(mapper=('get_geometry_convergence', ['.@']))))

# workflow.geometry_optimization.GeometryOptimizationMethod.single_point_workflows.m_annotations.setdefault(
#     MAPPING_ANNOTATION_KEY, {}
# ).update(geo_opt=Mapper(mapper='@'))

workflow.geometry_optimization.GeometryOptimizationMethod.single_point_convergence.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(geo_opt=Mapper(mapper=('get_single_point_convergence', ['.@']))))


# single point

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
    numerical_settings.KSpace.m_def.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
        bandstructure_xml=Mapper(mapper='.@')
    )


class KSpace(numerical_settings.KSpace):
    numerical_settings.KSpace.k_line_path.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
        bandstructure_xml=Mapper(mapper='.@')
    )


class KLinePath(numerical_settings.KLinePath):
    numerical_settings.KLinePath.high_symmetry_path_names.m_annotations[
        MAPPING_ANNOTATION_KEY
    ] = dict(bandstructure_xml=Mapper(mapper=r'bandstructure.vertex[*]."@label"'))
    numerical_settings.KLinePath.high_symmetry_path_values.m_annotations[
        MAPPING_ANNOTATION_KEY
    ] = dict(
        bandstructure_xml=Mapper(
            mapper=('reshape_coords', [r'bandstructure.vertex[*]."@coord"'])
        )
    )


#class DFT(model_method.DFT):
#    model_method.DFT.xc_functionals.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
#        info=Mapper(mapper=('get_xc_functionals', ['.type'])),
#        input_xml=Mapper(mapper=('get_xc_functionals', ['.libxc'])),
#    )


# class XCFunctional(model_method.XCFunctional):
#     model_method.XCFunctional.libxc_name.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
#         info=Mapper(mapper='.libxc'), input_xml=Mapper(mapper='.libxc')
#     )


class ModelSystem(model_system.ModelSystem):
    model_system.AtomicCell.m_def.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
        info=Mapper(mapper='.@')
    )
    model_system.ModelSystem.positions.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
        info=Mapper(mapper='.positions')
    )
    model_system.AtomsState.m_def.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
        info=Mapper(mapper='.atoms')
    )


class AtomicCell(model_system.AtomicCell):
    model_system.AtomicCell.lattice_vectors.m_annotations[MAPPING_ANNOTATION_KEY] = (
        dict(info=Mapper(mapper='.lattice_vectors'))
    )


class AtomsState(atoms_state.AtomsState):
    atoms_state.AtomsState.chemical_symbol.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
        info=Mapper(mapper='.symbol')
    )


class TotalEnergy(properties.TotalEnergy):
    properties.TotalEnergy.value.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
        info=Mapper(mapper='.final.energy_total || energy_total')
    )




class ElectronicEigenvalues(outputs.ElectronicEigenvalues):
    # outputs.ElectronicEigenvalues.n_bands.m_annotations[MAPPING_ANNOTATION_KEY] = dict(  # noqa: E501
    #     eigval=Mapper(mapper='.n_states')
    # )
    # we need to use setdefault when annot was defined previously
    outputs.ElectronicEigenvalues.value.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
        eigval=Mapper(mapper='.eigenvalues'),
    )
    outputs.ElectronicEigenvalues.occupation.m_annotations[MAPPING_ANNOTATION_KEY] = (
        dict(eigval=Mapper(mapper='.occupancies'))
    )


class ElectronicBandStructure(outputs.ElectronicBandStructure):
    # outputs.ElectronicBandStructure.n_bands.m_annotations.setdefault(
    #     MAPPING_ANNOTATION_KEY, {}
    # ).update(dict(bandstructure_xml=Mapper(mapper='.n_states')))
    outputs.ElectronicBandStructure.value.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(bandstructure_xml=Mapper(mapper='.energies')))


class ElectronicDensityOfStates(outputs.ElectronicDensityOfStates):
    ###### TODO read unit from axis
    outputs.ElectronicDensityOfStates.value.m_annotations[MAPPING_ANNOTATION_KEY] = (
        dict(
            dos_xml=Mapper(mapper=('to_float', [r'.point[*]."@dos"']), unit='1/hartree')
        )
    )
    outputs.ElectronicDensityOfStates.projected_dos.m_annotations[
        MAPPING_ANNOTATION_KEY
    ] = dict(dos_xml=Mapper(mapper='dos.partialdos.diagram'))


try:
    m_package.__init_metainfo__()
except Exception:
    pass
