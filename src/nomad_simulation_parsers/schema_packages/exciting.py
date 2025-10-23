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

#TODO Bring back the function that sets these annotations. Implement it as a class with a call method to be able to set defaults.  # noqa: E501


# simulation
add_mapping_annotation(general.Simulation.m_def, INFO_KEY, '@')
add_mapping_annotation(general.Simulation.m_def, INPUT_XML_KEY, '@')
add_mapping_annotation(general.Simulation.m_def, EIGVAL_KEY, '@')
add_mapping_annotation(general.Simulation.m_def, BANDSTRUCTURE_XML_KEY, '@')
add_mapping_annotation(general.Simulation.m_def, DOS_XML_KEY, '@')

# geometry optimization

workflow.geometry_optimization.GeometryOptimizationModel.m_def.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(dict(info=Mapper(mapper='@')))

workflow.GeometryOptimization.m_def.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(dict(info=Mapper(mapper='@')))
"""
workflow.geometry_optimization.GeometryOptimizationModel.optimization_method.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(out=Mapper(mapper=('get_workflow_method', []))))
"""

workflow.geometry_optimization.GeometryOptimizationModel.convergence.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(info=Mapper(mapper=('get_geometry_convergence', ['.@']))))

workflow.general.WorkflowConvergenceTarget.convergence_parameter_name.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(info=Mapper(mapper='.convergence_parameter_name')))

workflow.general.WorkflowConvergenceTarget.convergence_threshold.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(info=Mapper(mapper='.convergence_threshold')))

workflow.general.WorkflowConvergenceTarget.threshold_type.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(info=Mapper(mapper='.threshold_type')))

workflow.general.WorkflowConvergenceTarget.convergence_threshold_unit.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(info=Mapper(mapper='.convergence_threshold_unit')))

workflow.general.WorkflowConvergenceTarget.is_reached.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(info=Mapper(mapper='.is_reached')))



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
