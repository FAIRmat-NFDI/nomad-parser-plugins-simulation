from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass

from nomad.metainfo import SchemaPackage
from nomad_simulations.schema_packages import (
    general,
    model_method,
    model_system,
    outputs,
    properties,
    workflow,
)

from nomad_simulation_parsers.schema_packages.utils import add_mapping_annotation

m_package = SchemaPackage()

TEXT_KEY = 'text'
TEXT_DOS_KEY = 'text_dos'
TEXT_GW_KEY = 'text_gw'
SINGLE_POINT_KEY = 'single_point'
GEO_OPT_WORKFLOW_KEY = 'geo_opt_workflow'
MD_WORKFLOW_KEY = 'md_workflow'

# add_mapping_annotations(EntryArchive.m_def, TEXT_KEY, '@')
# add_mapping_annotations(EntryArchive.m_def, TEXT_DOS_KEY, '@')
# add_mapping_annotations(EntryArchive.m_def, TEXT_GW_KEY, '@')
# add_mapping_annotations(EntryArchive.m_def, SINGLE_POINT_KEY, '@')
# add_mapping_annotations(EntryArchive.m_def, GEO_OPT_WORKFLOW_KEY, '@')
# add_mapping_annotations(EntryArchive.m_def, MD_WORKFLOW_KEY, '@')

add_mapping_annotation(general.Simulation.m_def, TEXT_KEY, '@')
add_mapping_annotation(general.Simulation.m_def, TEXT_DOS_KEY, '@')
add_mapping_annotation(general.Simulation.m_def, TEXT_GW_KEY, '@')


class Simulation(general.Simulation):
    add_mapping_annotation(general.Simulation.program, TEXT_KEY, '.@')
    add_mapping_annotation(
        general.Simulation.model_system,
        TEXT_KEY,
        (
            'get_sections',
            ['.@'],
            dict(include=['lattice_vectors', 'structure', 'sub_structure']),
        ),
    )
    # DFT method
    add_mapping_annotation(model_method.DFT.m_def, TEXT_KEY, '.@')
    # gw method
    add_mapping_annotation(model_method.GW.m_def, TEXT_GW_KEY, '.@')
    # electronic structure outputs
    add_mapping_annotation(
        outputs.Outputs.m_def,
        TEXT_KEY,
        (
            'get_sections',
            ['.@'],
            dict(include=['energy', 'energy_components', 'forces', 'eigenvalues']),
        ),
    )
    add_mapping_annotation(
        outputs.Outputs.m_def,
        TEXT_DOS_KEY,
        (
            'get_sections',
            ['.@'],
            dict(include=['total_dos_files', 'species_projected_dos_files']),
        ),
    )


class Program(general.Program):
    add_mapping_annotation(general.Program.version, TEXT_KEY, '.version')


class DFT(model_method.DFT):
    add_mapping_annotation(model_method.DFT.xc, TEXT_KEY, '.@')


class XCFunctional(model_method.XCFunctional):
    add_mapping_annotation(
        model_method.XCFunctional.components,
        TEXT_KEY,
        ('get_xc_functionals', ['.controlInOut_xc']),
    )


class XCComponent(model_method.XCComponent):
    add_mapping_annotation(model_method.XCComponent.canonical_label, TEXT_KEY, '.name')


class GW(model_method.GW):
    add_mapping_annotation(
        model_method.GW.type, TEXT_GW_KEY, ('get_gw_flag', ['.gw_flag'])
    )


class ModelSystem(model_system.ModelSystem):
    add_mapping_annotation(model_system.AtomicCell.m_def, TEXT_KEY, '.@')
    add_mapping_annotation(
        model_system.ModelSystem.positions,
        TEXT_KEY,
        '.structure.positions',
        unit='angstrom',
    )
    add_mapping_annotation(model_system.AtomsState.m_def, TEXT_KEY, '.structure.labels')


class AtomicCell(model_system.AtomicCell):
    add_mapping_annotation(
        model_system.AtomicCell.lattice_vectors, TEXT_KEY, '.lattice_vectors'
    )


class AtomsState(model_system.AtomsState):
    add_mapping_annotation(model_system.AtomsState.chemical_symbol, TEXT_KEY, '.@')


class Outputs(outputs.Outputs):
    add_mapping_annotation(
        outputs.Outputs.total_energies, TEXT_KEY, ('get_energies', ['.@'])
    )
    add_mapping_annotation(
        outputs.Outputs.total_forces, TEXT_KEY, ('get_forces', ['.@'])
    )
    add_mapping_annotation(
        outputs.Outputs.electronic_eigenvalues,
        TEXT_KEY,
        ('get_eigenvalues', ['.eigenvalues', 'array_size_parameters']),
    )
    add_mapping_annotation(
        outputs.Outputs.electronic_dos,
        TEXT_DOS_KEY,
        (
            'get_dos',
            [
                '.total_dos_files',
                '.atom_projected_dos_files',
                '.species_projected_dos_files',
            ],
        ),
    )


class TotalEnergy(properties.energies.TotalEnergy):
    add_mapping_annotation(properties.energies.TotalEnergy.value, TEXT_KEY, '.value')
    add_mapping_annotation(
        properties.energies.TotalEnergy.contributions, TEXT_KEY, '.components'
    )


class BaseEnergy(properties.energies.BaseEnergy):
    add_mapping_annotation(properties.energies.BaseEnergy.name, TEXT_KEY, '.name')


class TotalForce(properties.forces.TotalForce):
    add_mapping_annotation(properties.forces.TotalForce.value, TEXT_KEY, '.forces')


class ElectronicEigenvalues(properties.ElectronicEigenvalues):
    add_mapping_annotation(
        properties.ElectronicEigenvalues.n_levels, TEXT_KEY, '.nbands'
    )
    add_mapping_annotation(
        properties.ElectronicEigenvalues.value, TEXT_KEY, '.eigenvalues'
    )
    add_mapping_annotation(
        properties.ElectronicEigenvalues.occupation, TEXT_KEY, '.occupations'
    )


class DOSProfile(properties.spectral_profile.DOSProfile):
    ### dos quantities
    add_mapping_annotation(
        properties.spectral_profile.DOSProfile.value, TEXT_DOS_KEY, '.values'
    )


class ElectronicDensityOfStates(properties.spectral_profile.ElectronicDensityOfStates):
    add_mapping_annotation(
        properties.spectral_profile.ElectronicDensityOfStates.value,
        TEXT_DOS_KEY,
        '.values',
    )
    ### projected dos
    add_mapping_annotation(
        properties.spectral_profile.ElectronicDensityOfStates.projected_dos,
        TEXT_DOS_KEY,
        '.projected_dos',
    )


class SimulationWorkflow(workflow.general.SimulationWorkflow):
    # TODO find a more elegant fix to not parse tasks recursively, this will be filled
    # in by workflow normalizer from outputs
    add_mapping_annotation(
        workflow.general.SimulationWorkflow.tasks, GEO_OPT_WORKFLOW_KEY, '.tasks'
    )
    add_mapping_annotation(
        workflow.general.SimulationWorkflow.tasks, MD_WORKFLOW_KEY, '.tasks'
    )


# workflow
add_mapping_annotation(workflow.single_point.SinglePoint.m_def, SINGLE_POINT_KEY, '.@')
# geometry optimization workflow
add_mapping_annotation(
    workflow.geometry_optimization.GeometryOptimization.m_def,
    GEO_OPT_WORKFLOW_KEY,
    '.@',
)
# molecular dynamics workflow
add_mapping_annotation(
    workflow.molecular_dynamics.MolecularDynamics.m_def, MD_WORKFLOW_KEY, '.@'
)


class MolecularDynamics(workflow.MolecularDynamics):
    add_mapping_annotation(
        workflow.molecular_dynamics.MolecularDynamicsMethod.m_def, MD_WORKFLOW_KEY, '.@'
    )
    add_mapping_annotation(
        workflow.molecular_dynamics.MolecularDynamicsResults.m_def,
        MD_WORKFLOW_KEY,
        '.@',
    )


class GeometryOptimization(workflow.GeometryOptimization):
    add_mapping_annotation(
        workflow.geometry_optimization.GeometryOptimizationMethod.m_def,
        GEO_OPT_WORKFLOW_KEY,
        '.@',
    )
    add_mapping_annotation(
        workflow.geometry_optimization.GeometryOptimizationResults.m_def,
        GEO_OPT_WORKFLOW_KEY,
        '.@',
    )


class GeometryOptimizationMethod(
    workflow.geometry_optimization.GeometryOptimizationMethod
):
    add_mapping_annotation(
        workflow.geometry_optimization.GeometryOptimizationMethod.optimization_method,
        GEO_OPT_WORKFLOW_KEY,
        '.geometry_relaxation_method',
    )


class MolecularDynamicsMethod(workflow.molecular_dynamics.MolecularDynamicsMethod):
    add_mapping_annotation(
        workflow.molecular_dynamics.MolecularDynamicsMethod.integration_timestep,
        MD_WORKFLOW_KEY,
        '.control_inout.md_timestep',
    )
    add_mapping_annotation(
        workflow.molecular_dynamics.MolecularDynamicsMethod.thermodynamic_ensemble,
        MD_WORKFLOW_KEY,
        '.control_inout.md_run[0].ensemble',
    )


class MolecularDynamicsResults(workflow.molecular_dynamics.MolecularDynamicsResults):
    pass


add_mapping_annotation(
    workflow.molecular_dynamics.MolecularDynamicsResults.temperatures,
    MD_WORKFLOW_KEY,
    'molecular_dynamics[*].md_calculation_info."Temperature (nuclei)"',
)


try:
    m_package.__init_metainfo__()
except Exception:
    pass
