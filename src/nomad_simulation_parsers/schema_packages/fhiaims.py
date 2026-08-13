from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass


from nomad.datamodel.metainfo.annotations import Mapper
from nomad.metainfo import SchemaPackage
from nomad.parsing.file_parser.mapping_parser import MAPPING_ANNOTATION_KEY
from nomad_simulations.schema_packages import (
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

TEXT_KEY = 'fhiaims_text'
TEXT_DOS_KEY = 'fhiaims_text_dos'
TEXT_GW_KEY = 'fhiaims_text_gw'
SINGLE_POINT_KEY = 'fhiaims_single_point'
GEO_OPT_WORKFLOW_KEY = 'fhiaims_geo_opt_workflow'
MD_WORKFLOW_KEY = 'fhiaims_md_workflow'

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
    # DFT method - only annotate DFT.m_def, not ModelMethod base class
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
            dict(
                include=[
                    'energy',
                    'energy_components',
                    'forces',
                    'eigenvalues',
                    'humo',
                    'lumo',
                    'self_consistency',
                ]
            ),
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


class SelfConsistency(numerical_settings.SelfConsistency):
    add_mapping_annotation(
        numerical_settings.SelfConsistency.threshold_change,
        TEXT_KEY,
        '.threshold_change',
    )
    add_mapping_annotation(numerical_settings.SelfConsistency.name, TEXT_KEY, '.name')

    add_mapping_annotation(
        numerical_settings.SelfConsistency.n_max_iterations,
        TEXT_KEY,
        '.n_max_iterations',
    )


class DFT(model_method.DFT):
    add_mapping_annotation(
        numerical_settings.SelfConsistency.m_def,
        TEXT_KEY,
        ('get_all_criteria', ['.@']),
        update_mode='append',
    )
    add_mapping_annotation(numerical_settings.KSpace.m_def, TEXT_KEY, '.@')
    # Materialize the `xc` subsection so its child `functional_key` mapper runs.
    add_mapping_annotation(model_method.DFT.xc, TEXT_KEY, '.@')


class XCFunctional(model_method.XCFunctional):
    # Set `functional_key` to the standard functional name from the FHI-aims XC
    # control string; the schema expands the name into components (family/kind)
    # and derives `jacobs_ladder`.
    add_mapping_annotation(
        model_method.XCFunctional.functional_key,
        TEXT_KEY,
        ('get_functional_key', ['.controlInOut_xc']),
    )


# class DFT(model_method.DFT):
#     model_method.DFT.xc_functionals.m_annotations.setdefault(
#         MAPPING_ANNOTATION_KEY, {}
#     ).update(dict(text=Mapper(mapper=('get_xc_functionals', ['.controlInOut_xc']))))


# class XCFunctional(model_method.XCFunctional):
#     model_method.XCFunctional.libxc_name.m_annotations.setdefault(
#         MAPPING_ANNOTATION_KEY, {}
#     ).update(dict(text=Mapper(mapper='.name')))


class GW(model_method.GW):
    add_mapping_annotation(
        model_method.GW.type, TEXT_GW_KEY, ('get_gw_flag', ['.gw_flag'])
    )


class ModelSystem(model_system.ModelSystem):
    add_mapping_annotation(model_system.Representation.m_def, TEXT_KEY, '.@')
    add_mapping_annotation(
        model_system.ModelSystem.positions,
        TEXT_KEY,
        '.structure.positions',
        unit='angstrom',
    )
    add_mapping_annotation(model_system.AtomsState.m_def, TEXT_KEY, '.structure.labels')


class Representation(model_system.Representation):
    add_mapping_annotation(
        model_system.Representation.lattice_vectors, TEXT_KEY, '.lattice_vectors'
    )
    add_mapping_annotation(
        model_system.Representation.periodic_boundary_conditions,
        TEXT_KEY,
        ('get_periodic_boundary_conditions', ['.@']),
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
    add_mapping_annotation(
        outputs.Outputs.electronic_band_gaps,
        TEXT_KEY,
        ('get_band_gaps', ['.@']),
    )
    add_mapping_annotation(
        outputs.Outputs.scf_steps,
        TEXT_KEY,
        ('get_scf_steps', ['.@']),
    )


class SCFSteps(outputs.SCFSteps):
    add_mapping_annotation(
        outputs.SCFSteps.delta_energies_total, TEXT_KEY, '.delta_energies_total'
    )
    add_mapping_annotation(
        outputs.SCFSteps.delta_charge_abs, TEXT_KEY, '.delta_charge_abs'
    )
    add_mapping_annotation(outputs.SCFSteps.durations, TEXT_KEY, '.durations')
    add_mapping_annotation(
        outputs.SCFSteps.code_specific_quantities, TEXT_KEY, '.code_specific_quantities'
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
        properties.ElectronicEigenvalues.value, TEXT_KEY, '.eigenvalues'
    )
    add_mapping_annotation(
        properties.ElectronicEigenvalues.occupation, TEXT_KEY, '.occupations'
    )
    add_mapping_annotation(
        properties.ElectronicEigenvalues.spin_channel, TEXT_KEY, '.spin_channel'
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
    add_mapping_annotation(variables.Energy2.m_def, TEXT_DOS_KEY, '.@')


class ElectronicBandGap(properties.ElectronicBandGap):
    add_mapping_annotation(properties.ElectronicBandGap.value, TEXT_KEY, '.value')
    add_mapping_annotation(
        properties.ElectronicBandGap.spin_channel, TEXT_KEY, '.spin_channel'
    )


class Energy2(variables.Energy2):
    add_mapping_annotation(
        variables.Energy2.points, TEXT_DOS_KEY, '.energies', unit='eV'
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
    # workflow.molecular_dynamics.MolecularDynamicsModel.m_def.m_annotations.setdefault(
    #     MAPPING_ANNOTATION_KEY, {}
    # ).update(dict(md_workflow=Mapper(mapper='.@')))
    workflow.molecular_dynamics.MolecularDynamicsResults.m_def.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(md_workflow=Mapper(mapper='.@')))


class GeometryOptimization(workflow.GeometryOptimization):
    workflow.geometry_optimization.GeometryOptimizationMethod.m_def.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(geo_opt_workflow=Mapper(mapper='.@')))
    workflow.geometry_optimization.GeometryOptimizationResults.m_def.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(geo_opt_workflow=Mapper(mapper='.@')))


class GeometryOptimizationMethod(
    workflow.geometry_optimization.GeometryOptimizationMethod
):
    workflow.geometry_optimization.GeometryOptimizationMethod.optimization_method.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(geo_opt_workflow=Mapper(mapper='.geometry_relaxation_method')))


class KSpace(numerical_settings.KSpace):
    add_mapping_annotation(numerical_settings.KSpace.k_mesh, TEXT_KEY, '.@')


class KMesh(numerical_settings.KMesh):
    add_mapping_annotation(numerical_settings.KMesh.grid, TEXT_KEY, '.k_grid')
    add_mapping_annotation(
        numerical_settings.KMesh.offset,
        TEXT_KEY,
        ('get_k_offset_with_default', ['.k_offset']),
    )


try:
    m_package.__init_metainfo__()
except Exception:
    pass
