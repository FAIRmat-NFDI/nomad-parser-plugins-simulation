from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass

from nomad.datamodel import EntryArchive
from nomad.datamodel.metainfo.annotations import Mapper
from nomad.metainfo import SchemaPackage
from nomad.parsing.file_parser.mapping_parser import MAPPING_ANNOTATION_KEY
from nomad_simulations.schema_packages import (
    general,
    model_method,
    model_system,
    outputs,
    properties,
    workflow,
)

m_package = SchemaPackage()

EntryArchive.m_def.m_annotations.setdefault(MAPPING_ANNOTATION_KEY, {}).update(
    dict(
        text=Mapper(mapper='@'),
        text_dos=Mapper(mapper='@'),
        text_gw=Mapper(mapper='@'),
        single_point=Mapper(mapper='@'),
        geo_opt_workflow=Mapper(mapper='@'),
        md_workflow=Mapper(mapper='@'),
    )
)

general.Simulation.m_def.m_annotations.setdefault(MAPPING_ANNOTATION_KEY, {}).update(
    dict(
        text=Mapper(mapper='@'),
        text_dos=Mapper(mapper='@'),
        text_gw=Mapper(mapper='@'),
    )
)


class Simulation(general.Simulation):
    general.Simulation.program.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(text=Mapper(mapper='.@')))
    general.Simulation.model_system.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(
        dict(
            text=Mapper(
                mapper=(
                    'get_sections',
                    ['.@'],
                    dict(include=['lattice_vectors', 'structure', 'sub_structure']),
                )
            )
        )
    )
    # DFT method
    model_method.DFT.m_def.m_annotations.setdefault(MAPPING_ANNOTATION_KEY, {}).update(
        dict(text=Mapper(mapper='.@'))
    )
    # gw method
    model_method.GW.m_def.m_annotations.setdefault(MAPPING_ANNOTATION_KEY, {}).update(
        dict(text_gw=Mapper(mapper='.@'))
    )
    # electronic structure outputs
    outputs.Outputs.m_def.m_annotations.setdefault(MAPPING_ANNOTATION_KEY, {}).update(
        dict(
            text=Mapper(
                mapper=(
                    'get_sections',
                    ['.@'],
                    dict(
                        include=[
                            'energy',
                            'energy_components',
                            'forces',
                            'eigenvalues',
                        ]
                    ),
                )
            ),
            text_dos=Mapper(
                mapper=(
                    'get_sections',
                    ['.@'],
                    dict(
                        include=[
                            'total_dos_files',
                            'species_projected_dos_files',
                        ]
                    ),
                )
            ),
        )
    )


class Program(general.Program):
    general.Program.version.m_annotations.setdefault(MAPPING_ANNOTATION_KEY, {}).update(
        dict(text=Mapper(mapper='.version'))
    )


class DFT(model_method.DFT):
    model_method.DFT.xc_functionals.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(text=Mapper(mapper=('get_xc_functionals', ['.controlInOut_xc']))))


class XCFunctional(model_method.XCFunctional):
    model_method.XCFunctional.libxc_name.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(text=Mapper(mapper='.name')))


class GW(model_method.GW):
    model_method.GW.type.m_annotations.setdefault(MAPPING_ANNOTATION_KEY, {}).update(
        dict(text_gw=Mapper(mapper=('get_gw_flag', ['.gw_flag'])))
    )


class ModelSystem(model_system.ModelSystem):
    model_system.AtomicCell.m_def.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(text=Mapper(mapper='.@')))
    model_system.ModelSystem.positions.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(text=Mapper(mapper='.structure.positions', unit='angstrom')))
    model_system.AtomsState.m_def.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(text=Mapper(mapper='.structure.labels')))


class AtomicCell(model_system.AtomicCell):
    model_system.AtomicCell.lattice_vectors.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(text=Mapper(mapper='.lattice_vectors')))


class AtomsState(model_system.AtomsState):
    model_system.AtomsState.chemical_symbol.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(text=Mapper(mapper='.@')))


class Outputs(outputs.Outputs):
    outputs.Outputs.total_energies.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(text=Mapper(mapper=('get_energies', ['.@']))))
    outputs.Outputs.total_forces.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(text=Mapper(mapper=('get_forces', ['.@']))))
    outputs.Outputs.electronic_eigenvalues.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(
        dict(
            text=Mapper(
                mapper=('get_eigenvalues', ['.eigenvalues', 'array_size_parameters'])
            )
        )
    )
    outputs.Outputs.electronic_dos.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(
        dict(
            text_dos=Mapper(
                mapper=(
                    'get_dos',
                    [
                        '.total_dos_files',
                        '.atom_projected_dos_files',
                        '.species_projected_dos_files',
                    ],
                )
            )
        )
    )


class TotalEnergy(properties.energies.TotalEnergy):
    properties.energies.TotalEnergy.value.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(text=Mapper(mapper='.value')))
    properties.energies.TotalEnergy.contributions.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(text=Mapper(mapper='.components')))


class BaseEnergy(properties.energies.BaseEnergy):
    properties.energies.BaseEnergy.name.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(text=Mapper(mapper='.name')))


class TotalForce(properties.forces.TotalForce):
    properties.forces.TotalForce.value.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(text=Mapper(mapper='.forces')))


class ElectronicEigenvalues(properties.ElectronicEigenvalues):
    properties.ElectronicEigenvalues.n_bands.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(text=Mapper(mapper='.nbands')))
    properties.ElectronicEigenvalues.value.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(text=Mapper(mapper='.eigenvalues')))
    properties.ElectronicEigenvalues.occupation.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(text=Mapper(mapper='.occupations')))


"""class DOSProfile(properties.spectral_profile.DOSProfile):
    ### dos quantities
    properties.spectral_profile.DOSProfile.value.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(text_dos=Mapper(mapper='.values')))"""


"""class ElectronicDensityOfStates(properties.spectral_profile.ElectronicDensityOfStates):
    properties.spectral_profile.ElectronicDensityOfStates.value.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(text_dos=Mapper(mapper='.values')))
    ### projected dos
    properties.spectral_profile.ElectronicDensityOfStates.projected_dos.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(text_dos=Mapper(mapper='.projected_dos')))"""


class SimulationWorkflow(workflow.general.SimulationWorkflow):
    # TODO find a more elegant fix to not parse tasks recursively, this will be filled
    # in by workflow normalizer from outputs
    workflow.general.SimulationWorkflow.tasks.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(
        dict(
            geo_opt_workflow=Mapper(mapper='.tasks'),
            md_workflow=Mapper(mapper='.tasks'),
        )
    )


# workflow
workflow.single_point.SinglePoint.m_def.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(dict(single_point=Mapper(mapper='.@')))
# geometry optimization workflow
workflow.geometry_optimization.GeometryOptimization.m_def.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(dict(geo_opt_workflow=Mapper(mapper='.@')))
# molecular dynamics workflow
workflow.molecular_dynamics.MolecularDynamics.m_def.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(dict(md_workflow=Mapper(mapper='.@')))


class MolecularDynamics(workflow.MolecularDynamics):
    workflow.molecular_dynamics.MolecularDynamicsModel.m_def.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(md_workflow=Mapper(mapper='.@')))
    workflow.molecular_dynamics.MolecularDynamicsResults.m_def.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(md_workflow=Mapper(mapper='.@')))


class GeometryOptimization(workflow.GeometryOptimization):
    workflow.geometry_optimization.GeometryOptimizationModel.m_def.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(geo_opt_workflow=Mapper(mapper='.@')))
    workflow.geometry_optimization.GeometryOptimizationResults.m_def.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(geo_opt_workflow=Mapper(mapper='.@')))


class GeometryOptimizationModel(
    workflow.geometry_optimization.GeometryOptimizationModel
):
    workflow.geometry_optimization.GeometryOptimizationModel.optimization_method.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(geo_opt_workflow=Mapper(mapper='.geometry_relaxation_method')))


class MolecularDynamicsModel(workflow.molecular_dynamics.MolecularDynamicsModel):
    workflow.molecular_dynamics.MolecularDynamicsModel.integration_timestep.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(md_workflow=Mapper(mapper='.control_inout.md_timestep')))
    workflow.molecular_dynamics.MolecularDynamicsModel.thermodynamic_ensemble.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(md_workflow=Mapper(mapper='.control_inout.md_run[0]')))


class MolecularDynamicsResults(workflow.molecular_dynamics.MolecularDynamicsResults):
    workflow.molecular_dynamics.MolecularDynamicsResults.temperature.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(
        dict(
            md_workflow=Mapper(
                mapper='molecular_dynamics[*].md_calculation_info.'
                '"Temperature (nuclei)"'
            )
        )
    )


try:
    m_package.__init_metainfo__()
except Exception:
    pass
