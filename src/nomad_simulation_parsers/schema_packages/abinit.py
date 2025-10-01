from nomad.datamodel.metainfo.annotations import Mapper
from nomad.metainfo import SchemaPackage
from nomad.parsing.file_parser.mapping_parser import MAPPING_ANNOTATION_KEY
from nomad_simulations.schema_packages import (
    general,
    model_method,
    model_system,
    outputs,
    workflow,
)

m_package = SchemaPackage()


class GeometryOptimizationMethod(
    workflow.geometry_optimization.GeometryOptimizationMethod
):
    workflow.geometry_optimization.GeometryOptimizationMethod.optimization_method.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(out=Mapper(mapper=('get_workflow_method', []))))
    workflow.geometry_optimization.GeometryOptimizationMethod.convergence_tolerance_energy_difference.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(
        dict(
            out=Mapper(
                mapper=(
                    'get_input_var',
                    [],
                    dict(name='tolmxde', n_dataset=1, default=0.0),
                ),
                unit='hartree',
            )
        )
    )
    workflow.geometry_optimization.GeometryOptimizationMethod.convergence_tolerance_force_maximum.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(
        dict(
            out=Mapper(
                mapper=(
                    'get_input_var',
                    [],
                    dict(name='tolmxf', n_dataset=1, default=0.0),
                ),
                unit='hartree/bohr',
            )
        )
    )


workflow.geometry_optimization.GeometryOptimizationMethod.m_def.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(dict(out=Mapper(mapper='@')))

workflow.GeometryOptimization.m_def.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(dict(out=Mapper(mapper='@')))


class Program(general.Program):
    general.Program.version.m_annotations.setdefault(MAPPING_ANNOTATION_KEY, {}).update(
        dict(out=Mapper(mapper='.program_version'))
    )


class AtomsState(model_system.AtomsState):
    model_system.AtomsState.chemical_symbol.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(out=Mapper(mapper='.label')))


class AtomicCell(model_system.AtomicCell):
    model_system.AtomicCell.lattice_vectors.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(out=Mapper(mapper='dataset[0].x_abinit_vprim', unit='bohr')))


class ModelSystem(model_system.ModelSystem):
    model_system.AtomicCell.m_def.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(out=Mapper(mapper='.@')))
    model_system.AtomsState.m_def.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(out=Mapper(mapper=('get_atoms', []), cache=True)))
    model_system.ModelSystem.positions.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(out=Mapper(mapper='.cartesian_coordinates', unit='bohr')))


class XCFunctional(model_method.XCFunctional):
    model_method.XCFunctional.libxc_name.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(out=Mapper(mapper='.XC_functional_name')))


class DFT(model_method.DFT):
    model_method.DFT.xc_functionals.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(out=Mapper(mapper=('get_xc_functionals', []))))


class TotalEnergy(outputs.TotalEnergy):
    outputs.TotalEnergy.value.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(out=Mapper(mapper='.value || .energy_total')))
    outputs.TotalEnergy.name.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(out=Mapper(mapper='.name')))
    outputs.TotalEnergy.contributions.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(out=Mapper(mapper=('get_energy_contributions', ['.@']))))


class TotalForce(outputs.TotalForce):
    outputs.TotalForce.value.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(out=Mapper(mapper='.cartesian_forces')))


class ElectronicDensityOfStates(outputs.ElectronicDensityOfStates):
    outputs.ElectronicDensityOfStates.value.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(dos=Mapper(mapper='.value', unit='1 / hartree')))


class ElectronicBandStructure(outputs.ElectronicBandStructure):
    outputs.ElectronicBandStructure.value.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(out=Mapper(mapper='.energies', unit='hartree')))


class Outputs(outputs.Outputs):
    outputs.Outputs.total_energies.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(out=Mapper(mapper='.@')))
    outputs.Outputs.total_forces.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(out=Mapper(mapper='.@')))
    outputs.Outputs.electronic_dos.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(dos=Mapper(mapper=('get_dos', ['.data']))))
    outputs.Outputs.electronic_band_structures.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(
        dict(
            out=Mapper(
                mapper=('get_bandstructures', ['.eigenvalues', '.occupation_numbers'])
            )
        )
    )


class Simulation(general.Simulation):
    general.Simulation.program.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(out=Mapper(mapper='.@')))
    general.Simulation.datetime.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(
        dict(
            out=Mapper(
                mapper=('get_datetime', ['x_abinit_start_date', 'x_abinit_start_time'])
            )
        )
    )
    general.Simulation.model_system.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(out=Mapper(mapper=('get_systems', []))))
    model_method.DFT.m_def.m_annotations.setdefault(MAPPING_ANNOTATION_KEY, {}).update(
        dict(out=Mapper(mapper='.@'))
    )
    general.Simulation.outputs.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(out=Mapper(mapper=('get_outputs', [])), dos=Mapper(mapper='.@')))


general.Simulation.m_def.m_annotations.setdefault(MAPPING_ANNOTATION_KEY, {}).update(
    dict(out=Mapper(mapper='@'), dos=Mapper(mapper='.@'))
)


try:
    m_package.__init_metainfo__()
except Exception:
    pass
