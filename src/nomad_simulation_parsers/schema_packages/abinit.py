from nomad.datamodel.metainfo.annotations import Mapper
from nomad.metainfo import SchemaPackage
from nomad.parsing.file_parser.mapping_parser import MAPPING_ANNOTATION_KEY
from nomad_simulations.schema_packages import general, workflow

m_package = SchemaPackage()


class GeometryOptimizationModel(
    workflow.geometry_optimization.GeometryOptimizationModel
):
    workflow.geometry_optimization.GeometryOptimizationModel.optimization_method.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(out=Mapper(mapper=('get_workflow_method', []))))
    workflow.geometry_optimization.GeometryOptimizationModel.convergence_tolerance_energy_difference.m_annotations.setdefault(
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
    workflow.geometry_optimization.GeometryOptimizationModel.convergence_tolerance_force_maximum.m_annotations.setdefault(
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


workflow.geometry_optimization.GeometryOptimizationModel.m_def.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(dict(out=Mapper(mapper='@')))

workflow.GeometryOptimization.m_def.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(dict(out=Mapper(mapper='@')))


class Program(general.Program):
    general.Program.version.m_annotations.setdefault(MAPPING_ANNOTATION_KEY, {}).update(
        dict(out=Mapper(mapper='.program_version'))
    )


class Simulation(general.Simulation):
    general.Simulation.program.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(out=Mapper(mapper='.@')))


Simulation.m_def.m_annotations.setdefault(MAPPING_ANNOTATION_KEY, {}).update(
    dict(out=Mapper(mapper='@'))
)


try:
    m_package.__init_metaino__()
except Exception:
    pass
