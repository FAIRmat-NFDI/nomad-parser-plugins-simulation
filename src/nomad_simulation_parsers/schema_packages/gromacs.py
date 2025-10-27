from nomad.datamodel.metainfo.annotations import Mapper
from nomad.metainfo import SchemaPackage
from nomad.parsing.file_parser.mapping_parser import MAPPING_ANNOTATION_KEY
from nomad_simulations.schema_packages import (
    general,
    model_system,
    outputs,
)

m_package = SchemaPackage()


class Program(general.Program):
    general.Program.version.m_annotations.setdefault(MAPPING_ANNOTATION_KEY, {}).update(
        dict(log=Mapper(mapper=('get_version', ['.version'])))
    )


class AtomsState(model_system.AtomsState):
    model_system.AtomsState.label.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(tpr=Mapper(mapper='.@')))


class AtomicCell(model_system.AtomicCell):
    model_system.AtomicCell.lattice_vectors.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(tpr=Mapper(mapper='.lattice_vectors')))
    model_system.AtomicCell.periodic_boundary_conditions.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(log=Mapper(mapper='.pbc')))


class ModelSystem(model_system.ModelSystem):
    model_system.ModelSystem.positions.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(tpr=Mapper(mapper='.positions')))
    model_system.ModelSystem.velocities.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(tpr=Mapper(mapper='.velocities')))
    model_system.AtomsState.m_def.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(tpr=Mapper(mapper='.labels')))
    model_system.AtomicCell.m_def.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(log=Mapper(mapper='.@'), tpr=Mapper(mapper='.@')))


class TotalEnergy(outputs.TotalEnergy):
    outputs.TotalEnergy.name.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(log=Mapper(mapper='.label'), edr=Mapper(mapper='.label')))
    outputs.TotalEnergy.value.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(log=Mapper(mapper='.value'), edr=Mapper(mapper='.value')))
    outputs.TotalEnergy.contributions.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(
        dict(log=Mapper(mapper='.contributions'), edr=Mapper(mapper='.contributions'))
    )


class TotalForce(outputs.TotalForce):
    outputs.TotalForce.value.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(tpr=Mapper(mapper='.@')))


class Outpus(outputs.Outputs):
    outputs.Outputs.total_energies.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(log=Mapper(mapper='.energy'), edr=Mapper(mapper='.energy')))
    outputs.Outputs.model_system_ref.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(log=Mapper(mapper='.system_ref'), edr=Mapper(mapper='.system_ref')))
    outputs.Outputs.total_forces.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(tpr=Mapper(mapper='.forces')))


class Simulation(general.Simulation):
    general.Simulation.program.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(log=Mapper(mapper='.header')))
    general.Simulation.model_system.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(
        dict(
            log=Mapper(mapper=('get_configurations', [])),
            tpr=Mapper(mapper=('get_configurations', [])),
        )
    )
    general.Simulation.outputs.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(
        dict(
            log=Mapper(mapper=('get_outputs', [])),
            tpr=Mapper(mapper=('get_outputs', [])),
            edr=Mapper(mapper=('get_outputs', ['.@'])),
        )
    )


general.Simulation.m_def.m_annotations.setdefault(MAPPING_ANNOTATION_KEY, {}).update(
    dict(log=Mapper(mapper='@'), tpr=Mapper(mapper='@'), edr=Mapper(mapper='@'))
)

try:
    m_package.__init_metainfo__()
except Exception:
    pass
