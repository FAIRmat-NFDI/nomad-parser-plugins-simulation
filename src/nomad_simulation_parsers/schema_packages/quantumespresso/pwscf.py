from nomad.datamodel.metainfo.annotations import Mapper
from nomad.metainfo import SchemaPackage
from nomad.parsing.file_parser.mapping_parser import MAPPING_ANNOTATION_KEY
from nomad_simulations.schema_packages import general, outputs

m_package = SchemaPackage()


class TotalForce(outputs.TotalForce):
    outputs.TotalForce.value.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(out=Mapper(mapper='.value || .forces', unit='rydberg/bohr')))
    outputs.TotalForce.contributions.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(
        dict(
            out=Mapper(mapper=('get_force_contributions', ['.@']), unit='rydberg/bohr')
        )
    )


class ElectronicEigenvalues(outputs.ElectronicEigenvalues):
    outputs.ElectronicEigenvalues.value.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(out=Mapper(mapper='.eigenvalues')))


class Outputs(outputs.Outputs):
    outputs.Outputs.total_forces.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(out=Mapper(mapper='.@')))
    outputs.Outputs.electronic_eigenvalues.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(out=Mapper(mapper=('get_eigenvalues', ['.@']))))


class Simulation(general.Simulation):
    general.Simulation.model_system.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(out=Mapper(mapper=('get_configurations', ['.@']), cache=True)))
    general.Simulation.outputs.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(out=Mapper(mapper=('get_configurations', ['.@']), cache=True)))


try:
    m_package.__init_metainfo__()
except Exception:
    pass
