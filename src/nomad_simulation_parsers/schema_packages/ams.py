from nomad.metainfo import SchemaPackage
from nomad_simulations.schema_packages import general
from nomad.datamodel.metainfo.annotations import Mapper

from nomad.parsing.file_parser.mapping_parser import MAPPING_ANNOTATION_KEY

m_package = SchemaPackage()


class Program(general.Program):
    general.Program.version.m_annotations.setdefault(MAPPING_ANNOTATION_KEY, {}).update(
        dict(out=Mapper(mapper='.program_version'))
    )


class Simulation(general.Simulation):
    general.Simulation.program.m_annotations.setdefault(MAPPING_ANNOTATION_KEY, {}).update(
        dict(out=Mapper(mapper='.@'))
    )


Simulation.m_def.m_annotations.setdefault(MAPPING_ANNOTATION_KEY, {}).update(
    dict(
        out=Mapper(mapper='@')
    )
)


try:
    m_package.__init_metainfo__()
except Exception:
    pass
