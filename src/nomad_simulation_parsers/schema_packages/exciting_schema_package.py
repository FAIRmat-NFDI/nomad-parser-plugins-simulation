from nomad.datamodel.metainfo.annotations import Mapper
from nomad.metainfo import SchemaPackage
from nomad.parsing.file_parser.mapping_parser import MAPPING_ANNOTATION_KEY
from nomad_simulations.schema_packages import (
    general,
    model_method,
    outputs,
    properties,
    variables,
)

m_package = SchemaPackage()

general.Program.version.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    info=Mapper(mapper='.program_version')
)

general.Simulation.program.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    info=Mapper(mapper='.@')
)

model_method.XCFunctional.libxc_name.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    info=Mapper(mapper='.libxc'), input_xml=Mapper(mapper='.libxc')
)

model_method.DFT.xc_functionals.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    info=Mapper(mapper=('get_xc_functionals', ['.type'])),
    input_xml=Mapper(mapper=('get_xc_functionals', ['.libxc'])),
)


model_method.DFT.m_def.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    info=Mapper(mapper='.initialization.xc_functional'),
    input_xml=Mapper(mapper='.input.groundstate'),
)


properties.TotalEnergy.value.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    info=Mapper(mapper='.final.energy_total || energy_total')
)


variables.Variables.n_points.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    info=Mapper(mapper='.n_points'), eigval=Mapper(mapper='n_k_points')
)


properties.forces.TotalForce.variables.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    info=Mapper(mapper='.@')
)

properties.forces.TotalForce.rank.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    info=Mapper(mapper='.rank')
)
properties.forces.TotalForce.value.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    info=Mapper(mapper='.forces')
)


outputs.ElectronicEigenvalues.n_bands.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    eigval=Mapper(mapper='n_states')
)
outputs.ElectronicEigenvalues.variables.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(dict(eigval=Mapper(mapper='@')))
outputs.ElectronicEigenvalues.value.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    eigval=Mapper(mapper='.eigenvalues'),
)
outputs.ElectronicEigenvalues.occupation.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    eigval=Mapper(mapper='.occupancies')
)


outputs.Outputs.total_energies.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    info=Mapper(mapper='.@')
)
outputs.Outputs.total_forces.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    info=Mapper(mapper=('get_forces', ['.@']))
)
outputs.Outputs.electronic_eigenvalues.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    eigval=Mapper(mapper=('get_eigenvalues', ['.@']))
)


general.Simulation.outputs.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    info=Mapper(mapper=('get_configurations', ['.@'])), eigval=Mapper(mapper='.@')
)


general.Simulation.m_def.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    info=Mapper(mapper='@'), input_xml=Mapper(mapper='@'), eigval=Mapper(mapper='@')
)


m_package.__init_metainfo__()
