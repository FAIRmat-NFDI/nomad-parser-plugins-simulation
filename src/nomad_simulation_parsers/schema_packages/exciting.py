from nomad.datamodel.metainfo.annotations import Mapper
from nomad.metainfo import SchemaPackage
from nomad.parsing.file_parser.mapping_parser import MAPPING_ANNOTATION_KEY
from nomad_simulations.schema_packages import (
    atoms_state,
    general,
    model_method,
    model_system,
    numerical_settings,
    outputs,
    properties,
    variables,
)

m_package = SchemaPackage()


# simulation
general.Simulation.m_def.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    info=Mapper(mapper='@'),
    input_xml=Mapper(mapper='@'),
    eigval=Mapper(mapper='@'),
    bandstructure_xml=Mapper(mapper='@'),
    dos_xml=Mapper(mapper='@'),
)
## program
general.Simulation.program.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    info=Mapper(mapper='.@')
)
### program quantities
general.Program.version.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    info=Mapper(mapper='.program_version')
)
## model_method
model_method.DFT.m_def.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    info=Mapper(mapper='.initialization.xc_functional'),
    input_xml=Mapper(mapper='.input.groundstate'),
    bandstructure_xml=Mapper(mapper='.@'),
)
### numerical_settings
numerical_settings.KSpace.m_def.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    bandstructure_xml=Mapper(mapper='.@')
)
#### numerical_settings sub sections
numerical_settings.KSpace.k_line_path.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    bandstructure_xml=Mapper(mapper='.@')
)
##### k_line_path
numerical_settings.KLinePath.high_symmetry_path_names.m_annotations[
    MAPPING_ANNOTATION_KEY
] = dict(bandstructure_xml=Mapper(mapper=r'bandstructure.vertex[*]."@label"'))
numerical_settings.KLinePath.high_symmetry_path_values.m_annotations[
    MAPPING_ANNOTATION_KEY
] = dict(
    bandstructure_xml=Mapper(
        mapper=('reshape_coords', [r'bandstructure.vertex[*]."@coord"'])
    )
)
### xc_functionals
model_method.DFT.xc_functionals.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    info=Mapper(mapper=('get_xc_functionals', ['.type'])),
    input_xml=Mapper(mapper=('get_xc_functionals', ['.libxc'])),
)
#### xc_functional quantities
model_method.XCFunctional.libxc_name.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    info=Mapper(mapper='.libxc'), input_xml=Mapper(mapper='.libxc')
)
## model_system
general.Simulation.model_system.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    info=Mapper(mapper=('get_configurations', ['.@']), cache=True)
)
### cell
model_system.AtomicCell.m_def.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    info=Mapper(mapper=('get_atoms', ['.atomic_positions']))
)
#### cell quantities
model_system.AtomicCell.positions.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    info=Mapper(mapper='.positions')
)
model_system.AtomicCell.atoms_state.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    info=Mapper(mapper='.atoms')
)
##### atoms_state quantities
atoms_state.AtomsState.chemical_symbol.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    info=Mapper(mapper='.symbol')
)
## outputs
general.Simulation.outputs.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    info=Mapper(mapper=('get_configurations', ['.@'])),
    eigval=Mapper(mapper='.@'),
    bandstructure_xml=Mapper(mapper='.@'),
    dos_xml=Mapper(mapper='.@'),
)
### variables
variables.Variables.n_points.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    info=Mapper(mapper='.n_points'),
    eigval=Mapper(mapper='n_k_points'),
    bandstructure_xml=Mapper(mapper='.n_kpoints'),
)
### total_energies
outputs.Outputs.total_energies.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    info=Mapper(mapper='.@')
)
#### total_energies quantities
properties.TotalEnergy.value.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    info=Mapper(mapper='.final.energy_total || energy_total')
)
### total_forces
outputs.Outputs.total_forces.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    info=Mapper(mapper=('get_forces', ['.@']))
)
#### total_forces quantities
properties.forces.TotalForce.variables.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    info=Mapper(mapper='.@')
)
properties.forces.TotalForce.rank.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    info=Mapper(mapper='.rank')
)
properties.forces.TotalForce.value.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    info=Mapper(mapper='.forces')
)
### electronic_eigenvalues
outputs.Outputs.electronic_eigenvalues.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    eigval=Mapper(mapper=('get_eigenvalues', ['.@']))
)
#### electronic_eigenvalues quantities
outputs.ElectronicEigenvalues.n_bands.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    eigval=Mapper(mapper='.n_states')
)
##### we need to use setdefault when annot was defined previously
outputs.ElectronicEigenvalues.variables.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(dict(eigval=Mapper(mapper='@')))
outputs.ElectronicEigenvalues.value.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    eigval=Mapper(mapper='.eigenvalues'),
)
outputs.ElectronicEigenvalues.occupation.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    eigval=Mapper(mapper='.occupancies')
)
### electronic_band_structures
outputs.Outputs.electronic_band_structures.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    bandstructure_xml=Mapper(mapper=('get_bandstructures', ['.@']))
)
#### electronic_band_structures quantities
outputs.ElectronicBandStructure.n_bands.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(dict(bandstructure_xml=Mapper(mapper='.n_states')))
outputs.ElectronicBandStructure.variables.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(dict(bandstructure_xml=Mapper(mapper='.@')))
outputs.ElectronicBandStructure.value.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(dict(bandstructure_xml=Mapper(mapper='.energies')))
### electronic_dos
outputs.Outputs.electronic_dos.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    dos_xml=Mapper(mapper=('dos.totaldos.diagram'))
)
#### electronic_dos quantities
variables.Energy2.m_def.m_annotations.setdefault(MAPPING_ANNOTATION_KEY, {}).update(
    dict(dos_xml=Mapper(mapper='.@'))
)
variables.Energy2.n_points.m_annotations.setdefault(MAPPING_ANNOTATION_KEY, {}).update(
    dict(dos_xml=Mapper(mapper='length(.point)'))
)
variables.Energy2.points.m_annotations.setdefault(MAPPING_ANNOTATION_KEY, {}).update(
    dict(dos_xml=Mapper(mapper=('to_float', [r'.point[*]."@e"']), unit='hartree'))
)
###### TODO read unit from axis
outputs.ElectronicDensityOfStates.value.m_annotations[MAPPING_ANNOTATION_KEY] = dict(
    dos_xml=Mapper(mapper=('to_float', [r'.point[*]."@dos"']), unit='1/hartree')
)
outputs.ElectronicDensityOfStates.projected_dos.m_annotations[
    MAPPING_ANNOTATION_KEY
] = dict(dos_xml=Mapper(mapper='dos.partialdos.diagram'))
m_package.__init_metainfo__()
