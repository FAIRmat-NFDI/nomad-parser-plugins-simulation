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
    outputs,
    properties,
    variables,
)

m_package = SchemaPackage()


general.Simulation.m_def.m_annotations.setdefault(MAPPING_ANNOTATION_KEY, {}).update(
    dict(
        text=Mapper(mapper='@'),
        text_dos=Mapper(mapper='@'),
        text_gw=Mapper(mapper='@'),
        text_gw_workflow=Mapper(mapper='@'),
    )
)

# program
general.Simulation.program.m_annotations.setdefault(MAPPING_ANNOTATION_KEY, {}).update(
    dict(text=Mapper(mapper='.@'))
)
## program quantities
general.Program.version.m_annotations.setdefault(MAPPING_ANNOTATION_KEY, {}).update(
    dict(text=Mapper(mapper='.version'))
)

# dft method
model_method.DFT.m_def.m_annotations.setdefault(MAPPING_ANNOTATION_KEY, {}).update(
    dict(text=Mapper(mapper='.@'))
)
## xc functionals
model_method.DFT.xc_functionals.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(dict(text=Mapper(mapper=('get_xc_functionals', ['.controlInOut_xc']))))
### xc functionals quantities
model_method.XCFunctional.libxc_name.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(dict(text=Mapper(mapper='.name')))
# gw method
model_method.GW.m_def.m_annotations.setdefault(MAPPING_ANNOTATION_KEY, {}).update(
    dict(text_gw=Mapper(mapper='.@'))
)
## gw quantities
model_method.GW.type.m_annotations.setdefault(MAPPING_ANNOTATION_KEY, {}).update(
    dict(text_gw=Mapper(mapper=('get_gw_flag', ['.gw_flag'])))
)

# model_system
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
## cell
model_system.AtomicCell.m_def.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(dict(text=Mapper(mapper='.@')))
### cell quantities
model_system.AtomicCell.lattice_vectors.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(dict(text=Mapper(mapper='.lattice_vectors')))
model_system.AtomicCell.positions.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(
    dict(text=Mapper(mapper=('to_array', ['.structure']), unit='angstrom'))
)

# outputs
general.Simulation.outputs.m_annotations.setdefault(MAPPING_ANNOTATION_KEY, {}).update(
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
## variables
variables.Variables.n_points.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(dict(text=Mapper(mapper='.npoints')))
# TODO this does not work as points is scalar
# variables.Variables.points.m_annotations.setdefault(MAPPING_ANNOTATION_KEY, {}).update(
#     dict(text=Mapper(mapper='.points'))
# )
variables.Energy2.m_def.m_annotations.setdefault(MAPPING_ANNOTATION_KEY, {}).update(
    dict(text_dos=Mapper(mapper='.@'))
)
### energy variables
variables.Energy2.n_points.m_annotations.setdefault(MAPPING_ANNOTATION_KEY, {}).update(
    dict(text_dos=Mapper(mapper='.nenergies'))
)
variables.Energy2.points.m_annotations.setdefault(MAPPING_ANNOTATION_KEY, {}).update(
    dict(text_dos=Mapper(mapper='.energies'))
)
## total energy
outputs.Outputs.total_energies.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(dict(text=Mapper(mapper=('get_energies', ['.@']))))
### total energy quantities
properties.energies.TotalEnergy.value.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(dict(text=Mapper(mapper='.value')))
properties.energies.TotalEnergy.contributions.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(dict(text=Mapper(mapper='.components')))
#### total energy contributions
properties.energies.EnergyContribution.name.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(dict(text=Mapper(mapper='.name')))
## total forces
outputs.Outputs.total_forces.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(dict(text=Mapper(mapper=('get_forces', ['.@']))))
### total forces quantities
properties.forces.TotalForce.rank.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(dict(text=Mapper(mapper='.rank')))
properties.forces.TotalForce.value.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(dict(text=Mapper(mapper='.forces')))
## electronic eigenvalues
outputs.Outputs.electronic_eigenvalues.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(
    dict(
        text=Mapper(
            mapper=('get_eigenvalues', ['.eigenvalues', 'array_size_parameters'])
        )
    )
)
### electronic eigenvalues quantites
properties.ElectronicEigenvalues.n_bands.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(dict(text=Mapper(mapper='.nbands')))
properties.ElectronicEigenvalues.variables.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(dict(text=Mapper(mapper='.@', search='.@[0]')))
properties.ElectronicEigenvalues.value.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(dict(text=Mapper(mapper='.eigenvalues')))
properties.ElectronicEigenvalues.occupation.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(dict(text=Mapper(mapper='.occupations')))
## electronic dos
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
### dos quantities
properties.spectral_profile.DOSProfile.value.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(dict(text_dos=Mapper(mapper='.values')))
properties.spectral_profile.ElectronicDensityOfStates.value.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(dict(text_dos=Mapper(mapper='.values')))
### projected dos
properties.spectral_profile.ElectronicDensityOfStates.projected_dos.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(dict(text_dos=Mapper(mapper='.projected_dos')))


m_package.__init_metainfo__()
