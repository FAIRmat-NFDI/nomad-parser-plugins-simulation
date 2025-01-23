from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass

from nomad.datamodel.metainfo.annotations import Mapper as MapperAnnotation
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
)

m_package = SchemaPackage()


general.Simulation.m_def.m_annotations.setdefault(MAPPING_ANNOTATION_KEY, {}).update(
    dict(
        xml=MapperAnnotation(mapper='modeling'),
        xml2=MapperAnnotation(mapper='modeling'),
        outcar=MapperAnnotation(mapper='@'),
    )
)

# program
general.Simulation.program.m_annotations.setdefault(MAPPING_ANNOTATION_KEY, {}).update(
    dict(
        xml=MapperAnnotation(mapper='.generator'),
        outcar=MapperAnnotation(mapper='.header'),
    )
)
## program quantities
general.Program.name.m_annotations.setdefault(MAPPING_ANNOTATION_KEY, {}).update(
    dict(
        xml=MapperAnnotation(
            mapper='.i[?"@name"==\'program\'] | [0].__value',
        )
    )
)
general.Program.version.m_annotations.setdefault(MAPPING_ANNOTATION_KEY, {}).update(
    dict(
        xml=MapperAnnotation(
            mapper='.i[?"@name"==\'version\'] | [0].__value',
        ),
        outcar=MapperAnnotation(
            mapper=('get_version', ['.@']),
        ),
    )
)
general.Program.compilation_host.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(dict(xml=MapperAnnotation(mapper='.i[?"@name"==\'platform\'] | [0].__value')))
# dft method
model_method.DFT.m_def.m_annotations.setdefault(MAPPING_ANNOTATION_KEY, {}).update(
    dict(
        xml=MapperAnnotation(mapper='.parameters.separator[?"@name"==\'electronic\']'),
        outcar=MapperAnnotation(mapper='parameters'),
    )
)
## xc functional
model_method.DFT.xc_functionals.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(
    dict(
        xml=MapperAnnotation(
            mapper='.separator[?"@name"==\'electronic exchange-correlation\']'
        ),
        outcar=MapperAnnotation(mapper=('get_xc_functionals', ['.@'])),
    )
)
### xc functional quantities
model_method.XCFunctional.libxc_name.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(
    dict(
        xml=MapperAnnotation(
            # TODO add LDA & mGGA, convert_xc
            mapper='.i[?"@name"==\'GGA\'] | [0].__value'
        ),
        outcar=MapperAnnotation(mapper='.name'),
    )
)
model_method.DFT.exact_exchange_mixing_factor.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(
    dict(
        xml=MapperAnnotation(
            mapper=(
                'mix_alpha',
                [
                    '.i[?"@name"==\'HFALPHA\'] | [0].__value',
                    '.i[?"@name"==\'LHFCALC\'] | [0].__value',
                ],
            )
        )
    )
)  # TODO convert vasp bool
## numerical settings
numerical_settings.KSpace.m_def.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(dict(xml=MapperAnnotation(mapper='modeling.kpoints')))
### k-mesh
numerical_settings.KSpace.k_mesh.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(dict(xml=MapperAnnotation(mapper='.@')))
#### k-mesh quantities
numerical_settings.KMesh.grid.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(
    dict(
        xml=MapperAnnotation(
            mapper='.generation.v[?"@name"==\'divisions\'] | [0].__value'
        )
    )
)
numerical_settings.KMesh.offset.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(
    dict(
        xml=MapperAnnotation(mapper='.generation.v[?"@name"==\'shift\'] | [0].__value')
    )
)
numerical_settings.KMesh.points.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(dict(xml=MapperAnnotation(mapper='.varray[?"@name"==\'kpointlist\'].v | [0]')))

numerical_settings.KMesh.weights.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(dict(xml=MapperAnnotation(mapper='.varray[?"@name"==\'weights\'].v | [0]')))
# model system
general.Simulation.model_system.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(
    dict(
        xml=MapperAnnotation(mapper='.calculation'),
        outcar=MapperAnnotation(mapper='.calculation'),
    )
)
## cell
model_system.AtomicCell.m_def.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(
    dict(
        xml=MapperAnnotation(mapper='.structure'),
        outcar=MapperAnnotation(mapper='.@'),
    )
)
### cell quantities
model_system.AtomicCell.positions.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(
    dict(
        xml=MapperAnnotation(mapper='.varray.v', unit='angstrom'),
        outcar=MapperAnnotation(
            mapper='.positions_forces', unit='angstrom', search='@ | [0]'
        ),
    )
)
model_system.AtomicCell.lattice_vectors.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(
    dict(
        xml=MapperAnnotation(
            mapper='.crystal.varray[?"@name"==\'basis\'] | [0].v', unit='angstrom'
        ),
        outcar=MapperAnnotation(
            mapper='.lattice_vectors', unit='angstrom', search='@ | [0]'
        ),
    )
)
# outputs
general.Simulation.outputs.m_annotations.setdefault(MAPPING_ANNOTATION_KEY, {}).update(
    dict(
        xml=MapperAnnotation(mapper='.calculation'),
        xml2=MapperAnnotation(mapper='.calculation'),
        outcar=MapperAnnotation(mapper='.calculation'),
    )
)
## variables
variables.Variables.n_points.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(
    dict(
        xml=MapperAnnotation(mapper='.npoints'),
        outcar=MapperAnnotation(mapper='.npoints'),
    )
)

## total energies
outputs.Outputs.total_energies.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(
    dict(
        xml=MapperAnnotation(mapper='.energy'),
        outcar=MapperAnnotation(mapper='.energies'),
    )
)
### total energies quantities
# value is already defined in TotalEnergy since they use the same def
# get_energy function should be able to handle extraction from both sources
properties.energies.TotalEnergy.value.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(
    dict(
        xml=MapperAnnotation(
            mapper=(
                'get_data',
                ['.@'],
                dict(path='.i[?"@name"==\'e_fr_energy\'] | [0].__value'),
            ),
            unit='eV',
        ),
        outcar=MapperAnnotation(
            mapper=('get_data', ['.@'], dict(path='.energy_total')), unit='eV'
        ),
    )
)
properties.energies.TotalEnergy.contributions.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(
    dict(
        xml=MapperAnnotation(
            mapper=(
                'get_energy_contributions',
                ['.i'],
                dict(exclude=['e_fr_energy']),
            )
        ),
        outcar=MapperAnnotation(
            mapper=(
                'get_energy_contributions',
                ['.@'],
                dict(exclude=['energy_total']),
            )
        ),
    )
)
#### contributions
properties.energies.EnergyContribution.name.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(
    dict(
        xml=MapperAnnotation(mapper='."@name"'),
        outcar=MapperAnnotation(mapper='.name'),
    )
)
## total force
outputs.Outputs.total_forces.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(
    dict(
        xml=MapperAnnotation(mapper=('get_forces', ['.@'])),
        outcar=MapperAnnotation(mapper=('get_forces', ['.@'])),
    )
)
### total forces quantities
properties.forces.TotalForce.rank.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(
    dict(
        xml=MapperAnnotation(mapper='.rank'),
        outcar=MapperAnnotation(mapper='.rank'),
    )
)
properties.forces.TotalForce.variables.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(dict(xml=MapperAnnotation(mapper='.@'), outcar=MapperAnnotation(mapper='.@')))
properties.forces.TotalForce.value.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(
    dict(
        xml=MapperAnnotation(mapper='.forces', unit='eV/angstrom'),
        outcar=MapperAnnotation(mapper='.forces', unit='eV/angstrom'),
    )
)
## eigenvalues
outputs.Outputs.electronic_eigenvalues.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(
    dict(
        xml=MapperAnnotation(mapper=('get_eigenvalues', ['eigenvalues'])),
        xml2=MapperAnnotation(mapper=('get_eigenvalues', ['eigenvalues'])),
        outcar=MapperAnnotation(
            mapper=('get_eigenvalues', ['.eigenvalues', 'parameters'])
        ),
    )
)
### eigenvalues quantities
outputs.ElectronicEigenvalues.n_bands.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(
    dict(
        xml=MapperAnnotation(mapper='length(.array.set.set.set[0].r)'),
        xml2=MapperAnnotation(mapper='length(.array.set.set.set[0].r)'),
        outcar=MapperAnnotation(mapper='.n_bands'),
    )
)
outputs.ElectronicEigenvalues.variables.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(dict(outcar=MapperAnnotation(mapper='.@')))

# TODO This only works for non-spin pol
outputs.ElectronicEigenvalues.occupation.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(
    dict(
        outcar=MapperAnnotation(mapper='.occupations'),
        xml2=MapperAnnotation(mapper='.occupations'),
    )
)
outputs.ElectronicEigenvalues.value.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(
    dict(
        outcar=MapperAnnotation(mapper='.eigenvalues'),
        xml2=MapperAnnotation(mapper='.eigenvalues'),
    )
)


m_package.__init_metainfo__()
