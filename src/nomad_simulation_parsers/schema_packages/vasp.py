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
)

m_package = SchemaPackage()


general.Simulation.m_def.m_annotations.setdefault(MAPPING_ANNOTATION_KEY, {}).update(
    dict(
        xml=MapperAnnotation(mapper='modeling'),
        xml2=MapperAnnotation(mapper='modeling'),
        outcar=MapperAnnotation(mapper='@'),
    )
)


class Simulation(general.Simulation):
    general.Simulation.program.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(
        dict(
            xml=MapperAnnotation(mapper='.generator'),
            outcar=MapperAnnotation(mapper='.header'),
        )
    )
    # dft method
    model_method.DFT.m_def.m_annotations.setdefault(MAPPING_ANNOTATION_KEY, {}).update(
        dict(
            xml=MapperAnnotation(
                mapper='.parameters.separator[?"@name"==\'electronic\']'
            ),
            outcar=MapperAnnotation(mapper='parameters'),
        )
    )
    general.Simulation.model_system.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(
        dict(
            xml=MapperAnnotation(mapper='.calculation'),
            outcar=MapperAnnotation(mapper='.calculation'),
        )
    )
    general.Simulation.outputs.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(
        dict(
            xml=MapperAnnotation(mapper='.calculation'),
            xml2=MapperAnnotation(mapper='.calculation'),
            outcar=MapperAnnotation(mapper='.calculation'),
        )
    )


class Program(general.Program):
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
    ).update(
        dict(xml=MapperAnnotation(mapper='.i[?"@name"==\'platform\'] | [0].__value'))
    )


class DFT(model_method.DFT):
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


class XCFunctional(model_method.XCFunctional):
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


class ModelMethod(model_method.ModelMethod):
    # kspace numerical settings
    numerical_settings.KSpace.m_def.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(xml=MapperAnnotation(mapper='modeling.kpoints')))


class KSpace(numerical_settings.KSpace):
    numerical_settings.KSpace.k_mesh.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(dict(xml=MapperAnnotation(mapper='.@')))


class KMesh(numerical_settings.KMesh):
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
            xml=MapperAnnotation(
                mapper='.generation.v[?"@name"==\'shift\'] | [0].__value'
            )
        )
    )
    numerical_settings.KMesh.points.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(
        dict(
            xml=MapperAnnotation(
                mapper=(
                    'reshape_array',
                    ['.varray[?"@name"==\'kpointlist\'].v | [0]'],
                    dict(shape_rest=(3)),
                )
            )
        )
    )

    numerical_settings.KMesh.weights.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(
        dict(
            xml=MapperAnnotation(
                mapper=(
                    'reshape_array',
                    ['.varray[?"@name"==\'weights\'].v | [0]'],
                    dict(shape_rest=()),
                )
            )
        )
    )


class ModelSystem(model_system.ModelSystem):
    # atomic cell
    model_system.AtomicCell.m_def.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(
        dict(
            xml=MapperAnnotation(mapper='.structure'),
            outcar=MapperAnnotation(mapper='.@'),
        )
    )
    model_system.ModelSystem.positions.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(
        dict(
            xml=MapperAnnotation(
                mapper=(
                    'reshape_array',
                    ['.structure.varray.v'],
                    dict(shape_rest=(3,)),
                ),
                unit='angstrom',
            ),
            outcar=MapperAnnotation(
                mapper='.positions_forces', unit='angstrom', search='@ | [0]'
            ),
        )
    )


class AtomicCell(model_system.AtomicCell):
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


class Outputs(outputs.Outputs):
    outputs.Outputs.total_energies.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(
        dict(
            xml=MapperAnnotation(mapper='.energy'),
            outcar=MapperAnnotation(mapper='.energies'),
        )
    )
    outputs.Outputs.total_forces.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(
        dict(
            xml=MapperAnnotation(mapper=('get_forces', ['.@'])),
            outcar=MapperAnnotation(mapper=('get_forces', ['.@'])),
        )
    )
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
    outputs.Outputs.electronic_dos.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(
        dict(
            xml=MapperAnnotation(mapper=('get_dos', ['.dos.total.array.set.set.r'])),
        )
    )


class TotalEnergy(properties.energies.TotalEnergy):
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


class EnergyContribution(properties.energies.EnergyContribution):
    properties.energies.EnergyContribution.name.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(
        dict(
            xml=MapperAnnotation(mapper='."@name"'),
            outcar=MapperAnnotation(mapper='.name'),
        )
    )


class TotalForce(properties.forces.TotalForce):
    properties.forces.TotalForce.value.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(
        dict(
            xml=MapperAnnotation(mapper='.forces', unit='eV/angstrom'),
            outcar=MapperAnnotation(mapper='.forces', unit='eV/angstrom'),
        )
    )


class ElectronicEigenvalues(outputs.ElectronicEigenvalues):
    outputs.ElectronicEigenvalues.n_bands.m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(
        dict(
            xml=MapperAnnotation(mapper='length(.array.set.set.set[0].r)'),
            xml2=MapperAnnotation(mapper='length(.array.set.set.set[0].r)'),
            outcar=MapperAnnotation(mapper='.n_bands'),
        )
    )

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


outputs.ElectronicDensityOfStates.value.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(
    dict(
        xml=MapperAnnotation(mapper='.value', unit='1/eV'),
    )
)
outputs.ElectronicDensityOfStates.energies.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(
    dict(
        xml=MapperAnnotation(mapper='.energies', unit='eV'),
    )
)

try:
    m_package.__init_metainfo__()
except Exception:
    pass
