from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass

from nomad.datamodel import ArchiveSection
from nomad.datamodel.hdf5 import HDF5Dataset
from nomad.metainfo import Quantity, SchemaPackage, SubSection
from nomad_simulations.schema_packages import (
    general,
    model_method,
    model_system,
    numerical_settings,
    outputs,
    properties,
)

from nomad_simulation_parsers.schema_packages.utils import add_mapping_annotation

m_package = SchemaPackage()

XML_KEY = 'vasp_xml'
XML2_KEY = 'vasp_xml2'
OUTCAR_KEY = 'vasp_outcar'
CHGCAR_KEY = 'vasp_chgcar'


add_mapping_annotation(general.Simulation.m_def, XML_KEY, 'modeling')
add_mapping_annotation(general.Simulation.m_def, XML2_KEY, 'modeling')
add_mapping_annotation(general.Simulation.m_def, OUTCAR_KEY, '@')
add_mapping_annotation(general.Simulation.m_def, CHGCAR_KEY, '@')


class ChargeDensity(ArchiveSection):
    value_h5_dataset = Quantity(type=HDF5Dataset)
    add_mapping_annotation(value_h5_dataset, CHGCAR_KEY, '.@')


class VASPOutputs(outputs.Outputs):
    charge_density = SubSection(sub_section=ChargeDensity.m_def, repeats=True)
    add_mapping_annotation(charge_density, CHGCAR_KEY, '.values')


class Simulation(general.Simulation):
    add_mapping_annotation(general.Simulation.program, XML_KEY, '.generator')
    add_mapping_annotation(general.Simulation.program, OUTCAR_KEY, '.header')
    # dft method
    add_mapping_annotation(
        model_method.DFT.m_def,
        XML_KEY,
        '.parameters.separator[?"@name"==\'electronic\']',
    )
    add_mapping_annotation(model_method.DFT.m_def, OUTCAR_KEY, 'parameters')
    add_mapping_annotation(general.Simulation.model_system, XML_KEY, '.calculation')
    add_mapping_annotation(general.Simulation.model_system, OUTCAR_KEY, '.calculation')
    add_mapping_annotation(general.Simulation.outputs, XML_KEY, '.calculation')
    add_mapping_annotation(general.Simulation.outputs, XML2_KEY, '.calculation')
    add_mapping_annotation(general.Simulation.outputs, OUTCAR_KEY, '.calculation')
    # TODO: make update_mode merge@last when mapping parser is updated
    add_mapping_annotation(VASPOutputs.m_def, CHGCAR_KEY, '.@', update_mode='append')


class Program(general.Program):
    add_mapping_annotation(
        general.Program.name, XML_KEY, '.i[?"@name"==\'program\'] | [0].__value'
    )
    add_mapping_annotation(
        general.Program.version, XML_KEY, '.i[?"@name"==\'version\'] | [0].__value'
    )
    add_mapping_annotation(general.Program.version, OUTCAR_KEY, ('get_version', ['.@']))
    add_mapping_annotation(
        general.Program.compilation_host,
        XML_KEY,
        '.i[?"@name"==\'platform\'] | [0].__value',
    )


# class DFT(model_method.DFT):
#     model_method.DFT.xc_functionals.m_annotations.setdefault(
#         MAPPING_ANNOTATION_KEY, {}
#     ).update(
#         dict(
#             xml=MapperAnnotation(
#                 mapper='.separator[?"@name"==\'electronic exchange-correlation\']'
#             ),
#             outcar=MapperAnnotation(mapper=('get_xc_functionals', ['.@'])),
#         )
#     )
#     model_method.DFT.exact_exchange_mixing_factor.m_annotations.setdefault(
#         MAPPING_ANNOTATION_KEY, {}
#     ).update(
#         dict(
#             xml=MapperAnnotation(
#                 mapper=(
#                     'mix_alpha',
#                     [
#                         '.i[?"@name"==\'HFALPHA\'] | [0].__value',
#                         '.i[?"@name"==\'LHFCALC\'] | [0].__value',
#                     ],
#                 )
#             )
#         )
#     )  # TODO convert vasp bool


# class XCFunctional(model_method.XCFunctional):
#     model_method.XCFunctional.libxc_name.m_annotations.setdefault(
#         MAPPING_ANNOTATION_KEY, {}
#     ).update(
#         dict(
#             xml=MapperAnnotation(
#                 # TODO add LDA & mGGA, convert_xc
#                 mapper='.i[?"@name"==\'GGA\'] | [0].__value'
#             ),
#             outcar=MapperAnnotation(mapper='.name'),
#         )
#     )


class ModelMethod(model_method.ModelMethod):
    # kspace numerical settings
    add_mapping_annotation(numerical_settings.KSpace.m_def, XML_KEY, 'modeling.kpoints')


class KSpace(numerical_settings.KSpace):
    add_mapping_annotation(numerical_settings.KSpace.k_mesh, XML_KEY, '.@')


class KMesh(numerical_settings.KMesh):
    add_mapping_annotation(
        numerical_settings.KMesh.grid,
        XML_KEY,
        '.generation.v[?"@name"==\'divisions\'] | [0].__value',
    )
    add_mapping_annotation(
        numerical_settings.KMesh.offset,
        XML_KEY,
        '.generation.v[?"@name"==\'shift\'] | [0].__value',
    )
    add_mapping_annotation(
        numerical_settings.KMesh.points,
        XML_KEY,
        (
            'reshape_array',
            ['.varray[?"@name"==\'kpointlist\'].v | [0]'],
            dict(shape_rest=(3)),
        ),
    )
    add_mapping_annotation(
        numerical_settings.KMesh.weights,
        XML_KEY,
        (
            'reshape_array',
            ['.varray[?"@name"==\'weights\'].v | [0]'],
            dict(shape_rest=()),
        ),
    )


class ModelSystem(model_system.ModelSystem):
    # atomic cell
    add_mapping_annotation(model_system.Representation.m_def, XML_KEY, '.structure')
    add_mapping_annotation(model_system.Representation.m_def, OUTCAR_KEY, '.@')
    add_mapping_annotation(
        model_system.ModelSystem.positions,
        XML_KEY,
        (
            'reshape_array',
            ['.structure.varray.v'],
            dict(shape_rest=(3,)),
        ),
        unit='angstrom',
    )
    add_mapping_annotation(
        model_system.ModelSystem.positions,
        OUTCAR_KEY,
        '.positions_forces',
        unit='angstrom',
        search='@ | [0]',
    )


class Representation(model_system.Representation):
    add_mapping_annotation(
        model_system.Representation.lattice_vectors,
        XML_KEY,
        '.structure.varray[?"@name"==\'basis\'] | [0].v',
        unit='angstrom',
    )
    add_mapping_annotation(
        model_system.Representation.lattice_vectors,
        OUTCAR_KEY,
        '.lattice_vectors',
        unit='angstrom',
        search='@ | [0]',
    )


class Outputs(outputs.Outputs):
    add_mapping_annotation(outputs.Outputs.total_energies, XML_KEY, '.energy')
    add_mapping_annotation(outputs.Outputs.total_energies, OUTCAR_KEY, '.energies')
    add_mapping_annotation(
        outputs.Outputs.total_forces, XML_KEY, ('get_forces', ['.@'])
    )
    add_mapping_annotation(
        outputs.Outputs.total_forces, OUTCAR_KEY, ('get_forces', ['.@'])
    )
    add_mapping_annotation(
        outputs.Outputs.electronic_eigenvalues,
        XML_KEY,
        ('get_eigenvalues', ['eigenvalues']),
    )
    add_mapping_annotation(
        outputs.Outputs.electronic_eigenvalues,
        XML2_KEY,
        ('get_eigenvalues', ['eigenvalues']),
    )
    add_mapping_annotation(
        outputs.Outputs.electronic_eigenvalues,
        OUTCAR_KEY,
        ('get_eigenvalues', ['.eigenvalues', 'parameters']),
    )
    add_mapping_annotation(
        outputs.Outputs.scf_steps,
        XML_KEY,
        ('get_scf_steps', ['.@']),
    )
    add_mapping_annotation(
        outputs.Outputs.scf_steps,
        XML2_KEY,
        ('get_scf_steps', ['.@']),
    )
    add_mapping_annotation(
        outputs.Outputs.scf_steps, OUTCAR_KEY, ('get_scf_steps', ['.@'])
    )


class SCFSteps(outputs.SCFSteps):
    add_mapping_annotation(outputs.SCFSteps.energies_total, XML_KEY, '.energies_total')
    add_mapping_annotation(outputs.SCFSteps.energies_total, XML2_KEY, '.energies_total')
    add_mapping_annotation(
        outputs.SCFSteps.energies_total, OUTCAR_KEY, '.energies_total'
    )
    add_mapping_annotation(
        outputs.SCFSteps.delta_energies_total, XML_KEY, '.delta_energies_total'
    )
    add_mapping_annotation(
        outputs.SCFSteps.delta_energies_total, XML2_KEY, '.delta_energies_total'
    )
    add_mapping_annotation(
        outputs.SCFSteps.delta_energies_total, OUTCAR_KEY, '.delta_energies_total'
    )
    add_mapping_annotation(outputs.SCFSteps.durations, XML_KEY, '.durations')
    add_mapping_annotation(outputs.SCFSteps.durations, XML2_KEY, '.durations')
    add_mapping_annotation(outputs.SCFSteps.durations, OUTCAR_KEY, '.durations')


class TotalEnergy(properties.energies.TotalEnergy):
    # value is already defined in TotalEnergy since they use the same def
    # get_energy function should be able to handle extraction from both sources
    add_mapping_annotation(
        properties.energies.TotalEnergy.value,
        XML_KEY,
        (
            'get_data',
            ['.@'],
            dict(path='.i[?"@name"==\'e_fr_energy\'] | [0].__value'),
        ),
        unit='eV',
    )
    add_mapping_annotation(
        properties.energies.TotalEnergy.value,
        OUTCAR_KEY,
        ('get_data', ['.@'], dict(path='.energy_total')),
        unit='eV',
    )
    add_mapping_annotation(
        properties.energies.TotalEnergy.contributions,
        XML_KEY,
        (
            'get_energy_contributions',
            ['.i'],
            dict(exclude=['e_fr_energy']),
        ),
    )
    add_mapping_annotation(
        properties.energies.TotalEnergy.contributions,
        OUTCAR_KEY,
        (
            'get_energy_contributions',
            ['.@'],
            dict(exclude=['energy_total']),
        ),
    )


class BaseEnergy(properties.energies.BaseEnergy):
    add_mapping_annotation(properties.energies.BaseEnergy.name, XML_KEY, '.@name')
    add_mapping_annotation(properties.energies.BaseEnergy.name, OUTCAR_KEY, '.name')


class TotalForce(properties.forces.TotalForce):
    add_mapping_annotation(
        properties.forces.TotalForce.value, XML_KEY, '.forces', unit='eV/angstrom'
    )
    add_mapping_annotation(
        properties.forces.TotalForce.value, OUTCAR_KEY, '.forces', unit='eV/angstrom'
    )


# TODO: check whether this section is k-dependent
class ElectronicEigenvalues(outputs.ElectronicEigenvalues):
    # outputs.ElectronicEigenvalues.n_bands.m_annotations.setdefault(
    #     MAPPING_ANNOTATION_KEY, {}
    # ).update(
    #     dict(
    #         xml=MapperAnnotation(mapper='length(.array.set.set.set[0].r)'),
    #         xml2=MapperAnnotation(mapper='length(.array.set.set.set[0].r)'),
    #         outcar=MapperAnnotation(mapper='.n_bands'),
    #     )
    # )

    # TODO This only works for non-spin pol
    add_mapping_annotation(
        outputs.ElectronicEigenvalues.occupation, OUTCAR_KEY, '.occupations'
    )
    add_mapping_annotation(
        outputs.ElectronicEigenvalues.occupation, XML2_KEY, '.occupations'
    )
    add_mapping_annotation(
        outputs.ElectronicEigenvalues.value, OUTCAR_KEY, '.eigenvalues'
    )
    add_mapping_annotation(
        outputs.ElectronicEigenvalues.value, XML2_KEY, '.eigenvalues'
    )


try:
    m_package.__init_metainfo__()
except Exception:
    pass
