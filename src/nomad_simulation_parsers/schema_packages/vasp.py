from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass

from nomad.metainfo import SchemaPackage
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

XML_KEY = 'xml'
XML2_KEY = 'xml2'
OUTCAR_KEY = 'outcar'
KPOINTS_XML = 'kpoints_xml'
PP_XML = 'xml_pseudopotentials'
PP_OUT = 'outcar_pseudopotentials'


add_mapping_annotation(general.Simulation.m_def, XML_KEY, 'modeling')
add_mapping_annotation(general.Simulation.m_def, XML2_KEY, 'modeling')
add_mapping_annotation(general.Simulation.m_def, OUTCAR_KEY, '@')
add_mapping_annotation(general.Simulation.m_def, KPOINTS_XML, 'modeling')
add_mapping_annotation(general.Simulation.m_def, PP_XML, 'modeling')
add_mapping_annotation(general.Simulation.m_def, PP_OUT, '@')


class Simulation(general.Simulation):
    add_mapping_annotation(general.Simulation.program, XML_KEY, '.generator')
    add_mapping_annotation(general.Simulation.program, OUTCAR_KEY, '.header')
    add_mapping_annotation(
        model_method.DFT.m_def,
        XML_KEY,
        'modeling.parameters.separator[?"@name"==\'electronic\']',
    )
    add_mapping_annotation(
        model_method.DFT.m_def,
        KPOINTS_XML,
        '.parameters.separator[?"@name"==\'electronic\']',
    )
    add_mapping_annotation(
        model_method.DFT.m_def,
        PP_XML,
        '.parameters.separator[?"@name"==\'electronic\']',
    )
    add_mapping_annotation(general.ModelMethod.m_def, KPOINTS_XML, '.@')
    add_mapping_annotation(general.ModelMethod.m_def, PP_XML, '.@')
    add_mapping_annotation(general.ModelMethod.m_def, PP_OUT, '.@')
    add_mapping_annotation(general.Simulation.model_system, XML_KEY, '.calculation')
    add_mapping_annotation(general.Simulation.model_system, OUTCAR_KEY, '.calculation')
    add_mapping_annotation(general.Simulation.outputs, XML_KEY, '.calculation')
    add_mapping_annotation(general.Simulation.outputs, XML2_KEY, '.calculation')
    add_mapping_annotation(general.Simulation.outputs, OUTCAR_KEY, '.calculation')


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


class ModelMethod(general.ModelMethod):
    add_mapping_annotation(
        numerical_settings.Pseudopotential.m_def,
        PP_XML,
        '.atominfo.array[?"@name"==\'atomtypes\'] | [0].set.rc',
    )


class DFT(model_method.DFT):
    # KSpace: Changed from absolute path 'kpoints' to relative '.kpoints' to enable
    # mapper tree traversal. The mapper builder needs a connected chain of annotations
    # from root → ModelMethod → numerical_settings → KSpace → k_mesh
    add_mapping_annotation(numerical_settings.KSpace.m_def, KPOINTS_XML, '.kpoints')
    add_mapping_annotation(model_method.DFT.xc, XML_KEY, '.@')
    add_mapping_annotation(model_method.DFT.xc, OUTCAR_KEY, '.@')


class XCFunctional(model_method.XCFunctional):
    add_mapping_annotation(
        model_method.XCFunctional.components,
        XML_KEY,
        '.separator[?"@name"==\'electronic exchange-correlation\']',
    )
    add_mapping_annotation(
        model_method.XCFunctional.components, OUTCAR_KEY, ('get_xc_functionals', ['.@'])
    )


class XCComponent(model_method.XCComponent):
    add_mapping_annotation(
        model_method.XCComponent.canonical_label,
        XML_KEY,
        '.i[?"@name"==\'GGA\'] | [0].__value',
    )
    add_mapping_annotation(
        model_method.XCComponent.canonical_label, OUTCAR_KEY, '.name'
    )


class KSpace(numerical_settings.KSpace):
    add_mapping_annotation(numerical_settings.KSpace.k_mesh, KPOINTS_XML, '.@')


class KMesh(numerical_settings.KMesh):
    add_mapping_annotation(
        numerical_settings.KMesh.grid,
        KPOINTS_XML,
        '.generation.v[?"@name"==\'divisions\'] | [0].__value',
    )
    add_mapping_annotation(
        numerical_settings.KMesh.offset,
        KPOINTS_XML,
        '.generation.v[?"@name"==\'shift\'] | [0].__value',
    )
    add_mapping_annotation(
        numerical_settings.KMesh.points,
        KPOINTS_XML,
        (
            'reshape_array',
            ['.varray[?"@name"==\'kpointlist\'].v | [0]'],
            dict(shape_rest=(3)),
        ),
    )
    add_mapping_annotation(
        numerical_settings.KMesh.weights,
        KPOINTS_XML,
        (
            'reshape_array',
            ['.varray[?"@name"==\'weights\'].v | [0]'],
            dict(shape_rest=()),
        ),
    )


class ModelSystem(model_system.ModelSystem):
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


class TotalEnergy(properties.energies.TotalEnergy):
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
    add_mapping_annotation(
        outputs.ElectronicEigenvalues.n_levels,
        XML_KEY,
        'length(.array.set.set.set[0].r)',
    )
    add_mapping_annotation(
        outputs.ElectronicEigenvalues.n_levels,
        XML2_KEY,
        'length(.array.set.set.set[0].r)',
    )
    add_mapping_annotation(
        outputs.ElectronicEigenvalues.n_levels, OUTCAR_KEY, '.n_bands'
    )

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


class Pseudopotential(numerical_settings.Pseudopotential):
    """VASP-specific pseudopotential with POTCAR metadata."""

    import numpy as np
    from nomad.metainfo import Quantity

    sha256 = Quantity(
        type=str,
        description='SHA256 hash of POTCAR file for library identification',
    )
    add_mapping_annotation(sha256, PP_OUT, '.SHA256')

    lpaw = Quantity(type=str, description='LPAW flag (T=PAW, F=other)')
    add_mapping_annotation(lpaw, PP_OUT, '.LPAW')

    lultra = Quantity(type=str, description='LULTRA flag (T=ultrasoft, F=other)')
    add_mapping_annotation(lultra, PP_OUT, '.LULTRA')

    lexch = Quantity(type=str, description='LEXCH exchange-correlation code')
    add_mapping_annotation(lexch, PP_OUT, '.LEXCH')

    enmax = Quantity(
        type=np.float64, unit='eV', description='ENMAX cutoff recommendation'
    )
    add_mapping_annotation(enmax, PP_OUT, '.ENMAX')

    enmin = Quantity(
        type=np.float64, unit='eV', description='ENMIN cutoff recommendation'
    )
    add_mapping_annotation(enmin, PP_OUT, '.ENMIN')

    add_mapping_annotation(numerical_settings.Pseudopotential.name, PP_OUT, '.TITEL')
    add_mapping_annotation(numerical_settings.Pseudopotential.name, PP_XML, '.c[4]')
    add_mapping_annotation(
        numerical_settings.Pseudopotential.reference_configuration, PP_OUT, '.VRHFIN'
    )
    add_mapping_annotation(
        numerical_settings.Pseudopotential.n_valence_electrons,
        PP_OUT,
        '.ZVAL',
    )
    add_mapping_annotation(
        numerical_settings.Pseudopotential.n_valence_electrons,
        PP_XML,
        '.c[3]',
    )
    add_mapping_annotation(
        numerical_settings.Pseudopotential.r_core,
        PP_OUT,
        '.RCORE',
        unit='angstrom',
    )
    add_mapping_annotation(numerical_settings.Pseudopotential.l_max, PP_OUT, '.LMAX')
    add_mapping_annotation(numerical_settings.Pseudopotential.lm_max, PP_OUT, '.LMMAX')
    add_mapping_annotation(
        numerical_settings.Pseudopotential.type,
        PP_OUT,
        ('derive_pp_type', ['.LPAW', '.LULTRA', '.TITEL']),
    )
    add_mapping_annotation(
        numerical_settings.Pseudopotential.is_norm_conserving,
        PP_OUT,
        ('derive_is_norm_conserving', ['.type']),
    )
    add_mapping_annotation(
        numerical_settings.Pseudopotential.is_gw_optimized,
        PP_OUT,
        ('derive_is_gw_optimized', ['.TITEL']),
    )


try:
    m_package.__init_metainfo__()
except Exception as e:
    print(f'[ERROR] Failed to initialize VASP schema package: {e}', flush=True)
    import traceback

    traceback.print_exc()
