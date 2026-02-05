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


add_mapping_annotation(general.Simulation.m_def, XML_KEY, 'modeling')
add_mapping_annotation(general.Simulation.m_def, XML2_KEY, 'modeling')
add_mapping_annotation(general.Simulation.m_def, OUTCAR_KEY, '@')


class Simulation(general.Simulation):
    add_mapping_annotation(general.Simulation.program, XML_KEY, 'modeling.generator')
    add_mapping_annotation(general.Simulation.program, OUTCAR_KEY, '.header')
    # dft method
    add_mapping_annotation(
        model_method.DFT.m_def,
        XML_KEY,
        'modeling.parameters.separator[?"@name"==\'electronic\']',
    )
    add_mapping_annotation(model_method.DFT.m_def, OUTCAR_KEY, 'parameters')
    add_mapping_annotation(
        general.Simulation.model_system, XML_KEY, 'modeling.calculation'
    )
    add_mapping_annotation(general.Simulation.model_system, OUTCAR_KEY, '.calculation')
    add_mapping_annotation(general.Simulation.outputs, XML_KEY, 'modeling.calculation')
    add_mapping_annotation(general.Simulation.outputs, XML2_KEY, 'modeling.calculation')
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


class DFT(model_method.DFT):
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
    """
    VASP-specific pseudopotential metadata extracted from POTCAR headers in
    OUTCAR and vasprun.xml.

    Extends base Pseudopotential class with VASP-specific fields including
    SHA256 hash, LPAW/LULTRA flags, and LEXCH code. Basic fields are mapped
    directly from POTCAR data. Type and related properties are derived using
    field-level transformers.
    """

    import numpy as np
    from nomad.metainfo import Quantity

    sha256 = Quantity(
        type=str,
        description="""
        SHA256 hash of the POTCAR file content. Uniquely identifies the
        pseudopotential file and enables verification that the correct POTCAR
        was used. This hash can be matched against pseudopotential library
        databases for automatic library detection.
        """,
    )

    # VASP-specific raw fields from POTCAR headers
    lpaw = Quantity(
        type=str,
        description='LPAW flag from POTCAR (T=PAW, F=other)',
    )
    lultra = Quantity(
        type=str,
        description='LULTRA flag from POTCAR (T=ultrasoft, F=other)',
    )
    lexch = Quantity(
        type=str,
        description='LEXCH exchange-correlation code from POTCAR',
    )
    enmax = Quantity(
        type=np.float64,
        unit='eV',
        description='ENMAX: maximum kinetic energy cutoff recommendation',
    )
    enmin = Quantity(
        type=np.float64,
        unit='eV',
        description='ENMIN: minimum kinetic energy cutoff recommendation',
    )

    # Direct field mappings from raw POTCAR keys using relative paths
    # String fields - direct mapping
    add_mapping_annotation(numerical_settings.Pseudopotential.name, OUTCAR_KEY, '.TITEL')
    add_mapping_annotation(
        numerical_settings.Pseudopotential.reference_configuration,
        OUTCAR_KEY,
        '.VRHFIN',
    )
    add_mapping_annotation(sha256, OUTCAR_KEY, '.SHA256')

    # Numeric fields - use type conversion transformers
    add_mapping_annotation(
        numerical_settings.Pseudopotential.n_valence_electrons,
        OUTCAR_KEY,
        ('to_float', ['.ZVAL']),
    )
    add_mapping_annotation(
        numerical_settings.Pseudopotential.r_core,
        OUTCAR_KEY,
        ('to_float', ['.RCORE']),
        unit='angstrom',
    )
    add_mapping_annotation(
        numerical_settings.Pseudopotential.l_max,
        OUTCAR_KEY,
        ('to_int', ['.LMAX'])
    )
    add_mapping_annotation(
        numerical_settings.Pseudopotential.lm_max,
        OUTCAR_KEY,
        ('to_int', ['.LMMAX'])
    )

    # VASP-specific fields
    add_mapping_annotation(lpaw, OUTCAR_KEY, '.LPAW')
    add_mapping_annotation(lultra, OUTCAR_KEY, '.LULTRA')
    add_mapping_annotation(lexch, OUTCAR_KEY, '.LEXCH')
    add_mapping_annotation(
        enmax,
        OUTCAR_KEY,
        ('to_float', ['.ENMAX'])
    )
    add_mapping_annotation(
        enmin,
        OUTCAR_KEY,
        ('to_float', ['.ENMIN'])
    )

    # Derived fields using field-level transformers (per-instance)
    add_mapping_annotation(
        numerical_settings.Pseudopotential.type,
        OUTCAR_KEY,
        ('derive_pp_type', ['.LPAW', '.LULTRA', '.TITEL'])
    )
    add_mapping_annotation(
        numerical_settings.Pseudopotential.is_norm_conserving,
        OUTCAR_KEY,
        ('derive_is_norm_conserving', ['.type'])
    )
    add_mapping_annotation(
        numerical_settings.Pseudopotential.is_gw_optimized,
        OUTCAR_KEY,
        ('derive_is_gw_optimized', ['.TITEL'])
    )

    # XML field annotations for basic metadata (name and valence electrons only)
    add_mapping_annotation(numerical_settings.Pseudopotential.name, XML_KEY, '.name')
    add_mapping_annotation(
        numerical_settings.Pseudopotential.n_valence_electrons,
        XML_KEY,
        '.n_valence_electrons',
    )


# Collection annotations for Pseudopotential - added after class definition
# so that Pseudopotential.m_def references the custom VASP class
# From vasprun.xml
add_mapping_annotation(
    Pseudopotential.m_def,
    XML_KEY,
    ('get_pseudopotentials_xml', ['atominfo.array']),
)
# From OUTCAR - direct path, field annotations handle mapping
add_mapping_annotation(
    Pseudopotential.m_def,
    OUTCAR_KEY,
    '@.pseudopotentials',
)


try:
    m_package.__init_metainfo__()
except Exception:
    pass
