from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass

from nomad.metainfo import Quantity, SchemaPackage
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

# Method-specific keys
DFT_XML_KEY = 'dft_xml'
DFT_OUTCAR_KEY = 'dft_outcar'

# Supplemental key for OUTCAR pseudopotential data to merge into XML structure
OUTCAR_PSEUDOPOT_KEY = 'outcar_pseudopot'


add_mapping_annotation(general.Simulation.m_def, XML_KEY, 'modeling')
add_mapping_annotation(general.Simulation.m_def, XML2_KEY, 'modeling')
add_mapping_annotation(general.Simulation.m_def, OUTCAR_KEY, '@')

# DFT-specific keys use same root paths as generic keys
add_mapping_annotation(general.Simulation.m_def, DFT_XML_KEY, 'modeling')
add_mapping_annotation(general.Simulation.m_def, DFT_OUTCAR_KEY, '@')


class Simulation(general.Simulation):
    add_mapping_annotation(general.Simulation.program, XML_KEY, '.generator')
    add_mapping_annotation(general.Simulation.program, OUTCAR_KEY, '.header')
    # dft method - use both generic and DFT-specific keys
    # Generic keys allow numerical_settings to populate into the created model_method
    add_mapping_annotation(
        general.Simulation.model_method,
        XML_KEY,
        '.parameters',
    )
    add_mapping_annotation(general.Simulation.model_method, OUTCAR_KEY, '.parameters')
    # DFT-specific keys ensure proper DFT type and avoid circular references
    add_mapping_annotation(
        model_method.DFT.m_def,
        DFT_XML_KEY,
        '.parameters',
    )
    add_mapping_annotation(model_method.DFT.m_def, DFT_OUTCAR_KEY, '.parameters')
    add_mapping_annotation(model_method.DFT.m_def, OUTCAR_PSEUDOPOT_KEY, '@')
    # NOTE: Pseudopotential annotations registered after class definition (line 203)
    # Ensures proper parser hierarchy: Simulation -> ModelMethod -> NumericalSettings
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


class DFT(model_method.DFT):
    add_mapping_annotation(model_method.DFT.xc, DFT_XML_KEY, '.@')
    add_mapping_annotation(model_method.DFT.xc, DFT_OUTCAR_KEY, '.@')

    # pseudopotential numerical settings (OUTCAR only)
    # Note: Pseudopotentials extracted only from OUTCAR, not vasprun.xml
    # vasprun.xml lacks detailed POTCAR metadata (LPAW, LULTRA, LEXCH, cutoffs, etc.)
    # Collection-level annotations are at module level (after Pseudopotential class)


class XCFunctional(model_method.XCFunctional):
    add_mapping_annotation(
        model_method.XCFunctional.components,
        DFT_XML_KEY,
        '.separator[?"@name"==\'electronic exchange-correlation\']',
    )
    add_mapping_annotation(
        model_method.XCFunctional.components,
        DFT_OUTCAR_KEY,
        ('get_xc_functionals', ['.@']),
    )


class XCComponent(model_method.XCComponent):
    add_mapping_annotation(
        model_method.XCComponent.canonical_label,
        DFT_XML_KEY,
        ('normalize_xc_label', ['.i[?"@name"==\'GGA\'] | [0].__value']),
    )
    add_mapping_annotation(
        model_method.XCComponent.canonical_label, DFT_OUTCAR_KEY, '.name'
    )


class ModelMethod(model_method.ModelMethod):
    # kspace numerical settings
    add_mapping_annotation(numerical_settings.KSpace.m_def, XML_KEY, 'modeling.kpoints')
    # TODO: Add KSpace mapping for OUTCAR k-points
    # add_mapping_annotation(numerical_settings.KSpace.m_def, OUTCAR_KEY, '@')

    # Note: Pseudopotential parsing is done only from OUTCAR auxiliary file
    # XML doesn't contain complete POTCAR metadata (LPAW, LULTRA, LEXCH,
    # cutoffs, SHA256). So we only create pseudopotentials when OUTCAR is available


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


class Pseudopotential(numerical_settings.Pseudopotential):
    """
    VASP-specific pseudopotential extension with ENMAX/ENMIN cutoff metadata.
    """

    sha256 = Quantity(
        type=str,
        description="""
        SHA256 hash of the POTCAR file for unique identification.
        This allows verification of pseudopotential provenance and ensures
        reproducibility by confirming the exact pseudopotential file used.
        """,
    )

    # Note: lpaw, lultra, lexch are extracted during parsing but not stored in schema.
    # They're used internally by the parser to derive `type` and `is_norm_conserving`.

    # Mapping annotations
    # OUTCAR provides complete POTCAR metadata
    add_mapping_annotation(
        numerical_settings.Pseudopotential.name, OUTCAR_KEY, '.titel'
    )
    add_mapping_annotation(
        numerical_settings.Pseudopotential.n_valence_electrons, OUTCAR_KEY, '.zval'
    )
    # vasprun.xml provides basic pseudopotential info (name and valence only)
    add_mapping_annotation(numerical_settings.Pseudopotential.name, XML_KEY, '.name')
    add_mapping_annotation(
        numerical_settings.Pseudopotential.n_valence_electrons,
        XML_KEY,
        '.n_valence_electrons',
    )
    add_mapping_annotation(
        numerical_settings.Pseudopotential.reference_configuration,
        OUTCAR_KEY,
        '.vrhfin',
    )
    add_mapping_annotation(
        numerical_settings.Pseudopotential.r_core, OUTCAR_KEY, '.rcore', unit='angstrom'
    )
    add_mapping_annotation(
        numerical_settings.Pseudopotential.l_max, OUTCAR_KEY, '.lmax'
    )
    add_mapping_annotation(
        numerical_settings.Pseudopotential.lm_max, OUTCAR_KEY, '.lmmax'
    )
    # PPCutoff subsections - ENMAX and ENMIN will be populated via custom parser logic
    # in get_pseudopotentials() rather than direct mapping annotations
    add_mapping_annotation(sha256, OUTCAR_KEY, '.sha256')

    # OUTCAR_PSEUDOPOT_KEY: Supplemental annotations for merging OUTCAR data into
    # existing pseudopotentials created by XML parsing (field-level only, no collection)
    add_mapping_annotation(
        numerical_settings.Pseudopotential.name, OUTCAR_PSEUDOPOT_KEY, '.titel'
    )
    add_mapping_annotation(
        numerical_settings.Pseudopotential.n_valence_electrons,
        OUTCAR_PSEUDOPOT_KEY,
        '.zval',
    )
    add_mapping_annotation(
        numerical_settings.Pseudopotential.reference_configuration,
        OUTCAR_PSEUDOPOT_KEY,
        '.vrhfin',
    )
    add_mapping_annotation(
        numerical_settings.Pseudopotential.r_core,
        OUTCAR_PSEUDOPOT_KEY,
        '.rcore',
        unit='angstrom',
    )
    add_mapping_annotation(
        numerical_settings.Pseudopotential.l_max, OUTCAR_PSEUDOPOT_KEY, '.lmax'
    )
    add_mapping_annotation(
        numerical_settings.Pseudopotential.lm_max, OUTCAR_PSEUDOPOT_KEY, '.lmmax'
    )
    add_mapping_annotation(sha256, OUTCAR_PSEUDOPOT_KEY, '.sha256')


# Pseudopotential collection mapping - must be after Pseudopotential class definition
# These annotations are added in the context of ModelMethod to ensure the parser
# creates Pseudopotential subsections in model_method.numerical_settings
# Use both generic and DFT-specific keys so pseudopotentials populate into DFT objects
add_mapping_annotation(
    Pseudopotential.m_def,
    OUTCAR_KEY,
    ('get_pseudopotentials', ['@.pseudopotentials']),
)
add_mapping_annotation(
    Pseudopotential.m_def,
    DFT_OUTCAR_KEY,
    ('get_pseudopotentials', ['@.pseudopotentials']),
)
add_mapping_annotation(
    Pseudopotential.m_def,
    OUTCAR_PSEUDOPOT_KEY,
    ('get_pseudopotentials', ['@.pseudopotentials']),
)

# Note: XML pseudopotential collection annotation is in ModelMethod class above
# to avoid duplicate registration conflicts


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


try:
    m_package.__init_metainfo__()
except Exception:
    pass
