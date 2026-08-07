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
    variables,
)

from nomad_simulation_parsers.schema_packages.utils import add_mapping_annotation

m_package = SchemaPackage()

XML_KEY = 'vasp_xml'
XML2_KEY = 'vasp_xml2'
OUTCAR_KEY = 'vasp_outcar'
CHGCAR_KEY = 'vasp_chgcar'
DOSCAR_KEY = 'vasp_doscar'


add_mapping_annotation(general.Simulation.m_def, XML_KEY, 'modeling')
add_mapping_annotation(general.Simulation.m_def, XML2_KEY, 'modeling')
add_mapping_annotation(general.Simulation.m_def, OUTCAR_KEY, '@')
add_mapping_annotation(general.Simulation.m_def, CHGCAR_KEY, '@')
add_mapping_annotation(general.Simulation.m_def, DOSCAR_KEY, '@')


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
    add_mapping_annotation(
        general.Simulation.outputs, DOSCAR_KEY, '.@', update_mode='append'
    )


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


# TODO: map hybrid exact-exchange mixing into `XCFunctional.global_exact_exchange`
# (formerly `DFT.exact_exchange_mixing_factor`). VASP exposes it via HFALPHA /
# LHFCALC, and `VasprunParser.mix_alpha` is the intended transformer. Deferred:
# needs a hybrid test fixture (AgAc_relax is PBE).
class XCFunctional(model_method.XCFunctional):
    # Set only the canonical functional name; `XCFunctional.normalize` expands it
    # into complete LibXC `components` (exchange + correlation, with family/kind)
    # from the shared alias/registry tables, and `DFT.normalize` derives
    # `jacobs_ladder`. The parser therefore does not build components itself. The
    # `get_functional_key` transformer gathers the VASP XC tags itself (flat
    # OUTCAR parameters, or by walking the vasprun `<parameters>` tree).
    add_mapping_annotation(
        model_method.XCFunctional.functional_key,
        XML_KEY,
        ('get_functional_key', ['.@']),
    )
    add_mapping_annotation(
        model_method.XCFunctional.functional_key,
        OUTCAR_KEY,
        ('get_functional_key', ['.@']),
    )


class DFT(model_method.DFT):
    # The binding only needs to materialize the `xc` subsection so the child
    # `functional_key` mapper runs -- `get_functional_key` ignores the bound node
    # and gathers the XC tags itself (walking the full vasprun `<parameters>`
    # tree, or the flat OUTCAR parameters). So bind to the current node (`.@`,
    # always present) for both sources rather than a specific nested separator
    # that some vasprun variants omit.
    add_mapping_annotation(model_method.DFT.xc, XML_KEY, '.@')
    add_mapping_annotation(model_method.DFT.xc, OUTCAR_KEY, '.@')


class ModelMethod(model_method.ModelMethod):
    # kspace numerical settings
    add_mapping_annotation(numerical_settings.KSpace.m_def, XML_KEY, 'modeling.kpoints')
    add_mapping_annotation(
        numerical_settings.SelfConsistency.m_def,
        XML_KEY,
        'modeling.parameters.separator[?"@name"==\'electronic\'] | [0]',
    )
    add_mapping_annotation(
        numerical_settings.SelfConsistency.m_def, OUTCAR_KEY, 'parameters'
    )


class KSpace(numerical_settings.KSpace):
    add_mapping_annotation(numerical_settings.KSpace.k_mesh, XML_KEY, '.@')


class SelfConsistency(numerical_settings.SelfConsistency):
    add_mapping_annotation(
        numerical_settings.SelfConsistency.threshold_change,
        XML_KEY,
        (
            'modeling.parameters.separator[?"@name"==\'electronic\'] '
            '| [0].i[?"@name"==\'EDIFF\'] | [0].__value'
        ),
        unit='eV',
    )
    add_mapping_annotation(
        numerical_settings.SelfConsistency.threshold_change,
        OUTCAR_KEY,
        'parameters.EDIFF',
        unit='eV',
    )


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
    add_mapping_annotation(
        model_system.AtomsState.m_def,
        XML_KEY,
        ('get_atoms', ['modeling.atominfo.array']),
    )
    add_mapping_annotation(
        model_system.AtomsState.m_def,
        OUTCAR_KEY,
        ('get_atoms', ['ions_per_type', 'species']),
    )
    add_mapping_annotation(
        model_system.ModelSystem.positions,
        XML_KEY,
        (
            'get_positions',
            [
                '.structure.varray.v',
                '.structure.crystal.varray[?"@name"==\'basis\'] | [0].v',
            ],
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
    add_mapping_annotation(
        model_system.ModelSystem.lattice_vectors,
        XML_KEY,
        '.structure.crystal.varray[?"@name"==\'basis\'] | [0].v',
        unit='angstrom',
    )
    add_mapping_annotation(
        model_system.ModelSystem.lattice_vectors,
        OUTCAR_KEY,
        '.lattice_vectors',
        unit='angstrom',
        search='@ | [0]',
    )
    add_mapping_annotation(
        model_system.ModelSystem.periodic_boundary_conditions,
        XML_KEY,
        ('get_periodic_boundary_conditions', []),
    )
    add_mapping_annotation(
        model_system.ModelSystem.periodic_boundary_conditions,
        OUTCAR_KEY,
        ('get_periodic_boundary_conditions', []),
    )


class AtomsState(model_system.AtomsState):
    add_mapping_annotation(model_system.AtomsState.chemical_symbol, XML_KEY, '.label')
    add_mapping_annotation(
        model_system.AtomsState.chemical_symbol, OUTCAR_KEY, '.label'
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
        ('get_eigenvalues', ['.eigenvalues']),
    )
    add_mapping_annotation(
        outputs.Outputs.electronic_eigenvalues,
        XML2_KEY,
        ('get_eigenvalues', ['.eigenvalues']),
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
    add_mapping_annotation(
        outputs.Outputs.electronic_band_gaps,
        XML_KEY,
        ('get_band_gaps', ['.eigenvalues']),
    )
    add_mapping_annotation(
        outputs.Outputs.electronic_band_gaps,
        XML2_KEY,
        ('get_band_gaps', ['.eigenvalues']),
    )
    add_mapping_annotation(
        outputs.Outputs.electronic_band_gaps,
        OUTCAR_KEY,
        ('get_band_gaps', ['.eigenvalues', 'parameters']),
    )
    add_mapping_annotation(
        outputs.Outputs.electronic_dos,
        XML_KEY,
        ('get_total_dos', ['.dos']),
    )
    add_mapping_annotation(
        outputs.Outputs.electronic_dos,
        XML2_KEY,
        ('get_total_dos', ['.dos']),
    )
    add_mapping_annotation(
        outputs.Outputs.electronic_dos,
        DOSCAR_KEY,
        ('get_dos', ['.total_dos', '.projected_dos']),
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

    add_mapping_annotation(outputs.ElectronicEigenvalues.n_levels, XML_KEY, '.n_levels')
    add_mapping_annotation(
        outputs.ElectronicEigenvalues.n_levels, XML2_KEY, '.n_levels'
    )
    add_mapping_annotation(
        outputs.ElectronicEigenvalues.value, XML_KEY, '.value', unit='eV'
    )
    add_mapping_annotation(
        outputs.ElectronicEigenvalues.value, XML2_KEY, '.value', unit='eV'
    )
    add_mapping_annotation(
        outputs.ElectronicEigenvalues.value, OUTCAR_KEY, '.value', unit='eV'
    )
    add_mapping_annotation(
        outputs.ElectronicEigenvalues.occupation, XML_KEY, '.occupation'
    )
    add_mapping_annotation(
        outputs.ElectronicEigenvalues.occupation, XML2_KEY, '.occupation'
    )
    add_mapping_annotation(
        outputs.ElectronicEigenvalues.occupation, OUTCAR_KEY, '.occupation'
    )
    add_mapping_annotation(
        outputs.ElectronicEigenvalues.spin_channel, XML_KEY, '.spin_channel'
    )
    add_mapping_annotation(
        outputs.ElectronicEigenvalues.spin_channel, XML2_KEY, '.spin_channel'
    )
    add_mapping_annotation(
        outputs.ElectronicEigenvalues.spin_channel, OUTCAR_KEY, '.spin_channel'
    )


class ElectronicBandGap(outputs.ElectronicBandGap):
    add_mapping_annotation(
        outputs.ElectronicBandGap.value, XML_KEY, '.value', unit='eV'
    )
    add_mapping_annotation(
        outputs.ElectronicBandGap.value, XML2_KEY, '.value', unit='eV'
    )
    add_mapping_annotation(
        outputs.ElectronicBandGap.value, OUTCAR_KEY, '.value', unit='eV'
    )
    add_mapping_annotation(
        outputs.ElectronicBandGap.spin_channel, XML_KEY, '.spin_channel'
    )
    add_mapping_annotation(
        outputs.ElectronicBandGap.spin_channel, XML2_KEY, '.spin_channel'
    )
    add_mapping_annotation(
        outputs.ElectronicBandGap.spin_channel, OUTCAR_KEY, '.spin_channel'
    )


class Energy2(variables.Energy2):
    add_mapping_annotation(variables.Energy2.points, XML_KEY, '.energies', unit='eV')
    add_mapping_annotation(variables.Energy2.points, XML2_KEY, '.energies', unit='eV')
    add_mapping_annotation(variables.Energy2.points, DOSCAR_KEY, '.energies', unit='eV')


class ElectronicDensityOfStates(outputs.ElectronicDensityOfStates):
    add_mapping_annotation(
        outputs.ElectronicDensityOfStates.value, XML_KEY, '.value', unit='1/eV'
    )
    add_mapping_annotation(
        outputs.ElectronicDensityOfStates.value, XML2_KEY, '.value', unit='1/eV'
    )
    add_mapping_annotation(
        outputs.ElectronicDensityOfStates.value, DOSCAR_KEY, '.value', unit='1/eV'
    )
    add_mapping_annotation(
        outputs.ElectronicDensityOfStates.spin_channel, XML_KEY, '.spin_channel'
    )
    add_mapping_annotation(
        outputs.ElectronicDensityOfStates.spin_channel, XML2_KEY, '.spin_channel'
    )
    add_mapping_annotation(
        outputs.ElectronicDensityOfStates.spin_channel, DOSCAR_KEY, '.spin'
    )
    add_mapping_annotation(
        outputs.ElectronicDensityOfStates.energies_origin,
        XML_KEY,
        '.energy_fermi',
        unit='eV',
    )
    add_mapping_annotation(
        outputs.ElectronicDensityOfStates.energies_origin,
        XML2_KEY,
        '.energy_fermi',
        unit='eV',
    )
    add_mapping_annotation(
        outputs.ElectronicDensityOfStates.energies_origin,
        DOSCAR_KEY,
        'e_fermi',
        unit='eV',
    )
    add_mapping_annotation(variables.Energy2.m_def, XML_KEY, '.@')
    add_mapping_annotation(variables.Energy2.m_def, XML2_KEY, '.@')
    add_mapping_annotation(variables.Energy2.m_def, DOSCAR_KEY, '@')
    add_mapping_annotation(
        outputs.ElectronicDensityOfStates.projected_dos, DOSCAR_KEY, '.projected'
    )


try:
    m_package.__init_metainfo__()
except Exception:
    pass
