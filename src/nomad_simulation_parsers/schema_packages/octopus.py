from nomad.metainfo import SchemaPackage
from nomad_simulations.schema_packages import (
    general,
    model_method,
    model_system,
    outputs,
)

from nomad_simulation_parsers.schema_packages.utils import add_mapping_annotation

m_package = SchemaPackage()

OUT_KEY = 'octopus_out'
INFO_KEY = 'octopus_info'
EIGENVALUES_KEY = 'octopus_eigenvalues'


class Program(general.Program):
    add_mapping_annotation(general.Program.version, OUT_KEY, '.Version')


class XCComponent(model_method.XCComponent):
    add_mapping_annotation(model_method.XCComponent.canonical_label, OUT_KEY, '.@')


class XCFunctional(model_method.XCFunctional):
    add_mapping_annotation(
        model_method.XCFunctional.components,
        OUT_KEY,
        ('get_xc_functionals', ['.theory_level']),
    )


class DFT(model_method.DFT):
    add_mapping_annotation(model_method.DFT.xc, OUT_KEY, '.@')


class ModelSystem(model_system.ModelSystem):
    add_mapping_annotation(model_system.ModelSystem.positions, OUT_KEY, '.positions')
    add_mapping_annotation(model_system.Representation.m_def, OUT_KEY, '.@')
    add_mapping_annotation(model_system.AtomsState.m_def, OUT_KEY, '.labels')


class Representation(model_system.Representation):
    add_mapping_annotation(
        model_system.Representation.lattice_vectors, OUT_KEY, '.lattice_vectors'
    )
    add_mapping_annotation(
        model_system.Representation.periodic_boundary_conditions, OUT_KEY, '.pbc'
    )


class AtomsState(model_system.AtomsState):
    add_mapping_annotation(model_system.AtomsState.chemical_symbol, OUT_KEY, '.@')


class TotalEnergy(outputs.TotalEnergy):
    add_mapping_annotation(outputs.TotalEnergy.value, OUT_KEY, '.energy')
    add_mapping_annotation(outputs.TotalEnergy.value, INFO_KEY, '.Total || .value')


class TotalForce(outputs.TotalForce):
    add_mapping_annotation(outputs.TotalForce.value, INFO_KEY, '.forces')


class ElectronicEigenvalues(outputs.ElectronicEigenvalues):
    add_mapping_annotation(
        outputs.ElectronicEigenvalues.value, INFO_KEY, '.eigenvalues'
    )
    add_mapping_annotation(
        outputs.ElectronicEigenvalues.value, EIGENVALUES_KEY, '.eigenvalues'
    )
    add_mapping_annotation(
        outputs.ElectronicEigenvalues.occupation, INFO_KEY, '.occupations'
    )
    add_mapping_annotation(
        outputs.ElectronicEigenvalues.occupation, EIGENVALUES_KEY, '.occupations'
    )
    add_mapping_annotation(
        outputs.ElectronicEigenvalues.highest_occupied,
        INFO_KEY,
        '.highest_occupied',
    )
    add_mapping_annotation(
        outputs.ElectronicEigenvalues.highest_occupied,
        EIGENVALUES_KEY,
        '.highest_occupied',
    )


class ElectronicBandStructure(outputs.ElectronicBandStructure):
    add_mapping_annotation(outputs.ElectronicBandStructure.value, INFO_KEY, '.value')
    add_mapping_annotation(
        outputs.ElectronicBandStructure.value, EIGENVALUES_KEY, '.value'
    )
    add_mapping_annotation(
        outputs.ElectronicBandStructure.highest_occupied,
        INFO_KEY,
        '.highest_occupied',
    )
    add_mapping_annotation(
        outputs.ElectronicBandStructure.highest_occupied,
        EIGENVALUES_KEY,
        '.highest_occupied',
    )


class ElectronicBandGap(outputs.ElectronicBandGap):
    add_mapping_annotation(outputs.ElectronicBandGap.value, INFO_KEY, '.value')
    add_mapping_annotation(outputs.ElectronicBandGap.value, EIGENVALUES_KEY, '.value')
    add_mapping_annotation(
        outputs.ElectronicBandGap.spin_channel, INFO_KEY, '.spin_channel'
    )
    add_mapping_annotation(
        outputs.ElectronicBandGap.spin_channel, EIGENVALUES_KEY, '.spin_channel'
    )


class Outputs(outputs.Outputs):
    add_mapping_annotation(outputs.Outputs.total_energies, OUT_KEY, '.@')
    add_mapping_annotation(outputs.Outputs.total_energies, INFO_KEY, '.energies')
    add_mapping_annotation(outputs.Outputs.total_forces, INFO_KEY, '.@')
    add_mapping_annotation(
        outputs.Outputs.electronic_eigenvalues,
        INFO_KEY,
        ('get_eigenvalues', ['eigenvalues.eigenvalues']),
    )
    add_mapping_annotation(
        outputs.Outputs.electronic_eigenvalues,
        EIGENVALUES_KEY,
        ('get_eigenvalues', ['eigenvalues']),
    )
    add_mapping_annotation(
        outputs.Outputs.electronic_band_structures,
        INFO_KEY,
        ('get_band_structures', ['eigenvalues.eigenvalues']),
    )
    add_mapping_annotation(
        outputs.Outputs.electronic_band_structures,
        EIGENVALUES_KEY,
        ('get_band_structures', ['eigenvalues']),
    )
    add_mapping_annotation(
        outputs.Outputs.electronic_band_gaps,
        INFO_KEY,
        ('get_band_gaps', ['eigenvalues.eigenvalues']),
    )
    add_mapping_annotation(
        outputs.Outputs.electronic_band_gaps,
        EIGENVALUES_KEY,
        ('get_band_gaps', ['eigenvalues']),
    )
    # TODO(legacy-parity): legacy Octopus parser did not populate explicit
    # electronic DOS result sections.


class Simulation(general.Simulation):
    add_mapping_annotation(
        general.Simulation.program, OUT_KEY, ('get_header', ['.header'])
    )
    add_mapping_annotation(model_method.DFT.m_def, OUT_KEY, '.@')
    add_mapping_annotation(
        general.Simulation.model_system, OUT_KEY, ('get_systems', ['.minimization'])
    )
    add_mapping_annotation(
        general.Simulation.outputs,
        OUT_KEY,
        ('get_outputs', ['.minimization || .time_dependent[].iteration']),
    )
    add_mapping_annotation(general.Simulation.outputs, INFO_KEY, '.@')
    add_mapping_annotation(general.Simulation.outputs, EIGENVALUES_KEY, '.@')


add_mapping_annotation(general.Simulation.m_def, OUT_KEY, '@')
add_mapping_annotation(general.Simulation.m_def, INFO_KEY, '@')
add_mapping_annotation(general.Simulation.m_def, EIGENVALUES_KEY, '@')


try:
    m_package.__init_metainfo__()
except Exception:
    pass
