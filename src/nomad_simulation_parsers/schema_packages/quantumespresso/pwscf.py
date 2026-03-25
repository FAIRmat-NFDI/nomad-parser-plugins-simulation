from nomad.metainfo import SchemaPackage
from nomad_simulations.schema_packages import general, outputs, variables

from nomad_simulation_parsers.schema_packages.utils import add_mapping_annotation

from .common import OUT_KEY, XML_KEY

m_package = SchemaPackage()


class TotalForce(outputs.TotalForce):
    add_mapping_annotation(
        outputs.TotalForce.value, OUT_KEY, '.value || .forces', unit='rydberg/bohr'
    )
    add_mapping_annotation(
        outputs.TotalForce.value, XML_KEY, ('get_forces', ['.__value'])
    )
    add_mapping_annotation(
        outputs.TotalForce.contributions,
        OUT_KEY,
        ('get_force_contributions', ['.@']),
        unit='rydberg/bohr',
    )


class ElectronicEigenvalues(outputs.ElectronicEigenvalues):
    add_mapping_annotation(outputs.ElectronicEigenvalues.n_levels, OUT_KEY, '.n_levels')
    add_mapping_annotation(outputs.ElectronicEigenvalues.value, OUT_KEY, '.eigenvalues')
    add_mapping_annotation(
        outputs.ElectronicEigenvalues.occupation, OUT_KEY, '.occupations'
    )


class ElectronicBandStructure(outputs.ElectronicBandStructure):
    add_mapping_annotation(outputs.ElectronicBandStructure.value, OUT_KEY, '.value')
    add_mapping_annotation(
        outputs.ElectronicBandStructure.occupation, OUT_KEY, '.occupation'
    )
    add_mapping_annotation(
        outputs.ElectronicBandStructure.n_levels, OUT_KEY, '.n_levels'
    )
    add_mapping_annotation(
        outputs.ElectronicBandStructure.spin_channel, OUT_KEY, '.spin_channel'
    )
    add_mapping_annotation(
        outputs.ElectronicBandStructure.highest_occupied,
        OUT_KEY,
        '.highest_occupied',
    )


class Energy2(variables.Energy2):
    add_mapping_annotation(variables.Energy2.points, OUT_KEY, '.points')


class ElectronicDensityOfStates(outputs.ElectronicDensityOfStates):
    add_mapping_annotation(outputs.ElectronicDensityOfStates.value, OUT_KEY, '.value')
    add_mapping_annotation(
        outputs.ElectronicDensityOfStates.energies_origin,
        OUT_KEY,
        '.energies_origin',
    )
    add_mapping_annotation(variables.Energy2.m_def, OUT_KEY, '.energies')


class Outputs(outputs.Outputs):
    add_mapping_annotation(outputs.Outputs.total_forces, OUT_KEY, '.@')
    add_mapping_annotation(outputs.Outputs.total_forces, XML_KEY, '.forces')
    add_mapping_annotation(
        outputs.Outputs.electronic_eigenvalues,
        OUT_KEY,
        '.electronic_eigenvalues',
    )
    add_mapping_annotation(
        outputs.Outputs.scf_steps,
        OUT_KEY,
        ('get_scf_steps', ['.@']),
    )
    add_mapping_annotation(
        outputs.Outputs.scf_steps,
        XML_KEY,
        ('get_scf_steps', ['.@']),
    )
    add_mapping_annotation(
        outputs.Outputs.electronic_band_structures,
        OUT_KEY,
        '.electronic_band_structures',
    )
    add_mapping_annotation(
        outputs.Outputs.electronic_dos,
        OUT_KEY,
        '.electronic_dos',
    )


class SCFSteps(outputs.SCFSteps):
    add_mapping_annotation(outputs.SCFSteps.energies_total, OUT_KEY, '.energies_total')
    add_mapping_annotation(outputs.SCFSteps.energies_total, XML_KEY, '.energies_total')
    add_mapping_annotation(
        outputs.SCFSteps.delta_energies_total, OUT_KEY, '.delta_energies_total'
    )
    add_mapping_annotation(
        outputs.SCFSteps.delta_energies_total, XML_KEY, '.delta_energies_total'
    )
    add_mapping_annotation(outputs.SCFSteps.durations, OUT_KEY, '.durations')
    add_mapping_annotation(
        outputs.SCFSteps.code_specific_quantities, OUT_KEY, '.code_specific_quantities'
    )
    add_mapping_annotation(
        outputs.SCFSteps.code_specific_quantities, XML_KEY, '.code_specific_quantities'
    )
    add_mapping_annotation(
        outputs.Outputs.scf_steps,
        OUT_KEY,
        ('get_scf_steps', ['.@']),
    )
    add_mapping_annotation(
        outputs.Outputs.scf_steps,
        XML_KEY,
        ('get_scf_steps', ['.@']),
    )


class SCFSteps(outputs.SCFSteps):
    add_mapping_annotation(outputs.SCFSteps.energies_total, OUT_KEY, '.energies_total')
    add_mapping_annotation(outputs.SCFSteps.energies_total, XML_KEY, '.energies_total')
    add_mapping_annotation(
        outputs.SCFSteps.delta_energies_total, OUT_KEY, '.delta_energies_total'
    )
    add_mapping_annotation(
        outputs.SCFSteps.delta_energies_total, XML_KEY, '.delta_energies_total'
    )
    add_mapping_annotation(outputs.SCFSteps.durations, OUT_KEY, '.durations')
    add_mapping_annotation(
        outputs.SCFSteps.code_specific_quantities, OUT_KEY, '.code_specific_quantities'
    )
    add_mapping_annotation(
        outputs.SCFSteps.code_specific_quantities, XML_KEY, '.code_specific_quantities'
    )


class Simulation(general.Simulation):
    add_mapping_annotation(
        general.Simulation.model_system,
        OUT_KEY,
        ('get_configurations', ['.@']),
        cache=True,
    )
    add_mapping_annotation(
        general.Simulation.model_system,
        XML_KEY,
        ('get_configurations', ['.@']),
        cache=True,
    )
    add_mapping_annotation(
        general.Simulation.outputs, OUT_KEY, ('get_configurations', ['.@']), cache=True
    )
    add_mapping_annotation(
        general.Simulation.outputs, XML_KEY, ('get_configurations', ['.@']), cache=True
    )


try:
    m_package.__init_metainfo__()
except Exception:
    pass
