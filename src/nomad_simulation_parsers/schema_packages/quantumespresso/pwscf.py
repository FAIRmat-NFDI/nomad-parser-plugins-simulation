from nomad.metainfo import SchemaPackage
from nomad_simulations.schema_packages import general, outputs

from nomad_simulation_parsers.schema_packages.utils import add_mapping_annotations

from .common import OUT_KEY, XML_KEY

m_package = SchemaPackage()


class TotalForce(outputs.TotalForce):
    add_mapping_annotations(
        outputs.TotalForce.value, OUT_KEY, '.value || .forces', unit='rydberg/bohr'
    )
    add_mapping_annotations(
        outputs.TotalForce.value, XML_KEY, ('get_forces', ['.__value'])
    )
    add_mapping_annotations(
        outputs.TotalForce.contributions,
        OUT_KEY,
        ('get_force_contributions', ['.@']),
        unit='rydberg/bohr',
    )


class ElectronicEigenvalues(outputs.ElectronicEigenvalues):
    add_mapping_annotations(
        outputs.ElectronicEigenvalues.value, OUT_KEY, '.eigenvalues'
    )


class Outputs(outputs.Outputs):
    add_mapping_annotations(outputs.Outputs.total_forces, OUT_KEY, '.@')
    add_mapping_annotations(outputs.Outputs.total_forces, XML_KEY, '.forces')
    add_mapping_annotations(
        outputs.Outputs.electronic_eigenvalues, OUT_KEY, ('get_eigenvalues', ['.@'])
    )


class Simulation(general.Simulation):
    add_mapping_annotations(
        general.Simulation.model_system,
        OUT_KEY,
        ('get_configurations', ['.@']),
        cache=True,
    )
    add_mapping_annotations(
        general.Simulation.model_system,
        XML_KEY,
        ('get_configurations', ['.@']),
        cache=True,
    )
    add_mapping_annotations(
        general.Simulation.outputs, OUT_KEY, ('get_configurations', ['.@']), cache=True
    )
    add_mapping_annotations(
        general.Simulation.outputs, XML_KEY, ('get_configurations', ['.@']), cache=True
    )


try:
    m_package.__init_metainfo__()
except Exception:
    pass
