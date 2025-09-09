from nomad.metainfo import SchemaPackage
from nomad_simulations.schema_packages import general, outputs

from nomad_simulation_parsers.schema_packages.utils import add_mapping_annotations

m_package = SchemaPackage()


class TotalForce(outputs.TotalForce):
    add_mapping_annotations(
        outputs.TotalForce.value, 'out', '.value || .forces', unit='rydberg/bohr'
    )
    add_mapping_annotations(
        outputs.TotalForce.contributions,
        'out',
        ('get_force_contributions', ['.@']),
        unit='rydberg/bohr',
    )


class ElectronicEigenvalues(outputs.ElectronicEigenvalues):
    add_mapping_annotations(outputs.ElectronicEigenvalues.value, 'out', '.eigenvalues')


class Outputs(outputs.Outputs):
    add_mapping_annotations(outputs.Outputs.total_forces, 'out', '.@')
    add_mapping_annotations(
        outputs.Outputs.electronic_eigenvalues, 'out', ('get_eigenvalues', ['.@'])
    )


class Simulation(general.Simulation):
    add_mapping_annotations(
        general.Simulation.model_system,
        'out',
        ('get_configurations', ['.@']),
        cache=True,
    )
    add_mapping_annotations(
        general.Simulation.outputs, 'out', ('get_configurations', ['.@']), cache=True
    )


try:
    m_package.__init_metainfo__()
except Exception:
    pass
