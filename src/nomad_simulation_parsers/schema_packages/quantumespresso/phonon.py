from nomad.metainfo import SchemaPackage
from nomad_simulations.schema_packages import general

from nomad_simulation_parsers.schema_packages.utils import add_mapping_annotations

m_package = SchemaPackage()


class Simulation(general.Simulation):
    add_mapping_annotations(general.Simulation.model_system, 'out', '.calculation')


try:
    m_package.__init_metainfo__()
except Exception:
    pass
