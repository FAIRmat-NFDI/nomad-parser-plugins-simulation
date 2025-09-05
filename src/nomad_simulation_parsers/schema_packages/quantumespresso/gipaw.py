from nomad.metainfo import SchemaPackage
from nomad_simulations.schema_packages import general, outputs

m_package = SchemaPackage()


class Outputs(outputs.Outputs):
    pass


class Simulation(general.Simulation):
    pass


try:
    m_package.__init_metainfo__()
except Exception:
    pass
