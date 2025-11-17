from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass

from nomad.datamodel.metainfo.annotations import Mapper
from nomad.metainfo import SchemaPackage
from nomad.parsing.file_parser.mapping_parser import MAPPING_ANNOTATION_KEY
from nomad_simulations.schema_packages import (
    atoms_state,
    general,
    model_system,
    numerical_settings,
)

m_package = SchemaPackage()

# Simulation
general.Simulation.m_def.m_annotations.setdefault(MAPPING_ANNOTATION_KEY, {}).update(
    dict(out=Mapper(mapper='@'))
)

# Program
general.Simulation.program.m_annotations.setdefault(MAPPING_ANNOTATION_KEY, {}).update(
    dict(out=Mapper(mapper='.@'))
)

# Program quantities
general.Program.version.m_annotations.setdefault(MAPPING_ANNOTATION_KEY, {}).update(
    dict(out=Mapper(mapper='.program_version'))
)

general.Simulation.model_system.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(dict(out=Mapper(mapper=('get_atoms', ['.@']))))


model_system.ModelSystem.positions.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(dict(out=Mapper(mapper='.positions')))


atoms_state.AtomsState.m_def.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(dict(out=Mapper(mapper='.particle_states')))

atoms_state.AtomsState.chemical_symbol.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(dict(out=Mapper(mapper='.chemical_symbol')))


############# DFT ###################


############# numerical settings ###################

numerical_settings.SelfConsistency.m_def.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(dict(out=Mapper(mapper=('get_numerical_settings', ['.@']))))

# individual quantities inside that section
numerical_settings.SelfConsistency.n_max_iterations.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(dict(out=Mapper(mapper='.n_max_iterations')))

numerical_settings.SelfConsistency.threshold_change.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(dict(out=Mapper(mapper='.threshold_change')))


try:
    m_package.__init_metainfo__()
except Exception:
    pass
