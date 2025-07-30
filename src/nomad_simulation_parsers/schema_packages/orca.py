from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass

from nomad.datamodel.metainfo.annotations import Mapper
from nomad.metainfo import SchemaPackage
from nomad.parsing.file_parser.mapping_parser import MAPPING_ANNOTATION_KEY
from nomad_simulations.schema_packages import (
    general,
    atoms_state,
    model_method,
    model_system,
    numerical_settings,
    outputs,
    properties,
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


# parent ModelSystem
general.Simulation.model_system.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(out=Mapper(mapper=('get_atoms', ['.@'])))

# # Basic fields of the parent ModelSystem
# model_system.ModelSystem.positions.m_annotations.setdefault(
#     MAPPING_ANNOTATION_KEY, {}
# ).update(out=Mapper(mapper='.positions'))

# atoms_state.AtomsState.m_def.m_annotations.setdefault(
#     MAPPING_ANNOTATION_KEY, {}
# ).update(out=Mapper(mapper='.particle_states'))

# atoms_state.AtomsState.chemical_symbol.m_annotations.setdefault(
#     MAPPING_ANNOTATION_KEY, {}
# ).update(out=Mapper(mapper='.chemical_symbol'))

# -----------------------------------------------------------------------------
# nested sub‑systems (hypothetical fragments) returned inside `.sub_systems`
# -----------------------------------------------------------------------------

model_system.ModelSystem.sub_systems.m_annotations.setdefault(
    MAPPING_ANNOTATION_KEY, {}
).update(out=Mapper(mapper='.sub_systems'))

# Map fields inside each child ModelSystem
for quantity, path in [
    ('name', '.name'),
    ('n_particles', '.n_particles'),
    ('particle_indices', '.particle_indices'),
]:
    getattr(model_system.ModelSystem, quantity).m_annotations.setdefault(
        MAPPING_ANNOTATION_KEY, {}
    ).update(out=Mapper(mapper=path))




# # # DFT annotations
# model_method.DFT.m_def.m_annotations.setdefault(MAPPING_ANNOTATION_KEY, {}).update(
#     dict(info=Mapper(mapper=('get_dft_data', ['.@'])))
# )

# model_method.DFT.contributions.m_annotations.setdefault(
#     MAPPING_ANNOTATION_KEY, {}
# ).update(dict(info=Mapper(mapper='.@')))

# model_method.DFT.jacobs_ladder.m_annotations.setdefault(
#     MAPPING_ANNOTATION_KEY, {}
# ).update(dict(info=Mapper(mapper='.jacobs_ladder')))

# model_method.DFT.exact_exchange_mixing_factor.m_annotations.setdefault(
#     MAPPING_ANNOTATION_KEY, {}
# ).update(dict(info=Mapper(mapper='.exact_exchange_mixing_factor')))

# model_method.DFT.xc_functionals.m_annotations.setdefault(
#     MAPPING_ANNOTATION_KEY, {}
# ).update(
#     dict(
#         info=Mapper(
#             mapper='.xc_functionals', sub_section=model_method.XCFunctional.m_def
#         )
#     )
# )

# model_method.XCFunctional.libxc_name.m_annotations.setdefault(
#     MAPPING_ANNOTATION_KEY, {}
# ).update(dict(info=Mapper(mapper='.libxc_name')))

# model_method.XCFunctional.name.m_annotations.setdefault(
#     MAPPING_ANNOTATION_KEY, {}
# ).update(dict(info=Mapper(mapper='.name')))

# model_method.XCFunctional.weight.m_annotations.setdefault(
#     MAPPING_ANNOTATION_KEY, {}
# ).update(dict(info=Mapper(mapper='.weight')))



try:
    m_package.__init_metainfo__()
except Exception:
    pass
