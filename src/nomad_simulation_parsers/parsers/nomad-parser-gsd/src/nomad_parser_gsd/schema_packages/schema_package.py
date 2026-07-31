from typing import (
    TYPE_CHECKING,
)

if TYPE_CHECKING:
    from nomad.datamodel.datamodel import (
        EntryArchive,
    )
    from structlog.stdlib import (
        BoundLogger,
    )

import numpy as np
from nomad.config import config
from nomad.datamodel.data import ArchiveSection, Schema
from nomad.datamodel.metainfo.annotations import ELNAnnotation, ELNComponentEnum
from nomad.metainfo import (
    Context,
    MEnum,
    Quantity,
    SchemaPackage,
    Section,
    SectionProxy,
    SubSection,
)

configuration = config.get_plugin_entry_point(
    'nomad_parser_gsd.schema_packages:schema_package_entry_point'
)

m_package = SchemaPackage()


class GSDSchemaPackage(Schema):
    name = Quantity(
        type=str, a_eln=ELNAnnotation(component=ELNComponentEnum.StringEditQuantity)
    )
    message = Quantity(type=str)

    def normalize(self, archive: 'EntryArchive', logger: 'BoundLogger') -> None:
        super().normalize(archive, logger)

        logger.info('GSDSchema.normalize', parameter=configuration.parameter)
        self.message = f'Hello {self.name}!'


m_package.__init_metainfo__()
