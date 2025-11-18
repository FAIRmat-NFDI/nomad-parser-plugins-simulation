from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass
from nomad.metainfo import SchemaPackage

m_package = SchemaPackage()
try:
    m_package.__init_metainfo__()
except Exception:
    pass
