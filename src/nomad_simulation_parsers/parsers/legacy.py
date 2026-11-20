from collections.abc import Callable
from importlib import import_module

from nomad.datamodel import EntryArchive
from nomad.parsing import MatchingParser
from nomad.utils import get_logger
from structlog.stdlib import BoundLogger

LOGGER = get_logger(__name__)


class LegacyParserRunner:
    """Run an optional legacy parser alongside the simulation archive writer."""

    def __init__(self, parser_path: str):
        self.parser_path = parser_path

    def run(
        self,
        mainfile: str,
        archive: EntryArchive,
        logger: BoundLogger | None,
        prepare: Callable[[MatchingParser], None] | None = None,
    ) -> MatchingParser | None:
        try:
            module, name = self.parser_path.rsplit('.', 1)
            parser = getattr(import_module(module), name)()
            if prepare is not None:
                prepare(parser)
            parser.parse(mainfile, archive, logger)
            return parser
        except Exception:
            (logger or LOGGER).warning(
                'Legacy parser failed',
                extra={'parser_path': self.parser_path, 'mainfile': mainfile},
                exc_info=True,
            )
            return None
