import numpy as np
import panedr
from nomad_file_parser import FileParser


class GromacsEDRParser(FileParser):
    @property
    def fileedr(self):
        if self._file_handler is None:
            try:
                self._file_handler = panedr.edr_to_df(self.mainfile)
            except Exception:
                self.logger.error('Error reading edr file.')

        return self._file_handler

    def parse(self, key: str):
        if self.fileedr is None:
            return

        val = self.fileedr.get(key, None)
        if self._results is None:
            self._results = dict()

        if val is not None:
            val = np.asarray(val)

        self._results[key] = val

    def keys(self) -> list[str]:
        if self.fileedr is None:
            return []
        return list(self.fileedr.keys())

    @property
    def length(self) -> int:
        return self.fileedr.shape[0]
