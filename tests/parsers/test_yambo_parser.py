from nomad.datamodel import EntryArchive
from nomad.utils import get_logger

from nomad_simulation_parsers.parsers.yambo  import YamboArchiveWriter   

LOGGER = get_logger(__name__)

def test_spectra(parser):
  archive = EntryArchive()
  parser = YamboArchiveWriter()
  parser.parse('tests/data/yambo/r_setup',archive,LOGGER)
  assert spectra_files[0] == 'o-R_methylox_TDLDA.alpha_q1_slepc_alda_bse'
  spectra = archive.data.outputs[0]
  assert spectra.sp_type == 'Polarizability'
  assert spectra_obj.shape == (4000,2)
