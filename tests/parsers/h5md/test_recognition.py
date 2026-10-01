import h5py
import pytest

from nomad_simulation_parsers.parsers import h5md_parser


@pytest.mark.unit
class TestH5MDRecognition:
    def test_recognizes_h5md_hdf5(self, tmp_path):
        mainfile = tmp_path / 'trajectory.h5'
        with h5py.File(mainfile, 'w') as handle:
            handle.create_group('h5md')

        parser = h5md_parser.load()
        with mainfile.open('rb') as stream:
            contents = stream.read()
        assert (
            parser.is_mainfile(
                str(mainfile), 'application/x-hdf', contents, contents.decode('latin1')
            )
            is True
        )

    @pytest.mark.parametrize('name', ['trajectory.txt', 'trajectory.h5.bak'])
    def test_rejects_non_h5md_filename(self, tmp_path, name):
        mainfile = tmp_path / name
        mainfile.write_bytes(b'not HDF5')
        parser = h5md_parser.load()
        assert not parser.is_mainfile(
            str(mainfile), 'text/plain', b'not HDF5', 'not HDF5'
        )
