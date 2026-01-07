#!/usr/bin/env python3

import signal
import sys

from nomad.datamodel import EntryArchive

from nomad_simulation_parsers.parsers.h5md.parser import H5MDParser


def timeout_handler(signum, frame):
    print('TIMEOUT: Parser hanging - aborting')
    sys.exit(1)


def minimal_test():
    # Set a 20 second timeout
    signal.signal(signal.SIGALRM, timeout_handler)
    signal.alarm(20)

    try:
        print('Creating parser...')
        parser = H5MDParser()
        print('Creating archive...')
        archive = EntryArchive()

        print('Starting parse... (timeout in 20s)')
        parser.parse(
            'tests/data/h5md/test_traj_openmm_5frames_08-08-25.h5', archive, None
        )

        signal.alarm(0)  # Cancel timeout
        print('Parse completed!')

        print(f'Archive has data: {archive.data is not None}')
        if archive.data:
            print(f'Data type: {type(archive.data)}')
            if hasattr(archive.data, 'model_system'):
                print(
                    f'Model systems: {len(archive.data.model_system) if archive.data.model_system else 0}'
                )

        return True

    except Exception as e:
        signal.alarm(0)  # Cancel timeout
        print(f'Parse failed: {e}')
        return False


if __name__ == '__main__':
    success = minimal_test()
    sys.exit(0 if success else 1)
