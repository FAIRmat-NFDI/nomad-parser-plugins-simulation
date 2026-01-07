#!/usr/bin/env python3

import sys

from nomad.datamodel import EntryArchive

from nomad_simulation_parsers.parsers.h5md.parser import H5MDParser


def debug_bond_parsing():
    parser = H5MDParser()
    archive = EntryArchive()

    print('Starting parse...')
    try:
        parser.parse(
            'tests/data/h5md/test_traj_openmm_5frames_08-08-25.h5', archive, None
        )
        print('Parse completed successfully')

        if hasattr(archive, 'data') and archive.data:
            if hasattr(archive.data, 'model_system') and archive.data.model_system:
                model_system = archive.data.model_system[0]
                print('Root model system found')
                print(f'bond_list present: {model_system.bond_list is not None}')

                if model_system.bond_list is not None:
                    print(f'Number of bonds: {len(model_system.bond_list)}')
                    return True
                else:
                    print('bond_list is None - this is the issue!')
                    return False
            else:
                print('No model_system found in archive.data')
        else:
            print('No data found in archive')

    except Exception as e:
        print(f'Parse failed with error: {e}')
        import traceback

        traceback.print_exc()
        return False

    return False


if __name__ == '__main__':
    success = debug_bond_parsing()
    sys.exit(0 if success else 1)
