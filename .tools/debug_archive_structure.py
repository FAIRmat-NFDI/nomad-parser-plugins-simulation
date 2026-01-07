#!/usr/bin/env python3


from nomad.datamodel import EntryArchive

from nomad_simulation_parsers.parsers.h5md.parser import H5MDParser


def debug_archive_structure():
    parser = H5MDParser()
    archive = EntryArchive()

    print('Starting parse...')
    try:
        parser.parse(
            'tests/data/h5md/test_traj_openmm_5frames_08-08-25.h5', archive, None
        )
        print('Parse completed successfully')

        print('Archive structure:')
        print(f'  archive.data: {type(archive.data) if archive.data else None}')
        print(
            f'  archive.run: {type(archive.run) if hasattr(archive, "run") and archive.run else None}'
        )
        print(
            f'  archive.workflow2: {type(archive.workflow2) if hasattr(archive, "workflow2") and archive.workflow2 else None}'
        )
        print(
            f'  archive.metadata: {type(archive.metadata) if hasattr(archive, "metadata") and archive.metadata else None}'
        )

        if archive.data:
            data = archive.data
            print('\narchive.data structure:')
            print(
                f'  model_system: {len(data.model_system) if hasattr(data, "model_system") and data.model_system else "None"}'
            )
            print(
                f'  model_method: {len(data.model_method) if hasattr(data, "model_method") and data.model_method else "None"}'
            )
            print(
                f'  computation: {len(data.computation) if hasattr(data, "computation") and data.computation else "None"}'
            )

            if (
                hasattr(data, 'model_system')
                and data.model_system
                and len(data.model_system) > 0
            ):
                ms = data.model_system[0]
                print('\nFirst model_system:')
                print(
                    f'  positions: {"Present" if ms.positions is not None else "None"}'
                )
                print(
                    f'  bond_list: {"Present" if ms.bond_list is not None else "None"}'
                )
                print(
                    f'  sub_systems: {len(ms.sub_systems) if hasattr(ms, "sub_systems") and ms.sub_systems else "None"}'
                )

    except Exception as e:
        print(f'Parse failed with error: {e}')
        import traceback

        traceback.print_exc()
        return False

    return True


if __name__ == '__main__':
    debug_archive_structure()
