#!/usr/bin/env python3

import sys

from nomad.datamodel import EntryArchive

from nomad_simulation_parsers.parsers.h5md.parser import H5MDParser


def debug_archive_before_normalize():
    parser = H5MDParser()
    archive = EntryArchive()

    print('Starting parse...')
    try:
        parser.parse(
            'tests/data/h5md/test_traj_openmm_5frames_08-08-25.h5', archive, None
        )
        print('Parse completed successfully')

        print('\n=== ARCHIVE STRUCTURE BEFORE NORMALIZATION ===')
        print(f'archive type: {type(archive)}')
        print(f'archive.data: {type(archive.data) if archive.data else None}')
        print(
            f'archive.run: {type(archive.run) if hasattr(archive, "run") and archive.run else None}'
        )
        print(
            f'archive.workflow2: {type(archive.workflow2) if hasattr(archive, "workflow2") and archive.workflow2 else None}'
        )
        print(
            f'archive.metadata: {type(archive.metadata) if hasattr(archive, "metadata") and archive.metadata else None}'
        )

        if archive.data:
            data = archive.data
            print('\n=== archive.data STRUCTURE ===')
            print(f'Type: {type(data)}')
            print(
                f'model_system: {len(data.model_system) if hasattr(data, "model_system") and data.model_system else "None/Empty"}'
            )
            print(
                f'model_method: {len(data.model_method) if hasattr(data, "model_method") and data.model_method else "None/Empty"}'
            )
            print(
                f'computation: {len(data.computation) if hasattr(data, "computation") and data.computation else "None/Empty"}'
            )
            print(
                f'outputs: {len(data.outputs) if hasattr(data, "outputs") and data.outputs else "None/Empty"}'
            )

            if (
                hasattr(data, 'model_system')
                and data.model_system
                and len(data.model_system) > 0
            ):
                print('\n=== ROOT MODEL_SYSTEM ===')
                ms = data.model_system[0]
                print(f'Type: {type(ms)}')
                print(
                    f'positions: {type(ms.positions)} shape={ms.positions.shape if ms.positions is not None and hasattr(ms.positions, "shape") else "None"}'
                )
                print(
                    f'bond_list: {type(ms.bond_list)} len={len(ms.bond_list) if ms.bond_list is not None else "None"}'
                )
                print(
                    f'sub_systems: {len(ms.sub_systems) if hasattr(ms, "sub_systems") and ms.sub_systems else "None/Empty"}'
                )
                print(
                    f'n_particles: {ms.n_particles if hasattr(ms, "n_particles") else "Not set"}'
                )
                print(
                    f'particle_indices: {type(ms.particle_indices)} shape={ms.particle_indices.shape if hasattr(ms, "particle_indices") and ms.particle_indices is not None and hasattr(ms.particle_indices, "shape") else "None"}'
                )

                if hasattr(ms, 'sub_systems') and ms.sub_systems:
                    print('\n=== SUBSYSTEMS ===')
                    for i, subsys in enumerate(ms.sub_systems[:3]):  # Only show first 3
                        print(f'Subsystem {i}:')
                        print(f'  Type: {type(subsys)}')
                        print(
                            f'  name: {subsys.name if hasattr(subsys, "name") else "Not set"}'
                        )
                        print(
                            f'  positions: {type(subsys.positions)} shape={subsys.positions.shape if hasattr(subsys, "positions") and subsys.positions is not None and hasattr(subsys.positions, "shape") else "None"}'
                        )
                        print(
                            f'  bond_list: {type(subsys.bond_list)} len={len(subsys.bond_list) if hasattr(subsys, "bond_list") and subsys.bond_list is not None else "None"}'
                        )
                        print(
                            f'  particle_indices: {type(subsys.particle_indices)} shape={subsys.particle_indices.shape if hasattr(subsys, "particle_indices") and subsys.particle_indices is not None and hasattr(subsys.particle_indices, "shape") else "None"}'
                        )
                        print(
                            f'  n_particles: {subsys.n_particles if hasattr(subsys, "n_particles") else "Not set"}'
                        )
                        if i >= 2:
                            break

        return True

    except Exception as e:
        print(f'Parse failed with error: {e}')
        import traceback

        traceback.print_exc()
        return False


if __name__ == '__main__':
    success = debug_archive_before_normalize()
    sys.exit(0 if success else 1)
