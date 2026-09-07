"""Shared assertions for parser tests.

See FAIRmat-NFDI/nomad-simulations#474 for the storage convention these enforce.
"""

from nomad.datamodel import EntryArchive


def assert_identity_populated_once(
    archive: EntryArchive, topology_index: int = 0
) -> None:
    """Assert per-particle identity is stored on exactly one (topology) frame.

    In a multi-frame ``model_system`` sequence -- an MD trajectory or a
    geometry-optimization step series -- the per-particle identity
    (``particle_states``) is frame-independent and must be written only on the
    topology/representative frame, never duplicated per frame. Duplicating it
    scales the archive and the Elasticsearch index document with
    ``n_frames x n_particles`` and can push a single entry past the Elasticsearch
    payload limit, failing the whole upload (FAIRmat-NFDI/nomad-simulations#474).

    Meaningful only for a multi-frame fixture; a single-frame system passes
    trivially.
    """
    systems = archive.data.model_system
    populated = [
        i for i, system in enumerate(systems) if len(system.particle_states) > 0
    ]
    assert len(populated) == 1, (
        f'`particle_states` populated on {len(populated)}/{len(systems)} '
        f'`model_system` frames; expected exactly one (topology) frame. '
        f'Offending frame indices: {populated[:10]}'
    )
    assert populated[0] == topology_index, (
        f'`particle_states` found on frame {populated[0]}, '
        f'expected the topology frame {topology_index}'
    )
