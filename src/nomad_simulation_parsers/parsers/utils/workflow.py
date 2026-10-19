from collections.abc import Iterable

from nomad.datamodel import EntryArchive
from nomad.datamodel.metainfo.workflow import Link
from nomad_simulations.schema_packages.workflow import Phonon
from nomad_simulations.schema_packages.workflow.general import (
    SimulationTask,
    SimulationTaskReference,
)
from nomad_simulations.schema_packages.workflow.phonon import (
    PhononMethod,
    PhononResults,
)


def create_phonon_workflow(
    phonon_archive: EntryArchive,
    force_archives: Iterable[EntryArchive] = (),
) -> Phonon:
    """Build one workflow from a Phonopy result and its known force calculations."""
    simulation = phonon_archive.data
    unit_cell = simulation.model_system[0]
    output = simulation.outputs[0]
    method = simulation.model_method[0] if simulation.model_method else None
    input_link = Link(name='Input system', section=unit_cell)

    force_tasks = [
        SimulationTaskReference(task=archive.workflow2)
        for archive in force_archives
        if archive.workflow2 is not None
    ]
    if not force_tasks:
        force_tasks = [SimulationTask(name='Force calculation')]
    force_tasks[0].inputs.append(Link(name='Input system', section=unit_cell))

    phonopy_task = SimulationTask(
        name='Phonopy analysis',
        outputs=[Link(name='Phonon properties', section=output)],
    )
    return Phonon(
        inputs=[input_link],
        tasks=[*force_tasks, phonopy_task],
        method=PhononMethod(
            force_calculator=method.force_calculator if method else None,
            phonon_calculator='Phonopy',
            mesh_density=method.mesh_density if method else None,
            with_non_analytic_correction=(
                method.with_non_analytic_correction if method else None
            ),
        ),
        results=PhononResults(n_imaginary_frequencies=output.n_imaginary_frequencies),
    )
