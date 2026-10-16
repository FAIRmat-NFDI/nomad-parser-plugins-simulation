# GROMACS hierarchy disambiguation fixture

This is a minimal reproducer derived from the hierarchy in the Frame49 issue
data. It contains six distinct protein-chain molecule types. Each chain has two
`sA` residues, so equivalent residue groups occur below different parents.

The fixture was generated with GROMACS 2025.2 from `frame49-mini.top`,
`frame49-mini.gro`, and `frame49-mini.mdp`. The one-step trajectory has 12
particles and no interactions. It replaces the 23 MB Frame49 archive while
preserving the parent-scoped group-label behavior under test.
