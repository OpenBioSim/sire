=================
Merging Molecules
=================

This page is a WORK IN PROGRESS.

This page will talk about how we can create merge molecules.

This uses BioSimSpace, so examples will start with

>>> import BioSimSpace as BSS
>>> import sire as sr

This page will then detail additional functionality that can be used
to control or edit the merge molecule.

Examples include shrinking ghost atoms.

>>> sr.morph.shrink_ghost_atoms(mols, length="0.6 A")

Or doing extra work like adding restraints, or adding extra controls
or forcefield parameters that represent best practice for complex
morphs.

Also show how to visualise the morph, e.g.

>>> merged_mol.perturbation().view()

Also show the interface for checking if a molecule is perturbable

>>> merged_mol.is_perturbable()

And to search for perturbable molecules in a collection.

>>> perturbable_mols = mols.molecules("property is_perturbable")

or

>>> mol = mols["molecule property is_perturbable"]

Also show how to use links to link to the 0 or 1 properties

>>> mol = mol.edit().add_link("coordinates0", "coordinates").commit()

(would be good to add this to the Cursor API)
