=======================
Viewing Molecules in 2D
=======================

The ability to convert :mod:`sire` objects to native objects of other
molecular packages makes it easy to create convenience functions that
leverage the capabilities of those packages under the hood.

For example, RDKit can be used to create 2D views of molecules.
This functionality is used by the :func:`~sire.mol.view2d` function to
create a two-dimensional view of any :class:`~sire.mol.Molecule`,
molecule sub-view or collections of molecules.

In a jupyter notebook you can run;

>>> import sire as sr
>>> mols = sr.load(sr.expand(sr.tutorial_url, "ala.crd", "ala.top"))
>>> mols[0].view2d()

