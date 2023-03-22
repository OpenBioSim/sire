===================================
Molecular Dynamics and Minimisation
===================================

`OpenMM <https://openmm.org>`__ is an excellent package that provides
a framework for running GPU-accelerated molecular dynamics (and related)
simulations.

The :class:`~sire.mol.Molecule`, molecule view and container objects
provide convenience functions that make it easier to use OpenMM
to perform minimisation and molecular dynamics simulations.

Minimisation
------------

You can perform minimisation on a molecule or collection of molecules
using the :func:`~sire.mol.SelectorMol.minimisation` function.
This returns a :class:`~sire.mol.Minimisation` object that can be used
to control minimisation.

>>> import sire as sr
>>> mols = sr.load(sr.expand(sr.tutorial_url, "ala.top", "ala.crd"))
>>> mol = mols[0]
>>> print(mol.energy())
29.3786 kcal mol-1
>>> m = mol.minimisation()

You perform minimisation itself by calling the :func:`~sire.mol.Minimisation.run`
function.

>>> m.run()
minimisation ✔
Minimisation()

You can extract the results of minimisation, converted back into the
original view object by calling :func:`~sire.mol.Minimisation.commit`

>>> mol = m.commit()
>>> print(mol.energy())
17.7799 kcal mol-1

You could run all of these steps on a single line, e.g.

>>> mol = mols[0].minimisation().run().commit()
>>> print(mol.energy())
minimisation ✔
17.7799 kcal mol-1

In the above case we minimised just the first molecule that was loaded.
We can minimise all molecules by calling ``minimisation`` on the whole
collection.

>>> print(mols.energy())
-5855.24 kcal mol-1
>>> mols = mols.minimisation().run().commit()
>>> print(mols.energy())
minimisation ✔
-7971.79 kcal mol-1

Molecular Dynamics
------------------

You can perform molecular dynamics on a molecule or collection of molecules
using the :func:`~sire.mol.SelectorMol.dynamics` function. This returns
a :class:`~sire.mol.Dynamics` object that can be used to control dynamics.

>>> import sire as sr
>>> mols = sr.load(sr.expand(sr.tutorial_url, "ala.top", "ala.crd"))
>>> mol = mols[0]
>>> d = mol.dynamics()
>>> print(d)
Dynamics(completed=0)

You run dynamics by calling the :func:`~sire.mol.Dynamics.run` function.
You pass in the amount of time you want to simulate, and (optionally)
the amount of time between saved coordinate/velocity snapshots.

For example, here we will run 10 picoseconds of dynamics, saving a
frame every 0.5 picoseconds

>>> d.run(10*sr.units.picosecond, 0.5*sr.units.picosecond)
>>> print(d)
Dynamics(completed=10 ps, energy=-8.82721 kcal mol-1, speed=93.7 ns day-1)

.. note::

   The speed of the simulation will depend on whether or not you have a
   GPU, and how fast it is. Reduce the simulation time if you find the
   above example takes too long.

You can extract the results of dynamics by calling
:func:`~sire.mol.Dynamics.commit`.

>>> mol = d.commit()
>>> mol.view()

.. image:: images/05_04_01.jpg
   :alt: Image from the movie of the molecular dynamics trajectory on aladip

In this case, we performed molecular dynamics on just the first molecule
of the loaded system. We can perform dynamics on all the molecules by
calling :func:`~sire.mol.SelectorMol.dynamics` on the complete collection.

>>> d = mols.dynamics()
>>> d.run(10*sr.units.picosecond, 0.5*sr.units.picosecond)
>>> mols = d.commit()
>>> mols.trajectory().energy().pretty_plot()

.. image:: images/05_04_02.jpg
   :alt: Graph of the energies across the trajectory

The frames from dynamics are stored as a trajectory in the molecules.
They can be processed using the
:doc:`trajectory functions introduced previously <../part04/02_trajectory>`.
In this case we called ``pretty_plot`` on the ``energy`` to get a
nice graph of the component energies versus time.

Controlling Dynamics
--------------------

There are several parameters that are needed to control the molecular
dynamics simulation...
