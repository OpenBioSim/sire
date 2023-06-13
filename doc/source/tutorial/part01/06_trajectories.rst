===============================
Loading and saving trajectories
===============================

:mod:`sire` also supports loading and saving many popular trajectory
formats.

Simply pass a trajectory file into the list of files and it will be loaded.

>>> import sire as sr
>>> mols = sr.load(sr.expand(sr.tutorial_url, "ala.top", "ala.traj"))

You can find the number of frames that have been loaded using
the :func:`~sire.mol.SelectorMol.num_frames` function.

>>> print(mols.num_frames())
500

You can load a particular frame using the :func:`~sire.mol.SelectorMol.load_frame`
function.

>>> print(mols.time())
0.2 ps
>>> mols.load_frame(100)
>>> print(mols.time())
20.2 ps

.. note::

   The :func:`~sire.mol.SelectorMol.time` function returns the simulation
   time of a particular frame.

Joining multiple trajectories together
--------------------------------------

You can join load multiple trajectory files together and join them
into a single trajectory simply by passing them as additional files
in :func:`sire.load`, e.g.

>>> mols = sr.load(sr.expand(sr.tutorial_url, "ala.top", "ala.traj", "ala.traj"))
>>> print(mols.num_frames())
1000

has loaded ``ala.traj`` twice, joining it together into a single 1000 frame
trajectory.

.. note::

   You can load as many trajectory files as you want. They can be a mixture
   of file formats, e.g. DCD, RST, TRR etc. The trajectory frames will
   be loaded in order, i.e. the frames from the second trajectory file
   in the list will be loaded after the frames in the first file.

.. warning::

   The trajectory frames are streamed from the files. They are not read
   in one go, as typically trajectories are larger than
   available memory. Instead, the trajectory files are opened, with
   frames read dynamically as they are accessed. This means that you
   CANNOT change the trajectory files while they are open (e.g.
   renaming them, moving them, or adding or deleting data). Doing so
   will result in undefined behaviour.

Saving trajectories
-------------------

You can save trajectories using :func:`sire.save`. To do this, you
must pass in the :class:`sire.mol.TrajectoryIterator` that results
from calling the :func:`~sire.mol.SelectorMol.trajectory` function
on a molecule (or set of molecules).

For example,

>>> f = sr.save(mols.trajectory(), "output", format=["RST"])
>>> print(f)
["output.rst"]

saves the trajectory in Amber RST format, returning the
names of the files that were written.

You can save trajectories in any file supported by :mod:`sire`, and
can even save to multiple trajectory formats at the same time, e.g.

>>> f = sr.save(mols.trajectory(), "output", format=["RST", "DCD", "TRR"])
>>> print(f)
["output.rst", "output.dcd", "output.trr"]

saves the trajectory in Amber RST, DCD and Gromacs TRR formats.

You can save only selected frames from the trajectory by slicing the
:class:`sire.mol.TrajectoryIterator`, e.g.

>>> f = sr.save(mols.trajectory()[0:10], "output", format=["XTC"])
>>> print(f)
["output.xtc"]

saves the first ten frames of the trajectory to a Gromacs XTC format
file called ``output.xtc``.

Similarly,

>>> f = sr.save(mols.trajectory()[::100], "output", format=["RST"])
>>> print(f)
["output.rst"]

saves every 100 frames from the trajectory to an Amber RST format file
called ``output.rst``.

You can even save trajectory from only certain molecules from the system.

>>> f = sr.save(mols[0].trajectory(), "output", format=["PRM", "DCD"])
>>> mol = sr.load(f)
>>> print(mol, mol.num_frames())

has saved the trajectory of only the first molecule to a DCD file called
``output.dcd``, and also saved an Amber PRM file for that first molecule
to ``output.prm``. These files were immediately reloaded, showing
you that they contain just the first molecule, with 1000 frames of
trajectory.

Trajectories as frames
----------------------

:mod:`sire` also supports loading and saving trajectories as a sequence
of individual files. For example, here we will save the trajectory
of the first molecule as a sequence of PDB files.

>>> f = sr.save(mols[0].trajectory(), "output", format=["PRM", "PDB"])
>>> print(f)
["output.prm", "output.pdb"]

In this case, ``output.pdb`` is a directory containing 1000 PDB files, one
for each trajectory frame.

>>> import os
>>> os.listdir("output.pdb")
['frame_143_28-8.pdb',
 'frame_304_61.pdb',
 'frame_006_1-4.pdb',
 'frame_301_60-4.pdb',
 'frame_220_44-2.pdb',
...
 'frame_421_84-4.pdb',
 'frame_196_39-4.pdb',
 'frame_056_11-4.pdb',
 'frame_428_85-8.pdb']

We can load these files as a trajectory simply by passing in the directory
name, e.g.

>>> mol = sr.load("output.prm", "output.pdb")
>>> print(mol, mol.num_frames())
System( name=output num_molecules=1 num_residues=3 num_atoms=22 ) 1000

.. note::

   The individual frame files are named using the frame number plus
   the time in picoseconds for that frame.

:mod:`sire` has extensive support for trajectories. To learn more,
check out the :doc:`detailed guide <../../cheatsheet/trajectory>`.
