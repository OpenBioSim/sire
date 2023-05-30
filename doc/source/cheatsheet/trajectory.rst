============
Trajectories
============

:mod:`sire` supports loading, editing and saving of trajectories via multiple
file formats. Trajectories are represented as a series of frames, with
each frame containing (optionally) coordinate, velocity, force, space,
time and generic property data.

Supported file formats
----------------------

Table of formats


Loading trajectories
--------------------

Loading trajectories. Debugging via direct use of sire.io.parser objects.

Loading multiple frames via mulitple files. Joining multiple files together.

Saving trajectories
-------------------

Saving trajectories. Difference between sr.save(mols) and
sr.save(mols.trajectory())

Skipping frames. Specifying frames. Aligning, centering or smoothing frames.
Saving to muliple formats at once.

Skipping molecules - saving trajectories only for certain views etc

Saving to a directory of individual frame files (e.g. a directory of PDBs)

Visualising trajectories
------------------------

Link to view. Mention view automatically, plus control skipping,
subsetting, aligning, centering and smoothing.

Skipping molecules - viewing trajectories only for certain views etc

Machinery
---------

Go into a little more detail about the `Trajectory` object (and associated
property) and also the `Frame` object and `TrajectoryIterator`.


