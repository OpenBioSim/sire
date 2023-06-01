============
Trajectories
============

:mod:`sire` supports loading, editing and saving of trajectories via multiple
file formats. Trajectories are represented as a series of frames, with
each frame containing (optionally) coordinate, velocity, force, space,
time and generic property data.

Supported file formats
----------------------

:mod:`sire` natively supports the following trajectory formats:

+----------+----------------------------+------------------------------------+
| Format   | Parser                     | Description                        |
+==========+============================+====================================+
| DCD      | :cls:`sire.io.parser.DCD`  | DCD coordinate/velocity binary     |
|          |                            | trajectory files based             |
|          |                            | on charmm / namd / x-plor format.  |
+----------+----------------------------+------------------------------------+
| RST      | :cls:`sire.io.parser.RST`  | Amber coordinate/velocity binary   |
|          |                            | (netcdf) restart/trajectory files  |
|          |                            | supported since Amber 9, now       |
|          |                            | default since Amber 16.            |
+----------+----------------------------+------------------------------------+
| TRAJ     | :cls:`sire.io.parser.TRAJ` | Amber trajectory (ascii)           |
|          |                            | coordinate or velocity files       |
|          |                            | supported from Amber 7 upwards.    |
+----------+----------------------------+------------------------------------+
| TRR      | :cls:`sire.io.parser.TRR`  | Gromacs TRR (XDR file) coordinate  |
|          |                            | / velocity / force trajectory      |
|          |                            | file                               |
+----------+----------------------------+------------------------------------+
| XTC      | :cls:`sire.io.parser.XTC`  | Gromacs XTC (XDR file) compressed  |
|          |                            | coordinate trajectory file         |
+----------+----------------------------+------------------------------------+

The :cls:`~sire.io.parser.DCD`, :cls:`~sire.io.parser.RST` and
:cls:`~sire.io.parser.TRR` formats are binary, and so hold the trajectory
data at the highest precision and in the most compact format. They are
easily seekable, so give the best balance between speed and trajectory size.

The :cls:`~sire.io.parser.TRAJ` format is text-based, and so takes up a
lot of space and stores data at a fixed (lower) precision. It is included
for completeness, but is not recommended as a format for saving new
trajectories.

The :cls:`~sire.io.parser.XTC` is a compressed binary format that stores
trajectory data in binary in lower precision. This is very space efficient,
but at the cost of losing quite a bit of precision in the exact coordinate
data. It is only recommended if you want to save storage space and don't
need to recover, e.g. energies or other properties that depend
on exact coordinates.

Loading trajectories
--------------------

You load trajectories in the same way as loading any other file, e.g.
using :func:`sire.load`.

>>> import sire as sr
>>> mols = sr.load(sr.expand(sr.tutorial_url, "ala.top", "ala.traj"))
>>> print(mols.num_frames())
500

:mod:`sire` will automatically determine the format of the trajectory.
For example, here we have loaded the 500 frames from the file
``ala.traj``. This file is in :cls:`~sire.io.parser.RST` format, not
:cls:`sire.io.parser.TRAJ` format, despite the file extension.

An exception will be raised if any errors are detected in the trajectory
files. These errors can sometimes be difficult to debug. To help, you can
try to directly load the trajectory file using the appropriate parser, e.g.

>>> t = sr.io.parser.TRAJ("ala.traj")
UserWarning: SireIO::parse_error: Disagreement over the number of read bytes... 128 vs -13.
This indicates a program bug or IO error. (call sire.error.get_last_error_details()
for more info)

Trying using the correct parser gives us

>>> t = sr.io.parser.RST("ala.traj")
>>> t
AmberRst( nAtoms() = 1912, nFrames() = 500 )

You can use the parsers to read frames individually into
:cls:`sire.mol.Frame` objects, e.g.

>>> f = t.get_frame(0)
>>> print(f)

The :cls:`~sire.mol.Frame` object has functions like ``.coordinates()``,
``.velocities()`` and ``.forces()``, which can be used to get the
raw coordinates, velocities and forces data from the frame.

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


