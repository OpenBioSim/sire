==========
Restraints
==========

It can be useful to add restraints to the system to hold things in place.  For
example, you might want to hold the protein fixed while allowing the ligand to
move around. Or you may want to hold a ligand in place while it is being
decoupled from the simulation. Or you may want to use restraints to control
the distances between atoms or pull a molecule into a particular conformation.

You can specify the restraints you want to add to a system via the functions
in the :mod:`sire.restraints` module. These functions return
:class:`~sire.mm.Restraints` objects that contain all of the information that
needs to be passed to OpenMM to add the restraints to the system.
