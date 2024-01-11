======================
Supported file formats
======================

:mod:`sire` supports reading and writing to many common molecular file formats.
You can print the list of supported formats using

.. code-block:: python

   >>> print(sr.supported_formats())

   ## Parser dcd ##
   Supports files: DCD
   DCD coordinate/velocity binary trajectory files based on charmm / namd / x-plor format.
   ################

   ## Parser gro87 ##
   Supports files: gro, g87
   Gromacs Gro87 structure format files.
   ##################

   ## Parser grotop ##
   Supports files: top, grotop, gtop
   Gromacs Topology format files.
   ###################

   ## Parser mol2 ##
   Supports files: mol2
   Sybyl Mol2 format files.
   #################

   ## Parser pdb ##
   Supports files: pdb
   Protein Data Bank (PDB) format files.
   ################

   ## Parser pdbx ##
   Supports files: pdbx, cif
   Protein Data Bank PDBx/mmCIF format files.
   #################

   ## Parser prm7 ##
   Supports files: prm7, prm, top7, top, prmtop7, prmtop
   Amber topology/parameter format files supported from Amber 7 upwards.
   #################

   ## Parser psf ##
   Supports files: psf
   Charmm PSF format files.
   ################

   ## Parser rst ##
   Supports files: rst, crd, trj, traj
   Amber coordinate/velocity binary (netcdf) restart/trajectory files supported since Amber 9, now default since Amber 16.
   ################

   ## Parser rst7 ##
   Supports files: rst7, rst, crd7, crd
   Amber coordinate/velocity text (ascii) restart files supported from Amber 7 upwards.
   #################

   ## Parser sdf ##
   Supports files: sdf, mol
   Structure Data File (SDF) format files.
   ################

   ## Parser traj ##
   Supports files: traj, trj, crd
   Amber trajectory (ascii) coordinate or velocity files supported from Amber 7 upwards.
   #################

   ## Parser trr ##
   Supports files: trr
   Gromacs TRR (XDR file) coordinate / velocity / force trajectory file
   ################

   ## Parser xtc ##
   Supports files: xtc
   Gromacs XTC (XDR file) compressed coordinate trajectory file
   ################

Symmetric Input / Output
------------------------

One of our design principles is that molecular file input and output
is symmetrical. This means that :mod:`sire` can read in and write out the same
amount of information from a file (i.e. it can always read what it writes).

Another design principle is that information should not be lost. As much
as possible, :mod:`sire` will load and preserve all
molecular-level information it can read from a molecular file.
