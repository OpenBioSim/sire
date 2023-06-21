============================================
Part 5 - Interconverting with other Packages
============================================

Another powerful feature of :mod:`sire` is the ability to interconvert
and work with other molecular modelling python packages.

You can also convert molecules between the :mod:`sire` format and
the format of other popular molecular packages, e.g.
`rdkit <https://www.rdkit.org>`__, `openmm <https://openmm.org>`_ and
`BioSimSpace <https://biosimspace.openbiosim.org>`__.

This chapter will teach you how to do the conversion, and also how
this conversion is used in :mod:`sire` to add functionality such as creation
of molecules from smiles strings, rendering two-dimensional structure
views of molecules, performing molecular dynamics and minimisation,
and parameterising and solvating molecules.

.. toctree::
   :maxdepth: 1

   part05/01_convert
   part05/02_view
   part05/03_smiles
   part05/04_smarts
   part05/05_dynamics
