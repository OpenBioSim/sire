=======================================
Part 5 - Creating and Editing Molecules
=======================================

Another powerful feature of :mod:`sire` is the ability to create and
edit molecules. For example, you can add or remove atoms, bonds etc
from existing molecules, or you can create molecules from, e.g.
`smiles strings <https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system>`__,
or assemble molecules from fragments.

You can also convert molecules between the :mod:`sire` format and
the format of other popular molecular packages, e.g.
`rdkit <https://www.rdkit.org>`__, `openmm <https://openmm.org>`_ and
`BioSimSpace <https://biosimspace.openbiosim.org>`__, and can use
the functions in those packages for molecular editing.

This chapter will teach you how to do all of the above, starting from the
simplest point, of creating a molecule from a
`smiles string <https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system>`__.

.. toctree::
   :maxdepth: 1

   part05/01_smiles
