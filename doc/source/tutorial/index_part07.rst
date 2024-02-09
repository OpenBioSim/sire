===========================================
Part 7 - Merged Molecules and Lambda Levers
===========================================

In the last chapter we learned how to use :mod:`sire` to set up and
run alchemical free energy calculations. We saw how the concept of a
merged molecule, together with a lambda lever, enabled us to define
a perturbation between two end states, and then run a set of
OpenMM GPU-accelerated molecular dynamics simulations to calculate
the free energy between those states.

In this chapter we will dive deeper into the machinery of :mod:`sire`
that supports the creation and management of merged molecules. We will
see how this, plus more advanced use of lambda levers enables the
running of more complex free energy simulations. These include
those to calculate absolute binding free energies, calculate free
energies of perturbing individual residues in proteins, and
calculating free energies that involve breaking rings in ligands.

.. toctree::
   :maxdepth: 1

   part07/01_perturbation
   part07/02_levers
