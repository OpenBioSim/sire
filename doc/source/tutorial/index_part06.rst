============================================
Part 6 - Morphs and Alchemical Free Energies
============================================

One of the original design goals of :mod:`sire` was to support the
prototyping, setup and running of alchemical free energy simulations.

:mod:`sire` has a lot of in-built support for perturbing (morphing)
between different molecules, and to simulate and then calculate
alchemical free energy differences along those morphs.

This is supported by integration with `OpenMM <https://openmm.org>`__,
thereby letting you set up and run molecular dynamics simulations that
morph one molecular system into another.

This chapter will teach you how to set up a morph, how to run
alchemical molecular dynamics simulations, and how to calculate
alchemical free energies along morph coordinates. This includes how
to add restraints to atoms, fix atoms in space, and change these
restraints as you perturb between different molecules.

.. toctree::
   :maxdepth: 1

   part06/01_merge
   part06/02_alchemical_dynamics
   part06/03_restraints
   part06/04_alchemical_restraints
   part06/05_free_energy_perturbation
