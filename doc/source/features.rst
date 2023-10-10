========
Features
========

:mod:`sire` is a powerful Python/C++ module for loading and manipulating
molecular (predominantly biomolecular) systems.

* :doc:`Calculate alchemical free energies <tutorial/index_part06>` using
  GPU-accelerated molecular dynamics simulations.
* :doc:`Load and save <tutorial/part01/02_loading_a_molecule>`
  molecules from a number of
  :doc:`molecular file formats <tutorial/part01/06_supported_file_formats>`.
  This includes Amber, CHARMM and GROMACS files.
* This includes :doc:`loading, editing and saving trajectory files <cheatsheet/trajectory>`!
* :doc:`View molecules <tutorial/part04/04_energy_trajectories>`
  (or subsets) using `nglviewer <https://nglviewer.org>`__. This
  includes creating complex :doc:`multi-representation views <cheatsheet/view>`.
* :doc:`Create 2D structure views of molecules <cheatsheet/view>`,
  including automatically deriving bond orders, formal charges and
  stereochemistries from molecules where that information is not
  present.
* Search for molecules, atoms etc using a powerful
  :doc:`search engine <cheatsheet/search>`.
* :doc:`Convert molecules to and from RDKit <tutorial/part05/01_convert>`,
  and easily :doc:`create molecules from smiles strings <tutorial/part05/03_smiles>`
  or :doc:`generate smiles strings from molecules <tutorial/part05/03_smiles>`.
* :doc:`Convert molecules to and from smarts strings too! <tutorial/part05/04_smarts>`.
  Use smarts (and smiles) :doc:`in your searches <cheatsheet/search>`.
* :doc:`Convert molecules to OpenMM <tutorial/part05/01_convert>`, and
  run molecular dynamics and minimisation
  :doc:`via an intuitive interface <tutorial/part05/05_dynamics>`.
* :doc:`Convert molecules to and from BioSimSpace <tutorial/part05/01_convert>`,
  thereby opening the way for automatic parameterisation of ligands,
  solvation of molecules, and running of complex simulation workflows.
* :doc:`Edit molecules <tutorial/part03/02_cursors>`,
  :doc:`move them <tutorial/part04/05_movement>`, recombine them into
  new molecular systems.
* :doc:`Calculate energies <tutorial/part04/03_energies>`
  of molecular systems (or subsets) using the forcefield
  parameters read from the input files. Calculate
  :doc:`energies across trajectories <tutorial/part04/04_energy_trajectories>`
  (or subsets of trajectories).
* :doc:`Measure distances, angles, dihedrals etc <tutorial/part04/01_measure>`,
  including across
  :doc:`frames of trajectories <tutorial/part04/02_trajectory>`
  (or subsets of trajectories).
* A powerful :doc:`property system <tutorial/index_part03>` that lets you
  associate nearly any data with nearly any molecular view or sub-view.
  For example, you can assign multiple different coordinate properties
  for atoms based on their different frames or perturbation states.
* A comprehensive :doc:`units system <cheatsheet/units>` that lets you
  choose you own physical units, and that interpret units from strings
  and convert from popular units packages.
* and lots more! Take a look at the :doc:`quick overview <quickstart/index>`
  for a flavour of what :mod:`sire` can do, or take a look at
  :doc:`the tutorial <tutorial/index>` for more detail.

