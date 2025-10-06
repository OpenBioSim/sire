=========
Changelog
=========

There have been four main phases of development of :mod:`sire`
since the project started in 2005.

GitHub (OpenBioSim): January 2023 - current
-------------------------------------------

Development was migrated into the
`OpenBioSim <https://github.com/openbiosim>`__
organisation on `GitHub <https://github.com/openbiosim/sire>`__.

`2025.2.0 <https://github.com/openbiosim/sire/compare/2025.1.0...2025.2.0>`__ - October 2025
--------------------------------------------------------------------------------------------

* Please add an item to this CHANGELOG for any new features or bug fixes when creating a PR.

* Add end-state coordinates properties to decoupled and annihilated molecules.

* Map the end state coordinates property when constructing an ``AmberParams`` object when
  creating a molecule from a ``SOMD`` perturbation file.

* Fix evaluation of custom force energies from OpenMM XML files by correctly looping over
  interaction groups.

* Switch to new Boost Conda package naming scheme for host requirements.

* Use C++20 standard for SireRDKit plugin.

* Fixed build issues on MacOS with Clang and LTO (incorrectly specified LTO library version).

* Added a feature for applying RMSD restraints using the OpenMM CustomCVForce/RMSDForce functionality.

* Added support for reading and writing CRYST1 records in PDB files.

* Added a ``save_velocities`` kwarg to :func:`sire.save()` to allow a user to control whether
  velocities are saved to the output file. This is useful when the magnitude of the velocities
  overflows the output file format precision, or when the velocities are not needed

* Use force field bonding when working out connectivity for QM link atoms.

* Fix handling of 4- and 5-point water models when writing GROMACS topology files.

* Fix handling of link atoms for non carbon-carbon bonds where there is no bonded
  term in the MM force field.

* Added alchemical Morse potentials for experimental use in scaffold-hopping transformations.

`2025.1.0 <https://github.com/openbiosim/sire/compare/2024.4.2...2025.1.0>`__ - June 2025
-----------------------------------------------------------------------------------------

* Fixed unbound local variable when simulation lambda value is not in ``lambda_windows`` array.

* Allow negative residue numbers.

* Make sure box vectors are in reduced form before setting via the ``OpenMM`` C++ API.

* Add isobaric and grand canonical terms to reduced potential.

* Keep kappa term fixed during decoupling schedules.

* Added CMAP support when reading / writing AMBER and GROMACS topology files. The parameters
  are in the ``cmap`` property.

* Fixed compile issues for MacOS Sequoia.

* Upgraded minimum ``CMake`` supported version to 3.5.

* Reset dynamics step counters and internal clock when a crash occurs so that energies
  and frames are save with the correct time stamp.

* Add OPC water model template for interconvesion between AMBER and GROMACS formats.

`2024.4.2 <https://github.com/openbiosim/sire/compare/2024.4.1...2024.4.2>`__ - Feb 2025
----------------------------------------------------------------------------------------

* Preserve molecule name when reading from SDF format.

* Handle missing velocity data for ``AMBER`` RST7 files with only a few atoms.

* Preserve SDF metadata when converting to ``RDKIt`` format.

`2024.4.1 <https://github.com/openbiosim/sire/compare/2024.4.0...2024.4.1>`__ - Feb 2025
----------------------------------------------------------------------------------------

* Allow user to force fresh inference of stereochemistry when converting to RDKit format.

* Fix setting of positive formal charge when reading SDF files.

* Only use ``atom->setNoImplicit(True)`` inside custom RDKit sterochemistry inference function.

* Fix redistribution of excess QM charge and make it the default behaviour.

`2024.4.0 <https://github.com/openbiosim/sire/compare/2024.3.1...2024.4.0>`__ - Feb 2025
----------------------------------------------------------------------------------------

* Fixed update of triclinic box vectors in ``SOMD`` following ``OpenMM`` bug fix.

* Don't automatically save energies and frames when ``dynamics.run()`` returns.

* Improved handling of ``OpenMM`` NaN errors during dynamics.

* Restore thread state before raising exceptions in the Sire OpenMM minimiser.

* Add support for Replica Exchange with Solute Tempering (REST2) simulations.

* Null SystemTrajectory pointer when all frames are deleted.

* Handle cis/trans double bond stereochemistry values in SDF bond blocks.

* Fix handling of positive formal charge when writing SDF files.

* Fix lambda schedule discussion and plots in the tutorial.

* Added support for angle and dihedral restraints which can be used in alchemical and standard simulations.

* Allow user to compute energy trajectory over a subset of the lambda windows for each lambda.

* Fix the ``hasForceSpecificEquation`` function in the ``LambdaLever`` class so that it returns true if
  there is a default equation for the force.

* Added support for the ``OpenMM`` ``MonteCarloMembraneBarostat`` to ``SOMD``.

`2024.3.1 <https://github.com/openbiosim/sire/compare/2024.3.0...2024.3.1>`__ - December 2024
--------------------------------------------------------------------------------------------

* Fixed instantiaton of ``QByteArray`` in ``Sire::Mol::Frame::toByteArray`` and count bytes with ``QByteArray::size``.

* Increase timeout before terminating ``QThread`` objects during ``PageCache`` cleanup.

* Expose missing ``timeout`` kwarg in :meth:`dynamics.minimise()` method.

* Expose missing ``include_constrained_energies`` kwarg in minimisation function.

* Make minimisation function settings consistent across API.

* Reload TorchQMForce module if OpenMM platform changes.

* Fix calculation of non-equilibrium work when starting from QM state.

* Fix stereochemistry in RDKit molecule conversion.

* Fixed :func:`sire.restraints.get_standard_state_correction` to be consistent with the definition of
  the "force constanst" as ``F = 2 ** k ** x`` (rather than ``F = k ** x``). Updated docstrings and
  restraints documentation to make it clear how the force constants are defined.

* Fix thread safety issue in Sire OpenMM minimiser.

`2024.3.0 <https://github.com/openbiosim/sire/compare/2024.2.0...2024.3.0>`__ - October 2024
--------------------------------------------------------------------------------------------

* Print residue indices of perturbed water molecules to SOMD1 log to allow
  for easier debugging.

* Add support for creating Na+ and Cl- ions as a means of generating templates for
  uses as alchemical ions.

* Fix ``sire.morph.merge`` function when one molecule is a monatomic ion. This prevents
  the attempted rigid-body alignment, which would fail due to there being too few
  degrees of freedom.

* Remove ``sire.move.OpenMMPMEFEP`` wrappers from build when OpenMM is not available.

* Set ``IFBOX`` pointer to 3 for general triclinic boxes in ``sire.IO.AmberPrm`` parser.

* Only exclude nonbonded interactions between ``from_ghost`` and ``to_ghost`` atoms
  if they are in the same molecule. This prevents spurious intermolcular interactions
  between molecules containing ghost atoms, e.g. a ligand and an alchemical water.

* Add Docker support for building wrappers on Linux x86.

* Port SOMD1 Boresch restraint implementation to PME code. (This feature was present
  in the reaction field implementation, but not for PME.)

* Port SOMD1 torsion fix to PME code. (This had been fixed for the reaction field implementation,
  but not for PME.)

* Fix issues with ``atomtype`` and ``atom`` records for dummy atoms in GROMACS topology files.

* Fixed buffer overflow when computing molecule indices to excluded to/from
  ghost atom interactions which caused corruption of the exclusion list.

* Fixed calculation of ``delta^2`` in soft-core Couloumb potential.

* Exclude to/from ghost atom interactions from the ``ghost_14ff``. Exclusions were
  already added to the ``ghost_ghostff``, but not the ``ghost_14ff``.

* Fixed description of soft-core alpha parameter in :doc:`tutorial <tutorial/part07/03_ghosts>`.

* Added debugging function to evaluate custom forces in OpenMM XML files. This
  allows a user to decompose the pair-wise contribtions to the custom OpenMM
  forces created by :mod:`sire`.

* Added a timeout to the OpenMM minimiser function. This gives the user a single tunable
  parameter to control roughly how long a minimisation should last before being aborted.

* Exposed missing pickle operator on the ``LambdaLever`` class.

* Fix bug setting custom nonbonded parameters for ghost atoms used in
  positional restraints in OpenMM.

* Fix exchange probability equations in ``sire.morph.replica_exchange`` function.

* Fix calculation of energy change following final constraint projection
  after energy minimisation. Previously the energy change was calculated from
  the final step of the minimisation, rather than the change in energy
  following the application of the constraints.

* Clear internal OpenMM state from dynamics object during minimisation,
  preventing the previous, pre-minimisation, state from being used when
  ``get_state()`` is called.

* Add support for QM/MM simulations using OpenMM. This uses the recent ``CustomCPPForceImpl``
  introduced in OpenMM 8.1 to allow an interface between OpenMM and external
  QM or ML codes. We support a generic Python callback interface and a ``Torch``
  based interface for ML models. This is documented in the new :doc:`tutorial <tutorial/part08>`.

* Reinitialise OpenMM context if constraints change when setting lambda. Updating
  constraints in an OpenMM system does not update the associated data structures
  in the context. A full reinitialiasation is required.

* Give custom OpenMM forces meaningful names. This makes it easier to parse OpenMM
  XML files and debug custom forces, particularly when multiple forces of the same
  type are present.

`2024.2.0 <https://github.com/openbiosim/sire/compare/2024.1.0...2024.2.0>`__ - June 2024
-----------------------------------------------------------------------------------------

* Correctly set the ``element1`` property in ``sire.morph.create_from_pertfile``.

* Added mising :meth:`~sire.vol.TriclinicBox.maximum_cutoff` method so that
  the cutoff is set correctly when creating a :obj:`~sire.system.ForceFieldInfo`
  object.

* Added a :class:`sire.base.PageCache` class which can be used to cache and
  restore objects to memory pages which are automatically paged to and from
  disk as needed. This lets you work on data that can't fit in memory.

* Updated the way that :class:`sire.system.System` objects hold the
  set of temporary frames in a trajectory. Rather than each molecule holding
  its own temporary frame, now the :class:`~sire.system.System` object holds
  a ``SystemTrajectory`` object. This holds the frame data for all molecules
  in the :class:`~sire.system.System` in a single binary array. The data
  for this array is paged to disk as needed via the above
  :class:`sire.base.PageCache` class. This both significantly speeds up
  processing of these temporary frames, and ensures that long simulations
  do not fill memory, causing the system to crash. In addition, the
  ``SystemTrajectory`` object is NOT streamed to a S3 file. This means that
  the S3 file (used normally for restarts) won't grow unbounded with
  temporary frames, meaning that it is safe to create restarts of
  long-running simulations. Note that this does mean that the temporary
  directory is lost. You **must** save the trajectory to a file at the
  end of your simulation or it will be lost. You can do this using the
  standard trajectory save functions, e.g.
  ``sire.save(mols.trajectory(), "output", format=["PRMTOP", "RST"])``.

* Added code that automatically excludes non-bonded interactions between
  from_ghost and to_ghost atoms in the OpenMM layer. This is to prevent
  crashes caused by poor interactions between from_ghost atoms appearing
  over the top of to_ghost atoms during a perturbation where one group
  is grown over another.

* Ignore BioSimSpace format position restraint include directives when
  parsing GROMACS topology files.

* Added a map option (fix_perturbable_zero_sigmas) to prevent perturbation of
  the Lennard-Jones sigma parameter for ghost atoms during alchemical free energy simulations.

* [CHANGE IN BEHAVIOUR] - added code that ensures that, when editing molecules,
  the CGAtomIdx order will always follow the AtomIdx order of atoms. This is
  because a lot of code had implicitly assumed this, and so it was a cause
  of bugs when this wasn't the case. Now, when you edit a molecule, on committing,
  the orders will be checked. If they don't agree, then the CutGroups will be
  reordered, with atoms reordered as necessary to make the CGAtomIdx order match
  the AtomIdx order. If this isn't possible (e.g. because atoms in CutGroups
  are not contiguous), then the molecule will be converted to a single-cutgroup
  molecule, with the atoms placed in AtomIdx order. As part of this change,
  the merge code will now also ensure that added atoms are added with the
  correct AtomIdx, rather than added as the last atoms in the molecule. This
  is also more natural. This fixes issue #202.

* Added the "center" keyword argument to the ``make_whole`` functions of
  :class:`~sire.mol.Cursors`, :class:`~sire.mol.CursorsM` and
  :class:`~sire.system.System` (as well as to the legacy System class).
  Also allowed the constructor of :class:`~sire.maths.Vector` to recognise
  ``origin`` and ``zero`` as arguments, meaning you can write
  ``cursor.make_whole(center="origin")``. This fixes issue #199.

`2024.1.0 <https://github.com/openbiosim/sire/compare/2023.5.2...2024.1.0>`__ - April 2024
------------------------------------------------------------------------------------------

* Dropped official builds and support for Python 3.9, and added official
  builds and support for Python 3.12. Note that MacOS builds are currently
  3.10 and 3.11 only, due to missing dependencies. This will be fixed
  in upcoming point releases.

* BREAKING CHANGE: Updated the API of :class:`sire.cas.LambdaSchedule` so that
  you have to use named arguments for many of the functions (e.g.
  :meth:`~sire.cas.LambdaSchedule.set_equation`). This is because the addition
  of force levers (as described below) made positional arguments ambiguous,
  and we wanted to make the API more consistent. This is a breaking change,

* Added the ability to customise the lambda schedule applied to a lambda lever
  so that you can use different equations for different molecules and
  different forces in the OpenMM context. This gives a lot of control over
  how forcefield parameters are scaled with lambda. Specifically, this is used
  to add support for calculating absolute binding free energies.
  This is described in the new :doc:`tutorial chapter <tutorial/index_part07>`.

* Exposed the underlying :class:`~sire.legacy.Convert.PerturbableOpenMMMolecule`
  class, which can be created from a merged molecule via
  ``mol.perturbation().to_openmm()``. This lets you easily see which parameters
  are changing between the reference and perturbed states. This is described
  in the :doc:`tutorial <tutorial/part07/01_perturbation>`.

* Added the ability for the lambda schedule to show how it will actually
  act to perturb the parameters of the
  :class:`~sire.legacy.Convert.PerturbableOpenMMMolecule` molecule.
  This is described in the :doc:`tutorial <tutorial/part07/02_levers>`.

* Added support for reading older somd-style pertfiles, and creating
  merged molecules from these. This is described in
  the :doc:`tutorial <tutorial/part07/05_pertfile>`.

* Added "not-perturbable" constraints so that bonds and angles that change
  with lambda are not perturbed. As part of this, have also added a
  ``dynamic_constraints`` option that lets constrained bonds update with
  lambda, so that they are set to the length corresponding to r0 at that
  lambda value. Have also changed the constraints so that bonds will be
  constrained to their r0 value, rather than their current length.
  These constraints are ``X-not-perturbed``, meaning that it constrains
  all ``X``, except for bonds or angles involving perturbed atoms. Or
  ``X-not-heavy-perturbed``, meaning that it constrains all ``X``, except
  for bonds or angles involving perturbed atoms, unless they involve a
  hydrogen in any end state. The code to detect hydrogens has been improved,
  now looking at mass, element and ambertype. There are options to control
  this, described in the :doc:`OpenMM detailed guide <cheatsheet/openmm>`.

* Added more automatic conversions, so that string will more readily auto-convert
  to units where possible. Also added a ``sire.v`` function to make it easier to
  create vectors of units, e.g. ``sire.v("1.0A", "2.0A", "3.0A")`` will create
  a ``sire.maths.Vector(1, 2, 3)``, while ``sire.v([3, 4, 5], units="A ps-1")``
  will create a ``Velocity3D``. This is documented in the units cheat sheet.

* You can now set the background color of a 3D view using the ``bgcolor="..."``
  keyword. This is documented in the view cheat sheet.

* MacOS/ARM64 now includes AmberTools and Gromacs dependencies when built
  for BioSimSpace (matching MacOS/X64 and Linux).

* Updated the electrostatic softening potential to have an additional
  ``shift_coulomb`` parameter, so that you can control how much the
  distance is increased by the alpha softening parameter. This was
  the equivalent of 10 Å, but has been set as default to 1 Å to match
  the value used in somd.

* Added support for LJ 12-6-4 potentials, plus the ability to read and write
  LJ parameter exceptions to Amber topology files. This fixes issue #125.

* Added peek support to the datastream reader, so that it can recover
  when it doesn't find the magic value it expects on reading.

* Added functionality to SparseMatrix to make it easier to detect when
  non-default values have been added, and also to set up a matrix which
  has a concept of unset values.

* Added a ``to_same_molecule`` argument to the ``mol.extract()`` function,
  so that it is possible to keep the same molecule number for the extracted
  molecule. As part of this, also relaxed the requirement that the
  ``mol.update(...)`` function can only be called if the molecule layout
  is not changed. You can now update even if you have changed the numbers
  of atoms, residues etc. The ``to_same_molecule`` argument is default False,
  so as not to change any existing behaviour.

* Added lots of convenience functions to ``sire.morph``, as described in the
  :doc:`new tutorial <tutorial/index_part07>`. Functions include
  linking to the reference or perturbed states for all molecules, or extracting
  all of the reference or perturbed states of all molecules. Also I've added
  functions for zeroing ghost torsions and creating molecules from pertfiles.
  As part of this, I added an ``auto_commit`` argument to the
  Perturbation ``link_to_reference`` and ``link_to_perturbed`` functions,
  which defaults to True. This is a change in behaviour, but it makes the
  API much easier to use. If you are affected by this, please let us know.
  It was a little-used part of the code, with the main use case being the
  replacement with the easier ``sire.morph.link_to_XXX`` functions.

* Added ability to create merge molecules for relative free energy calculations
  via :func:`sire.morph.merge` and :func:`sire.morph.match`. This is
  described in the :doc:`tutorial <tutorial/part07/04_merge>`.

* Added ability to create merge molecules for absolute free energy
  calculations via :func:`sire.morph.decouple` and
  :func:`sire.morph.annihilate`. This is described in the
  :doc:`tutorial <tutorial/part07/06_decouple>`.

* Added support for residue perturbations and also for mutating residues
  and parts of molecules using a new "copy and paste" algorithm. This is
  described in the :doc:`tutorial <tutorial/part07/07_residue>`.

* Exposed the ``SOMMContext``, ``PerturbableOpenMMMolecule``,
  ``OpenMMMetaData`` and ``LambdaLever`` classes to Python, as part of the
  new ``sire.convert.openmm`` module. These are useful if you want more
  control over OpenMM from Python. In particular, the ``PerturbableOpenMMMolecule``
  class lets you see all of the parameters that go into the OpenMM forces
  that involve perturbable molecules. There are useful functions that can
  be used to get changing parameters as dataframes, which is really useful
  for debugging. These are described in the :doc:`new tutorial <tutorial/index_part07>`.

* Preserve user atom names when writing to PDB format.

* Updated the :class:`~sire.mol.Cursor` so that it is easier to get and
  set the expression used for the potential energy (using the
  ``get_potential`` and ``set_potential`` functions).

* Fixed compile error using Python 3.12. This fixes issue #147.

* Optimised the OpenMM minimisation code and making it more robust.
  This includes vectorising for Apple Silicon and adding more tests for
  convergence so that we can have more confidence that the structures
  output are sensible. Also made sure that optimised compilation (-O3) is
  used for all of the plugins (SireOpenMM, SireGemmi and SireRDKit).
  They were previously compiled with wrapper options (e.g. -Os).
  Minimisation now gives better progress updates, using a progress
  bar to show progress towards the maximum number of iterations.
  This has been reduced to 1500 by default. Also, if the minimisation
  fails to create a structure that obeys constraints on the first pass,
  then the minimisation is repeated, with the maximum number of
  iterations reset. If it fails again, then this structure, with
  constraints re-applied, is returned.

* Added more support for Boresch restraints. Specifically, :func:`sire.restraints.boresch`
  now supports the specification of equilibrium values, uses different default force
  constants, and warns the user if the restraints are likely to be unstable.
  :func:`sire.restraints.get_standard_state_correction` was implemented for Boresch
  restraints. Tests were added for restraint creation and for the standard state
  correction. Boresch restraints were added to :doc:`tutorial <tutorial/part06/03_restraints>`.

* Added power support to GeneralUnit in python. You can now raise a unit value
  to a valid power, e.g. ``sr.u("5A")**2``. You can also square root via a new
  ``.sqrt()`` function on the unit. There is also a new ``sire.sqrt`` function
  that will automatically called ``obj.sqrt()`` if that exists, or will fall
  back to ``math.sqrt`` if not. This implements wishlist item #176.

* Conversion to and from RDKit now preserves atoms and residue names and
  numbers. This makes used of AtomPDBResidueInfo in RDKit to populate metadata
  when the RDKit molecule is created. On conversion to sire, the atom monomer
  info will be checked. If it is simple, then only the atom name will be
  obtained. If it is a AtomPDBResidueInfo, then the atom name and number,
  and residue name and number (plus any chain information) will be extracted.
  If no atom name is set, then the value of the property
  "molFileAlias" will be checked. This implements wishlist item #168.

* Implemented the ``auto-bonds`` constraint, which automatically chooses
  bonds to constrain by comparing the estimated vibrational frequency
  against the simulation timestep multiplied by the factor
  ``auto_bonds_factor`` (defaults to 10).
  This is described in the :doc:`OpenMM detailed guide <cheatsheet/openmm>`.
  This implements wishlist item #185.

* Fixed a bug in the algorithm used to infer bond order when converting to
  RDKit format. This fixes issue #177.

* Fixed a bug in the :class:`~sire.legacy.Convert.LambdaLever` class where
  it was not using the stage-specific value of lambda when using multiple
  stages where one or more stages contained a standard morph equation.

`2023.5.2 <https://github.com/openbiosim/sire/compare/2023.5.1...2023.5.2>`__ - March 2024
------------------------------------------------------------------------------------------

* Fix bug that disabled the ``DEBUG`` log level from the global logger.

* Fixed bug in :class`sire.legacy.Mol.ResIdxAtomCoordMatcher` by ensuring
  that we only compare residues with the same number of atoms.

* Added an ``AtomCoordMatcher`` to match atoms by coordinates in two selections.

* Added ``map`` support to writing perturbable Gromacs topology files. This
  enables the user to specify which perturbable properties to use,
  e.g. ``map={"dihedral0": "dihedral_a", "dihedral1": "dihedral_b"}``.

* Code can now detect when an Amber PRMTOP file has discontiguous molecules,
  and thus when atoms are reordered after load. This information is passed
  to subsequent frame file parsers that are loaded at the same time, so
  that they are able to reorder the frames before being added to the atoms.
  This happens transparently, so that the user doesn't have to worry about
  the reordering. This fixes issue #164.

* Fixed a bug where the SDF parser would wrongly try to parse Amber RST7 files that
  weren't immediately recognised as such. The fix adds ``.inpcrd`` as a recognised
  extension for Amber RST7 files, and changes the scoring logic of the SDF parser
  to equal the number of molecules times the number of atoms.

* Switched to using the SYBL atom type to infer the element of an atom
  when reading Mol2 files. This is more robust than using the atom name.
  Fixes issue #166.

* Made it easier to convert from strings to elements. Added the ability to
  customise the list of elements that are considered biological. This
  fixes issue #170.

`2023.5.1 <https://github.com/openbiosim/sire/compare/2023.5.0...2023.5.1>`__ - January 2024
--------------------------------------------------------------------------------------------

* Added a ``.dynamics().step(num_steps)`` function to make it easier to quickly run
  steps of OpenMM dynamics with minimal overhead (this directly called ``Integrator.step()``)

* Some optimisations to the OpenMM layer that make full use of the
  experimental "updateSomeParametersInContext" functions.

* Updated gemmi to 0.6.4, so that it can be default-enabled on all supported platforms.
  As part of this, had to change the version of the compilers used on Windows and Linux
  to make the conda packages. Windows now uses Visual Studio 2022 instead of 2017,
  and Linux now uses GCC 12.3.0 instead of GCC 13.

`2023.5.0 <https://github.com/openbiosim/sire/compare/2023.4.2...2023.5.0>`__ - December 2023
---------------------------------------------------------------------------------------------

* Added a new :mod:`sire.options` module that contains new
  :class:`sire.options.Option` objects to represent configurable options.
  These include documentation, and make it easier to validate and expose
  possible values of configurable options. The API docs for
  :class:`~sire.options.Option` shows how to create your own Option type.
  The unit test in ``tests/options/test_options.py`` show how to use
  the options. This is integrated into the sire/OpenMM layer.

* Extended the ``.atom(s)``, ``.residue(s)``, ``.bond(s)`` and all other
  indexing functions so that you can pass in an existing view or views as
  the key. This lets you look up views in a container by other views, e.g.
  ``mols.bond(mols.atoms()[0], mols.atoms()[1])`` would return the bond
  between the first two atoms in the container ``mols``. Also added
  a ``error_on_missing`` flag to the ``atoms``, ``residues``, ``bonds`` etc
  functions, so that you get a ``KeyError`` exception if there is no match,
  and ``error_on_missing`` is ``True``. For example,
  ``mols.atoms("element C", error_on_missing=True)`` would raise an exception
  if there are no carbon atoms in this container. This is default ``False``
  to keep existing behaviour, but we would recommend setting this to ``True``
  and would like to change the default in the future.

* Added :func:`sire.convert` support for converting between :mod:`sire`
  objects and `Gemmi <https://gemmi.readthedocs.io>`__ objects. This
  has allowed us to support reading and writing of PDBx/mmCIF files.
  We've updated :func:`sire.load` to automatically choose PDBs/mmCIF
  files if gemmi-support is available. We've also added support for the
  new-style PDB codes (e.g. "pdb_00003nss" instead of "3NSS"). Note that
  this needs a custom Gemmi package build, where "shared libraries" are
  turned on. This should be available from conda-forge in 2024, but for now,
  you will need to clone the `Gemmi feedstock <https://github.com/conda-forge/gemmi-feedstock>`__
  and build the conda package yourself. You will then need to recompile
  sire from source. We will release 2023.5.1 as a conda package once
  the conda-forge Gemmi package with shared library support is available.

* Optimised the ``LambdaLever`` class so that it caches the forcefield parameters
  calculated at different lambda values. This means that we don't have to
  re-calculate the parameters at each lambda update step. This is a
  significant speed-up for alchemical free energy simulations.

`2023.4.2 <https://github.com/openbiosim/sire/compare/2023.4.1...2023.4.2>`__ - December 2023
---------------------------------------------------------------------------------------------

* Fixed use of ``QString::compare`` when comparing molecular properties during
  a water topology swap.

* Have :class:`~sire.io.parser.RST7` return a list of angles from the
  ``box_angles()`` function, rather than a :class:`~sire.maths.Vector`.
  This prevents the confusing behaviour where the angles are wrongly
  shown in units of angstroms... This fixes issues #106.

* Added a new :func:`sire.maths.rotate` function, and added the option
  (default True) to rotate velocities as well as coordinates when usign
  a cursor to rotate molecule views. This fixes issue #103.

* Fix validation of ``perturbable_constraint`` dynamics option when the string
  includes hyphens. This fixes issue #130.

* Fix streaming of :class:`~sire.vol.TriclinicBox` objects. This fixes issue #128.

* Fix the sire to OpenMM conversion so that null LJ parameters will never have
  a zero sigma value. They will either be sigma=1/epsilon=0 for non-perturbable
  atoms, or sigma=1e-9/epsilon=1e-9 for perturbable atoms.

* Now catch ``std::bad_alloc`` and raise it as a ``MemoryError``. This
  means that we can catch out-of-memory errors and raise a more
  informative exception.

* Fixed the bug where the wrong return type from ``.minimisation()`` and
  ``.dynamics()`` was returned. This fixes issue #137.

* Fixed the bug where the cutoff would not be set correctly if a string
  was passed. You can now do ``mol.dynamics(cutoff="10A")`` or
  ``mol.dynamics(cutoff="infinite")`` and it will be processed correctly.
  This also required adding a ``map.unset("key")`` option to ``PropertyMap``,
  to make it easier to unset mapped properties.

`2023.4.1 <https://github.com/openbiosim/sire/compare/2023.4.0...2023.4.1>`__ - October 2023
--------------------------------------------------------------------------------------------

* Fixed regression introduced in 2023.4.0 that meant that removed the constraints
  from water molecules that had no internal bonds. These waters would blow up
  as there was nothing holding them together. The need for these constraints is
  now better detected and explicitly added.

* Significantly sped up the OpenMM layer by checking for similar constraint lengths
  and matching them all to be the same (within 0.05 A for calculated constraints,
  e.g. unbonded atoms or angle constraints) or to R0 for bonds where the bond
  length is within 0.1 A of R0 and the molecule isn't perturbable.

* Added a custom minimiser that is based on OpenMM's LocalEnergyMinimizer,
  but that copes better with exclusion errors, and that has deep integration
  with the progress bar / interuption system.

* Fixed a bug where the exclusions and exceptions were mismatched for the
  OpenMM CPU platform, leading to exclusion errors.

* Fixed an issue where the vacuum dynamics and minimisation simulations still
  had a spurious periodic box added when ``.commit()`` was called.

`2023.4.0 <https://github.com/openbiosim/sire/compare/2023.3.0...2023.4.0>`__ - October 2023
--------------------------------------------------------------------------------------------

* Added ``closest`` and ``furthest`` keywords to enable searching for the n closest
  or furthest views. This is very general, and is described in the
  :doc:`detailed search guide <cheatsheet/search>`. Searches such as
  ``closest 10 waters to protein`` or
  ``furthest (residue in protein) from water`` are supported.

* Added a :func:`sire.minimum_distance` function to calculate the minimum
  distance between atoms of two views.

* Added support for perturbable molecules to the OpenMM converter. Have addded
  ``LambdaLever`` and ``LambdaSchedule`` classes that can be used to control
  how forcefield parameters are changed with lambda. These levers change
  the parameters in the OpenMM context, enabling simulations at different
  values of lambda to be performed. This is initial functionality which
  will be documented and expanded by subsequent PRs.

* Added support for softening potentials used to smooth creation and
  deletion of ghost atoms during alchemical free energy simulations.
  Added a new ``sire.morph`` module that includes functions that should
  make it easier to set up, view and control morphs (perturbations).

* Forced all new-style modules to import when `sr.use_new_api()` is called.
  This will make it easier to use sire with multiprocessing.

* Added option to allow GROMACS water molecules to be flagged as crystal waters.
  This means that they will be ignored by ``gmx genion`` when choosing water
  molecules to replace with ions.

* Added the ability to align trajectories and views against molecule views
  or containers. Added the :class:`sire.mol.AtomMapping` class to control
  how to map from atoms in one group of molecules to another. This can
  be used to align trajectories and views against atoms / molecules that
  are not part of that trajectory.

* Added the :func:`sire.mol.TrajectoryIterator.rmsd` function to make it
  easier to calculate RMSDs across trajectories. The RMSD can be calculated
  against all atoms, a subset of atoms, or even against a different
  set of atoms that are matched via an :class:`~sire.mol.AtomMapping` object.
  Full details in the :doc:`tutorial <tutorial/part04/02_trajectory>`.

* Significantly optimised the loading of trajectory frames and of updating
  properties in molecules. Switched from ``CentralCache`` to a new
  ``LazyEvaluator`` class that uses ``tbb::collaborative_call`` to
  lazy-calculate the results of functions in a thread-safe and
  thread-cooperative manner. Moved ``PropertyMap`` to use a shared
  pointer to assigned properties (removing costs of unnecessary
  allocations and deallocations) and added ``update`` and ``updateFrom``
  functions to ``Properties`` and ``MoleculeData`` so that properties
  can be updated in place, thereby minimising new/free.

* Fixed a bug that prevented ``mols.trajectory().view()`` from working.
  You can now view trajectory subsets again, e.g. ``mols.trajectory()[0:5].view()``.

* Updated ``FreeEnergyAnalysis.py`` to be compatible with both the new pymbar 4 API
  and the old pymbar 3 API.

* Added support for restraints to the OpenMM dynamics layer. Initial tested
  support for positional and distance/bond restraints is included, as well
  as experimental support for Boresch restraints. The restraint are documented
  in the :doc:`tutorial <tutorial/part06/03_restraints>`. This also documents
  new code to let you specify atoms that should be fixed in space.

* Added support for alchemical restraints to the OpenMM dynamics layer.
  This lets you scale restraints as part of a λ-coordinate. This is
  documented in the :doc:`tutorial <tutorial/part06/04_alchemical_restraints>`.
  Restraints can be named, meaning that you can scale different restraints
  at different stages and by different values across the λ-coordinate.

* Added an :class:`~sire.maths.EnergyTrajectory` class that lets us record the
  energy trajectory along a dynamics simulation. This includes recording
  energies at different λ-windows to that being simulated, thereby providing
  the raw data for free energy calculations. By default the
  ``EnergyTrajectory`` is returned to the user as a pandas DataFrame.

* Added the ability to export an :class:`~sire.maths.EnergyTrajectory` as
  an alchemlyb-compatible data frame. Added :func:`sire.morph.to_alchemlyb`
  to convert lots of ``EnergyTrajectory`` objects (or files containing
  s3 streams) into a single alchemlyb-compatible data frame that is
  ready for analysis. You can now calculate relative hydration and binding
  free energies and analyse the results using alchemlyb. This is documented
  in the :doc:`tutorial <tutorial/part06/05_free_energy_perturbation>`.

* Added a :func:`sire.morph.repartition_hydrogen_masses` to make it easier to
  repartition hydrogen masses during alchemical free energy simulations.
  Set the default mass factor to 1.5 to support a 4 fs timestep with the
  default ``LangevinMiddleIntegrator``.

* Added support for an Andersen thermostat in the OpenMM dynamics layer.

* Added support for scaling intramolecular non-bonded scale factors to the
  ``LambdaLever``, so that we have rudimentary support for perturbations
  that involve bond breaking and forming.

* Added support to somd for one or more "permanent" distance restraints. These
  are distance restraints that are always applied, and are never scaled by λ.
  This allows the release of all other distance restraints to a single
  harmonic or flat-bottomed restraint. When the ligand is fully decoupled,
  the free energy of release of the single remaining restraint can be
  computed without simulation. See
  <https://pubs.acs.org/doi/10.1021/acs.jctc.3c00139> for more details.

`2023.3.2 <https://github.com/openbiosim/sire/compare/2023.3.1...2023.3.2>`__ - September 2023
----------------------------------------------------------------------------------------------

* Made sure that a title is written to an AmberRst file, even if the system
  has no name (issue #99).

* Modularise the :class:`~sire.vol.TriclinicBox` lattice rotation and reduction functionality
  and make both optional. (PR #102).

* Updated default units so that units of pressure default to printing out in units of atmospheres

`2023.3.1 <https://github.com/openbiosim/sire/compare/2023.2.3...2023.3.1>`__ - July 2023
-----------------------------------------------------------------------------------------

* Fixed a bug in ``analyse_freenrg`` which produced incorrect TI results
  when not all lambda windows were run for equal lengths of time.

* Make sure atom serial number in PDB files are capped when renumbering when
  TER records are present.

* Fixed a bug in the AmberRst parser where velocities were written with the wrong
  unit (A ps-1 instead of AKMA time). Also added the correct labels to the AmberRst file.

* Fixed a bug where outputs from legacy script would be written with base physical
  units, rather than prettier internal or SI units.

* Fixed a bug in the writing of DCD headers, meaning that the files couldn't be read
  by other DCD reader software (written non-compliant header)

* Fixed a bug in the trajectory measure code, where the ProgressBar class was
  not being properly imported (`fix_88 <https://github.com/OpenBioSim/sire/issues/88>`__).

* Fixed a deadlock in the file trajectory loading code. This was because multiple threads
  trying to read the same frame lead to starvation of the thread that had progressed to
  read the frame. Now a single thread loads the frame, with subsequent threads using
  this cached load (`fix_88 <https://github.com/OpenBioSim/sire/issues/88>`__).

* Optimised the speed of viewing large molecules in NGLView, plus of searching
  for water molecules. Added a new ``is_water`` function. Optimised the
  find function in ``SelectorM<T>`` so that it is not an O(N^2) search. It
  is now roughly O(N), using a hash to lookup at the molecule level, so that
  we don't have to test individual atoms.

* Fixed ``StandardStateCorrection``. This stopped working after
  the commit https://github.com/OpenBioSim/sire/commit/e2e370940894315838fb8f65e141baaf07050ce0,
  because not all required changes were included.

* Fix for crash when not passing a map to the SelectorImproper constructor

* Fix for crash when checking a list of atoms rather than a list of molecules

`2023.3.0 <https://github.com/openbiosim/sire/compare/2023.2.3...2023.3.0>`__ - June 2023
-----------------------------------------------------------------------------------------

* Added alignment and smoothing options to trajectory views (and trajectory processing).
  You can now align a trajectory against any search string, can wrap molecules into
  boxes, and can smooth coordinates across multiple frames. This is described in the
  new :doc:`detailed trajectory guide <cheatsheet/trajectory>`.

* Rewrote many of the "frame" trajectory parsers, and consolidated all of these
  parsers on top of the new "Trajectory" / "Frame" design. All trajectory frames
  are now streamed on demand from disk, and are not saved in memory (except
  for a small cache). Loading and scanning through the frames of a trajectory
  is massively optimised, and now quite fast :-)

* Used the same framework to all streamed saving of trajectory frames to disk.
  Trajectories can be written in parallel. Frame data comes either directly
  from the underlying molecular data, or can come from the result of aligning,
  wrapping or smoothing the trajectory. Because loading and saving is streamed,
  this means we can easily and quickly convert one trajectory format to another
  without consuming much memory. Indeed, parallel streaming means that we
  can write multiple new formats at the same time.

* As part of this, we now support a wider range of trajectory file formats.
  We support Amber RST (NetCDF), Amber TRAJ, Gromacs TRR, Gromacs XTC and DCD.

* We have also added code to allow any single-frame format to be used to load
  and save trajectories. This is a little experimental still, but supports
  writing out the frames of a trajectory to several individual files in
  a directory. Sire will automatically recognise these directories on load,
  and will stream-load the frames as needed.

* Added full smarts and smiles searching support, including proper returning
  and querying of sub-structure matches. This is described in the
  :doc:`search guide <cheatsheet/search>` and
  :doc:`new tutorial <tutorial/part05/04_smarts>`.

* Re-worked the progress bar, and how sire communicates from the C++ layer
  up to the Python layer. Progress bars are now created in C++ and are
  thread-safe. They propogate up to the Python layer, meaning that they
  render and update even though the C++ function is running (often in
  parallel). The Python GIL is correctly released and recovered around
  these functions and around progress bar updates. This makes it much easier
  to use progress bars, as well as making it easier to interupt long-running
  C++ functions (they catch and respond to the break signal in a signal
  handler that alerts the progress bar, so it raises an interupt_exception
  at the next update). This is all thread safe, meaning that child threads
  can create progress bars that become children of their parent's bars,
  with them all rendering correctly. The progress bars are currently used
  for the trajectory saving code, and the OpenMM MD and minimisation code.
  We will develop them further (and make them prettier) in future
  releases.

* Added better handling of :class:`~sire.system.System`, so that we don't
  lose system-level properties or the system name during manipulations.
  This was achieved by better calling these functions on the
  underlying :class:`sire.legacy.System.System` object, and not
  dropping straight to :class:`~sire.mol.SelectorMol`. Also added
  :func:`~sire.system.System.space` and :func:`~sire.system.System.time`
  functions (with appropriate ``set_space`` and ``set_time``) to more
  easily see and change the system space and time.

* Added "shared properties" to :class:`~sire.system.System`. These are
  properties which will be automatically copied into contained
  molecules (and kept up to date if they are changed). The
  ``space`` and ``time`` properties are default "shared". You can
  add or remove shared properties via new functions
  :func:`~sire.system.System.add_shared_property` and
  :func:`~sire.system.System.remove_shared_property`.

* Cleaned up the sanitisation of molecules generated by smiles strings.
  This now raises an exception if the molecule can't be sanitised. You
  can switch this off by passing ``must_sanitize=False`` to
  :func:`sire.smiles`, thereby only running the sanitisation steps
  that pass.

* Improved functionality of harmonic restraints in openMMMD. Each
  restrained atom will now have to a corrsponding dummy atom,
  with the location of this dummy atom restraining the real atom.
  This makes restrained systems consistent in NPT regimes. Provided
  that a modified system containing the dummy atoms is given, the argument
  ``use restraints = True`` can be added to a SOMD ``.cfg`` file, along with
  the argument ``restrained atoms`` containing a dictionary of dummy atom
  numbers along with the numbers of the corresponding real atoms
  (``{dummy_atom_num:real_atom_num}``).

* Added a new units grammar and parser, so that we can robustly
  read physical quantities (units) from strings. This is a complete
  grammar, meaning that the full range of physical units, SI prefixes,
  long and short forms, unicode and ASCII representations are supported.
  A convenient :func:`sire.u` function has been added to easily convert
  its arguments to :class:`sire.units.GeneralUnit`, e.g.
  ``timestep = sr.u("25ps")`` or ``m = sr.u("25km").to("miles")``.
  This even supports automatic conversion
  from other units systems (e.g. you can pass `pint` units to ``sr.u``
  to convert to ``GeneralUnit``). This is described in full in the
  new :doc:`units detailed guide <cheatsheet/units>`.

* We have begun to add automatic conversion from strings to unit-quantities
  (or from any unit system) to functions. Currently the dynamics functions
  are supported, e.g. you can type ``d = mols.dynamics(timestep="4fs")``
  and ``d.run("500ps")`` (or even, ``d = mols.dynamics(timestep={pint time})``).
  We will add more in the next release, with the ambition that every function
  that accepts a unit argument will automatically convert from ``pint``,
  a string, or any other supported units framework.

* As part of this, we have also updated the way physical units are printed.
  Units will now always be printed in the default format specified by
  the user, with default base units used for any composed unit that
  has not been specified. You set a default unit using
  :func:`sire.units.set_default_unit` or :func:`sire.units.set_default_units`,
  e.g. ``sr.units.set_default_unit("nm")`` would change the default
  length unit to ``nm``. You can set combined units, e.g.
  ``sr.units.set_default_unit("kcal mol-1 A-2")`` would set the default
  bond force constant units to ``kcal mol-1 A-2``. The framework automatically
  works out the unit, and will print out any value in that unit if the future
  in the default unit, with the string supplied by the user as the unit name.
  This is all described in the :doc:`units detailed guide <cheatsheet/units>`
  (including :func:`sire.units.set_si_units` and :func:`sire.units.set_internal_units`).

* Added :func:`~sire.mol.SelectorMol.make_whole` functions to all views,
  so that molecules can be recombined after being split across periodic
  boundaries. You can automatically make molecules whole on load by
  passing ``{"make_whole": True}`` as a ``map`` to :func:`sire.load` or
  the :func:`~sire.mol.SelectorMol.load_frame` functions. Or, you can
  manually make molecules whole by calling
  :func:`~sire.system.System.make_whole` on :class:`~sire.system.System`,
  or ``mol = mol.move().make_whole().commit()`` on any view.

* Significantly accelerated the reading and writing of files, especially Amber
  topology files.

* Enhanced the integration with NGLView by making it easy to choose colours
  and opacities of representations (e.g. see the :doc:`detailed guide <cheatsheet/view>`).

* Various compile fixes so that :mod:`sire` compiles and works well
  with GCC 13.

* Lots of bug fixes, including `fix_67 <https://github.com/OpenBioSim/sire/issues/67>`__,
  `fix_49 <https://github.com/OpenBioSim/sire/issues/49>`__ (Triclinic box angles
  flipping during a trajectory), and `fix_44 <https://github.com/OpenBioSim/sire/issues/44>`__.

`2023.2.3 <https://github.com/openbiosim/sire/compare/2023.2.2...2023.2.3>`__ - May 2023
----------------------------------------------------------------------------------------

* Fixed numerical precision issues caused by lattice reduction of triclinic
  lattice box vectors to prevent oscillation of the box angles. This is caused
  by the fixed-width format for box dimensions and angles used in the molecular
  input files. `PR 51 - fix_49_50 <https://github.com/OpenBioSim/sire/pull/51>`__

* Added a ``run_constrained`` entry for the optional ``rdkit`` dependency in our
  conda recipe using a minor level pin. This ensures that the correct version of
  ``rdkit`` is installed alongside ``sire``, i.e. one that is compatible with the
  version that ``sire`` was built against. `PR 51 - fix_49_50 <https://github.com/OpenBioSim/sire/pull/51>`__

`2023.2.2 <https://github.com/openbiosim/sire/compare/2023.2.1...2023.2.2>`__ - April 2023
------------------------------------------------------------------------------------------

* Fixed random crashes when loading Amber PRMTOP files when parallelisation
  was enabled. `PR 45 - fix_44 <https://github.com/OpenBioSim/sire/pull/45>`__

* Fixed failure to read an Amber PRMTOP file when no atom names or residues names
  are set. `PR 43 - fix_42 <https://github.com/OpenBioSim/sire/pull/43>`__

* Edited GitHub Actions workflow so that builds of ``devel`` automatically
  upload to the ``dev`` channel, while builds of ``main`` automatically
  upload to the ``test`` channel (for testing before being re-labelled
  to the ``main`` channel)

`2023.2.1 <https://github.com/openbiosim/sire/compare/2023.2.0...2023.2.1>`__ - April 2023
------------------------------------------------------------------------------------------

* Added in ``openmmtools`` as a host requirement. This allows it to be installed in the
  same environment as :mod:`sire`. Note that this changes the dependencies of :mod:`sire`
  to use an older version of ``libnetcdf``. `PR 34 <https://github.com/OpenBioSim/sire/pull/34>`__

* Reactivated the parallel processing code in the Amber parameter/topology parser.
  This significantly speeds up reading and writing of Amber parameter/topology files.

* Fixed compile issues with some MacOS compilers using the C++ 2017 standard, when
  ``std::unary_function`` has been removed.

* Fixed the lookup of Gromacs wildcard dihedrals of the form ``A-*-*-D``.

* Added full support for Urey-Bradley terms in the Gromacs topology parser.

* Added full support for harmonic improper angles in the Gromacs topology parser.
  Note that we don't yet have support for these in the molecular mechanic engine
  or the openmm converter, so they can only currently be read and written.

* Added a developer check for when the version number has changed, so that
  people compiling manually know when they have to rebuild from scratch.


`2023.2.0 <https://github.com/openbiosim/sire/compare/2023.1.3...2023.2.0>`__ - March 2023
------------------------------------------------------------------------------------------

* Completed the :mod:`sire.convert` framework for interconverting :mod:`sire`
  objects with `BioSimSpace <https://biosimspace.openbiosim.org>`__,
  `RDKit <https://rdkit.org>`__ and `OpenMM <https://openmm.org>`__.
  This is now :doc:`fully documented in a tutorial <tutorial/part05/01_convert>`.

* Added support for creating molecules from smiles strings, or generating
  smiles strings from molecules, based on the RDKit integration. Have
  also added a :func:`~sire.mol.SelectorMol.view2d` function that generates
  two-dimensional structure views of molecules. These have infered bond orders,
  formal charges and stereochemistries. This is documented in
  :doc:`two <tutorial/part05/02_view>` :doc:`tutorials <tutorial/part05/03_smiles>`.

* Added new support to the 3D view code to give control over the representation
  used to view the molecule (e.g. licorice, spacefill, cartoon etc). This is
  documented in full (together with more detail about 2D views) in
  a :doc:`detailed guide <cheatsheet/view>`.

* Added support for performing minimisation and molecular dynamics simulations
  based on the OpenMM integration. This is documented in full via both
  :doc:`a tutorial <tutorial/part05/05_dynamics>` and a
  :doc:`detailed guide <cheatsheet/openmm>`.

* Fixed the Amber PRMTOP `dihedral ring bug <https://github.com/OpenBioSim/sire/commit/397271f4229f3cbed6a4c3b425e4baaf4aae4ec5>`__.

* Fixed the bug regarding preservation of water properties when
  `changing topology <https://github.com/michellab/BioSimSpace/issues/247>`__.

* Fixed the bug that caused simulation restarts from short ``waterswap``
  jobs `to fail <https://github.com/OpenBioSim/sire/issues/11>`__.

* Added versioned package support to :func:`sire.utils.try_import`. Now the version
  of the package to be installed can be specified.

* Moved ``pymbar`` from a ``run`` to ``host`` dependency, and switched
  ``analyse_freenrg`` to use :func:`~sire.utils.try_import` to import
  the module. :mod:`sire` now doesn't depend on ``pymbar<4``. Instead,
  ``pymbar`` will be installed at run-time if ``analyse_freenrg`` is
  used in ``mbar`` mode.

* Updated the list of build, run and host dependencies to reduce the number
  of pinned dependencies for :mod:`sire`. This included fixing the way we
  specify ``blas`` so that we don't force a pin to ``openblas``,
  removing the requirement for ``watchdog`` as it is not used any more,
  removing ``pypdb`` from the BioSimSpace run requirements,
  and switching to ``qt-main`` rather than the entire ``qt`` package. Our run
  dependencies are now just ``boost``, ``gsl``, ``lazy_import``,
  ``libnetcdf``, ``openmm``, ``pandas``, ``qt-main``, ``rich`` and ``tbb``.

* Updated the name of the `TIP4P template <https://github.com/OpenBioSim/sire/commit/60cb5827635de0abc7f88419b596586c0e8c185f>`__
  to match convention.

* Added a utility function used by BioSimSpace to remove specified named
  properties from all molecules in a collection.

* Fixed `the bug in the Gro87 parser <https://github.com/OpenBioSim/sire/issues/21>`__
  whereby garbage velocities were written for molecules that didn't have
  a velocity property. These will now be given a default velocity of zero.

* Added an option that can be used to fix an
  `atom numbering issue <https://github.com/OpenBioSim/sire/issues/23>`__ when
  writing PDB files that involve ``TER`` records and multiple molecules.

* Added a fix to `replace spaces <https://github.com/OpenBioSim/sire/commit/6cb7df19721799ff771f235606350bba96bd6e4b>`__
  in GROMACS molecular topology names with underscores, so that topology files
  written by :mod:`sire` can be read by GROMACS.

* Added the :class:`sire.system.ForceFieldInfo` class to hold and report
  metadata related to the forcefields used to calculate energies and
  perform molecular dynamics. This is now used to parse and interpret
  this metadata, giving consistency between the new OpenMM-based
  dynamics code and the energy functions that used the
  in-built molecular mechanics engine.

* Added `a fix <https://github.com/OpenBioSim/sire/commit/71fcf9a0345f9e07b3ec9f56fe4f33b1aada6d4b>`__
  for better handling of :class:`~sire.mol.AtomRadii`-based properties.  This
  helps ensure that radii will be given lengths by default, even if they
  are initialised with zero values.

* Removed the global warnings filter as this was no longer needed.
  :mod:`sire` will now not automatically filter out all warnings.

* Updated :class:`~sire.utils.Console` to use the in-built spinner from
  `rich <https://rich.readthedocs.io>`__ rather than one based on ``yaspin``.
  This removes a dependency and also better integrates the spinner code.

* Added Python 3.10 support and now build Python 3.10 packages. This is now
  the default version of Python for :mod:`sire`, and the version we
  recommend for new workflows. Note that we will drop automatic building
  of Python 3.8 packages later this year (likely Q3 or Q4). This will be
  timed to co-incide with when we add Python 3.11 support, and when
  (we anticipate) conda-forge will drop Python 3.8. Our aim is to only
  build packages for a maximum of 3 Python versions at a time.

* Added the ``future`` branch for feature branches that are accepted,
  but not yet ready for the next release. Adopting a more
  :doc:`regular release and bugfix process <contributing/roadmap>`
  based on a quarterly release cycle.


`2023.1.3 <https://github.com/openbiosim/sire/compare/2023.1.2...2023.1.3>`__ - February 2023
---------------------------------------------------------------------------------------------

* Added the beginnings of the new :mod:`sire.convert` framework for converting
  between different molecule object formats. Created initial converters for RDKit,
  so that we can convert sire molecules to RDKit molecules. This is still considered
  experimental. It will be cleaned up fully for 2023.2.0. It has been added now
  to let others play with this code, to refine a workable API.

* Used the RDKit code to create a :func:`sire.smiles` function to create molecules
  from smiles strings. This is still considered experimental. It will be cleaned
  up fully for 2023.2.0. It has been added to let others begin to explore
  how this capability could be useful.

* Used the RDKit code to create a :func:`~sire.mol.SelectorMol.view2d` function for
  quickly creating 2D views of molecules (or all molecules in a container / system).
  Again, this is considered experimental. It will be cleaned up fully for 2023.2.0.
  It has been added to let others beging to explore how this capability could be
  useful.

* Fixed the SDF bug reported in `issue #8 <https://github.com/OpenBioSim/sire/issues/8>`__.

* Fixed a bug in writing Amber PRMTOP files, where atoms with index zero should not
  be written to the third or fourth column of dihedral / improper entries.

* Adjusted the cutoffs and schemes so that the `.energy()` function gives energies
  that closely agree with those reported by pmemd. Added a unit test that validates
  this.

* Added an :func:`~sire.mol.MoleculeView.extract` function so that it is easy
  to create a new molecule as a subset of another molecule (and the same for
  molecule containers)

* Switched fully to need a C++ 2017 compiler, and adapted the code to fully
  support C++ 2017. Added guards to reduce the number of spurious compiler
  warnings emitted by dependencies of sire during a compile.

* Fixed bugs related to null space parameters specified for triclinic spaces.

* Added classes at the C++ level to represent Stereochemistry, Hybridization
  Chirality, and BondOrder. These are used by the RDKit code and the SDF parser.
  These will be fully exposed in a later release.

`2023.1.2 <https://github.com/openbiosim/sire/compare/2023.1.1...2023.1.2>`__ - February 2023
---------------------------------------------------------------------------------------------

* Used clang-format to autoformat all the C++ files.
* Fixed SDF pickle bug (molecules read from SDF files could not be pickled / unpickled)
* Fixed the bugs in waterswap that led to incorrect energies being calculated.
* Fixed bugs in analyse_freenrg that prevented it from running on newly generated simfiles.
* Fixed a segfault when searching for non-existant atoms in a molecule editor.

`2023.1.1 <https://github.com/openbiosim/sire/compare/2023.1.0...2023.1.1>`__ - January 2023
--------------------------------------------------------------------------------------------

* Fix incompatibility between the updated code and the Boresch restraint code.
* Fixes try_import so that it works within a conda environment, and so that
  it only uses ``conda`` or ``mamba`` to install dependencies.
* Fixed ``NaN`` values of ``r0`` for null amber bonds and angles. Now the
  value of ``r0`` is taken from the current bond length, or else the
  options ``keep_null_bonds=False`` or ``keep_null_angles=False`` can be
  passed via a ``map`` to prevent the writing of null bonds and angles
  to amber parameter files.
* Fixed a bug in :func:`sire.save` that meant that the save directory was
  ignored when the format was specified. Files will now save into the correct
  directory.
* Updated the instructions for :doc:`writing unit tests <contributing/development>`
  to say how to use fixtures to load files, and how to use ``tmpdir`` to write
  files to a temporary directory during a test.
* Addition of lots of files, e.g. issue templates, pull request templates,
  security file etc to improve community engagement via GitHub.
* Created `sire_bigtests <https://github.com/openbiosim/sire_bigtests>`__ from
  `SireUnitTests <https://github.com/michellab/SireUnitTests>`__ and created
  an integration testing pipeline based on these tests. Now the latest ``devel``
  release can be tested via `sire_bigtests <https://github.com/openbiosim/sire_bigtests>`__
  as an extra validation check before creating a release. This release has
  been checked this way :-)
* Lots of minor bugfixes related to those checks, e.g. mostly relating
  to fixing paths on Windows. Now all the integration tests pass on Windows
  (something not before attained, as running the tests on Windows was
  not easy).

`2023.1.0 <https://github.com/openbiosim/sire/releases/tag/2023.1.0>`__ - January 2023
--------------------------------------------------------------------------------------

* Initial release of the OpenBioSim version of sire. The code has been completely
  refurbished using a tutorial-driven development process and has a new
  public API. This is now :mod:`sire`, rather than ``Sire``. The new
  API is activated when you import from this module. You can still use the
  old API by calling :func:`sire.use_old_api` or :func:`sire.use_mixed_api`.
  The new API is pythonic in style, with our aim to be fully PEP8 compliant.
  Functions are named in snake_case, with classes in CapitalCase. Modules
  are all in lowercase. Only a portion of the legacy Sire API has been
  exposed publicly. You can access unexposed classes / functions via
  ``sire.legacy.Module``, e.g. ``sire.legacy.Mol.Connectivity`` will
  get access to the ``Sire.Mol.Connectivity`` class.

* We have
  a `new website <https://sire.openbiosim.org>`__ with easy
  `install instructions <https://sire.openbiosim.org/install>`__, a
  `quickstart guide <https://sire.openbiosim.org/quickstart>`__ and
  a `comprehensive tutorial <https://sire.openbiosim.org/tutorial>`__.
  This is built using sphinx from the files in the ``doc`` directory.

* Migrated from `michellab/sire <https://github.com/michelllab/sire>`__
  to `openbiosim/sire <https://github.com/openbiosim/sire>`__. The new
  repo has had old large files removed, and so is much smaller,
  and so quicker and easier to clone.

* Added a :func:`~sire.mol.SelectorMol.find` function to all of the
  molecule view containers. This returns the index of the view(s)
  within the container. This can be used to quickly get the index
  of, e.g. atoms in a system via ``mols.atoms().find(atom)``.

* Made sure that all units and constants were exposed to the
  new public API, and that the constants were exposed with units, e.g.
  now ``sire.units.k_boltz * (25 * sire.units.celsius)`` gives
  ``0.592486 kcal mol-1`` (be careful to put brackets around the
  temperature, or it will be ``25*k_boltz`` multiplied by ``1 celsius``).

* Made sure that the Rich console is initialised at module import
  time if the new API is used.

* Moved ``show_warnings`` to default ``True`` when loading files. This
  now prints out the method to silence warnings. This is better for, e.g.
  loading gromacs topologies, which were too noisy when ``show_warning``
  was ``False`` and a message told you how to turn them on...

* Added `sse2neon <https://github.com/DLTcollab/sse2neon>`__ so
  that we can use the manually vectorised code
  on ARM64 systems. This fixed issues with Linux/ARM64. This is as fast,
  if not faster, than relying on openmp::simd as we did before.

* Cleaned up the new sire API
  via use of `__all__` in all of the new modules. The public API is
  very limited at the moment, but will grow as we port in more classes.
  However, the aim is that users will mostly not create classes directly,
  but will instead implicitly create them as they load molecular systems
  and call functions on those systems.

* Fully updated the search functionality, making it more robust, more consistent
  and more powerful. Added a
  `detailed guide <https://sire.openbiosim.org/cheatsheet/search.html>`__
  on the search grammar to the new website.

* Added a set of :class:`~sire.mol.Cursor` classes for editing, and made these
  work consistently with most of the property types. Getting and
  setting properties should now be easier, with auto-wrapping and
  expanding of properties.

* Made the AtomProperty classes behave more like standard python containers.
  This makes them easier to work with, and is the first step to hiding
  them completely (they will eventually be auto-converted to/from standard
  Python containers or NumPy arrays).

* Added :func:`~sire.mol.SelectorMol.apply` and
  :func:`~sire.mol.SelectorMol.apply_reduce` functions that let you map
  functions across all objects in a molecular container.

* Cleaned up the handling of units - now everything maps into
  :class:`~sire.units.GeneralUnit` and
  :class:`~sire.units.GeneralUnitProperty`, which are auto-converted when
  exposed to Python.
  Added Python wrapping and monkey-patching to
  :class:`~sire.maths.Vector` so that it
  has length units. Improved the printing of units to the screen (using
  the correct unicode). Added functions that empower the user to choose
  their own default units, e.g. changing angstroms to picometers, or
  switching to full SI units. This only impacts the Python layer when
  rendering the unit, or auto-converting numbers to units, so does
  not break or change the C++ layer. Any view can now be assigned a
  :class:`~sire.units.GeneralUnit` property.

* Added :class:`~sire.mm.Bond`, :class:`~sire.mm.Angle`,
  :class:`~sire.mm.Dihedral`, :class:`~sire.mm.Improper` and their related
  molecule view container
  classes (e.g. :class:`~sire.mm.SelectorBond`,
  :class:`~sire.mm.SelectorMBond` etc). This allows you to have
  molecule views that represent bonds, angles and dihedrals (or collections
  of these). Added measurement functions so that you can easily get their
  lengths or sizes.

* Added :func:`~sire.mol.SelectorMol.energy` to let you calculate
  energies of views (or views with views).
  This uses the parameters / forcefield loaded with the molecule(s). You can
  get energies of any views of sub-views. Also created an proper return type
  for energies that embeds the energy components. Now
  ``view.energy().components()``
  works as you would expect.

* Added :func:`~sire.mol.SelectorMol.energies` to molecule containers so that
  you can get the energies
  of each view in the container. Added support for progress bars using Rich so
  that the user has an indication of progress.

* Added initial support for trajectories. Reworked the molecular parser so that
  multiple "frame" types files will load multiple frames of a trajectory
  (e.g. so that a trajectory can be loaded from multiple PDB files, or
  from multiple DCD or traj files). Added a
  :class:`~sire.mol.TrajectoryIterator` class that
  lets you easily iterate over and query trajectories. Fully documented this
  in the tutorial. You can now do cool things like measure bonds over
  trajectories, or evaluate energies.

* Added a :func:`~sire.mol.SelectorMol.view` function based on
  NGLView that lets you easily see any
  molecule view (or collection of molecule views). Added a
  :func:`~sire.save_to_string` function
  that writes a text-based molecule file to an in-memory string rather than
  a file (so that you don't have to use temporary files with NGLView). Added
  support for viewing trajectories, so that trajectories that are loaded in
  sire are also playable in NGLView.

* Added movement functions to the Cursor classes so that you can more easily
  move molecules (or molecule views). Documented this in
  :doc:`the tutorial <tutorial/part04/05_movement>`. Re-worked
  the way PropertyMap is passed via Python. Now have a
  :func:`~sire.base.create_map` function
  that can create a PropertyMap from anything that is passed. This has some
  examples in its documentation that show how is can be used. Made sure that
  all of the new functionality can use PropertyMap and uses
  :func:`~sire.base.create_map`
  to support function calls like
  ``cursor.translate( (1,2,3), map={"coordinates":"coords2"} )``.

* Speaking of which, also updated :class:`~sire.maths.Vector`
  adding in functions that
  allow auto-conversion of list-like python objects to
  :class:`~sire.maths.Vector`.
  It should almost be the case that a user will not have to use this class
  directly themselves, as things should just auto-convert. Added support for
  creating Vectors from plain numbers or length units, using the default length
  unit if plain numbers are used.

* Removed lots of unnecessary files. Moved some files into the website docs
  so that there is a single source of truth. Updated paths
  and links to point to the new locations in OpenBioSim. Fixed CI build issues
  on Windows by building in the right directory. Updated the pythonizing framework
  so that we only pythonize the C++ layer, and avoid the circular dependencies
  that were causing random import errors (particularly on Windows).

* Fixed lots of bugs and expanded the unit test suite to test the above
  functionality.

GitHub (michellab): June 22nd 2015 - January 2023
-------------------------------------------------

Thanks to `@ppxasjsm <https://github.com/ppxasjsm>`__ and
`@jmichel80 <https://github.com/jmichel80>`__ development
was migrated into the `michellab <https://github.com/michellab>`__
organisation on `GitHub <https://github.com/michellab/sire>`__.

This comprised 2495 commits, from developers
`@lohedges <https://github.com/lohedges>`__,
`@chryswoods <https://github.com/chryswoods>`__,
`@ppxasjsm <https://github.com/ppxasjsm>`__,
`@halx <https://github.com/halx>`__,
`@jmichel80 <https://github.com/jmichel80>`__,
`@ptosco <https://github.com/ptosco>`__,
`@SofiaBariami <https://github.com/SofiaBariami>`__,
`@fjclark <https://github.com/fjclark>`__,
`@Steboss <https://github.com/Steboss>`__,
`@nigel-palmer <https://github.com/nigel-palmer>`__,
`@msuruzon <https://github.com/msuruzhon>`__ and
`@kexul <https://github.com/kexul>`__.

Here is the changelog for this stage of development.

..

    [2023.0.3] January 2023: Added the beginnings of a new sphinx-based website
               (in the `doc` folder), which includes the sire API documentation.
               Cleaned up the new sire API via use of `__all__` in all of the
               new modules. The public API is very limited at the moment, but
               will grow as we port in more classes.  However, the aim is that
               users will mostly not create classes directly, but will instead
               implicitly create them as they load molecular systems and call
               functions on those systems. Added a tutorial to this website
               that will be used to demonstrate and teach the new sire API.
               Fully updated the search functionality, making it more robust,
               more consistent and more powerful. Added a detailed guide on the
               search grammar to the new website. Added a set of Cursor classes
               for editing, and made these work consistently with most of the
               property types. Getting and setting properties should now be
               easier, with auto-wrapping and expanding of properties. Made
               the AtomProperty classes behave more like standard python
               containers.  This makes them easier to work with, and is the
               first step to hiding them completely (they will eventually be
               auto-converted to/from standard Python containers or NumPy
               arrays. Added `apply` and `apply_reduce` functions that let you
               map functions across all objects in a molecular container.
               Cleaned up the handling of units - now everything maps into
               GeneralUnit and GeneralUnitProperty, which are auto-converted
               when exposed to Python. Added Python wrapping and
               monkey-patching to sire.maths.Vector so that it has length units.
               Improved the printing of units to the screen (using the correct
               unicode). Added functions that empower the userto choose their
               own default units, e.g. changing angstroms to picometers, or
               switching to full SI units. This only impacts the Python layer
               when rendering the unit, or auto-converting numbers to units,
               so does not break or change the C++ layer. Any view can now be
               assigned a GeneralUnit property. Added Bond, Angle, Dihedral,
               Improper and their related molecule view container classes (e.g.
               SelectorBond, SelectorMBond etc). This allows you to have
               molecule views that represent bonds, angles and dihedrals (or
               collections of these). Added measurement functions so that you
               can easily get their lengths or sizes. Added `.energy()` to let
               you calculate energies of views (or views with views). This uses
               the parameters / forcefield loaded with the molecule(s). You can
               get energies of any views of sub-views. Also created an proper
               return type for energies that embeds the energy components.
               Now `view.energy().components()` works as you would expect.
               Added `.energies()` to molecule containers so that you can get
               the energies of each view in the container. Added support for
               progress bars using Rich so that the user has an indication of
               progress. Added initial support for trajectories. Reworked the
               molecular parser so that multiple "frame" types files will load
               multiple frames of a trajectory (e.g. so that a trajectory can
               be loaded from multiple PDB files, or from multiple DCD or traj
               files). Added a TrajectoryIterator class that lets you easily
               iterate over and query trajectories. Fully documented this in
               the tutorial. You can now do cool things like measure bonds over
               trajectories, or evaluate energies. Added a `.view()` function
               based on NGLView that lets you easily see any molecule view (or
               collection of molecule views). Added a `save_to_string` function
               that writes a text-based molecule file to an in-memory string
               rather than a file (so that you don't have to use temporary
               files with NGLView). Added support for viewing trajectories, so
               that trajectories that are loaded in sire are also playable in
               NGLView. Added movement functions to the Cursor classes so that
               you can more easily move molecules (or molecule views).
               Documented this in the tutorial. Re-worked the way PropertyMap is
               passed via Python. Now have a sire.base.create_map function that
               can create a PropertyMap from anything that is passed. This has
               some examples in its documentation that show how is can be used.
               Made sure that all of the new functionality can use PropertyMap
               and uses `create_map` to support function calls like
               `cursor.translate( (1,2,3), map={"coordinates":"coords2"} )`.
               Speaking of which, also updated `sire.maths.Vector` adding in
               functions that allow auto-conversion of list-like python objects
               to `sire.maths.Vector`. It should almost be the case that a user
               will not have to use this class directly themselves, as things
               should just auto-convert. Added support for creating Vectors
               from plain numbers or length units, using the default length
               unit if plain numbers are used. Fixed lots of bugs and expanded
               the unit test suite to test the above functionality. Removed
               lots of unnecessary files. Moved some files into the website
               docs so that there is a single source of truth. Began the process
               of updating paths and links to point to the new locations in
               OpenBioSim. Fixed CI build issues on Windows by building in the
               right directory. Updated the pythonizing framework so that we
               only pythonize the C++ layer, and avoid the circular dependencies
               that were causing random import errors (particularly on Windows).

    [2023.0.2] December 2022: Fix multiple distance restraint bug in SOMD
               (@fjclark). Add support for PME FEP with SOMD and fix
               associated bugs (@halx, @jmichel80). Fix CI issues so that
               PRs use the correct URL when triggered by external forks.
               Exclude dummy atoms when repartitioning hydrogen masses.
               Deprecate py37.

    [2023.0.1] November 2022: Improve handling of HETATM and TER records in
               PDB files. Fix SOMD selection issues following update to the
               2023 API. Fix writing of steps to SOMD simfile.dat (@fjclark).
               Throw exception when CHAMBER format AMBER topology files are
               detected. Expose toVector() method for the velocity property.
               Match against inverted dihedral records of for A-B-C-A when
               building GROMACS topologies. Fixed calling of static Py++
               functions. Build against conda-forge AmberTools and GROMACS
               packages as host requirements, allowing users to create
               BioSimSpace environments with or without these dependencies
               installed. Added the ability to search on whether or not a
               property exists.  Make sure searches are returned in MolIdx
               order. Ensure Sire is built against packages with the "dev"
               label.

    [2023.0.0] July 2022 - Updated Sire's API to a more pythonic style.
               Module names are in lower case, e.g. `import Sire` becomes
               `import sire`, or `import sire as sr`. Functions are in
               underscore_case. This change is not backwards compatible. To
               support old code, a `sire.use_old_api()` function has been added.
               New functions have been added that make it easier to load
               and save molecules. These can load from URLs. Tests have been
               updated to pytest and now load input data from the sire website.
               The search system has been overhauled, optimised and updated.
               This is described in the new tutorials that are in the process
               of being written in the `doc` directory. This also contains
               the new sphinx website. The `CMakeLists.txt` files and build
               system have been completely reworked. These now use more
               pythonic `setup.py` scripts. These have been updated to fully
               support MacOS M1 and Windows. The conda recipe has been
               updated to use these scripts. Conda packages are now built
               and supported across Linux, MacOS and Windows.

    [2022.3.0] June 2022 - Added support for parsing SDF files (@chryswoods).
               Move conda build process to Miniforge and mambdabuild (boa) to
               avoid timeouts and memory issues. Update GroTop parser to ensure
               new atom types are created when names match but parameters
               differ. Added additional BioSimSpace wrapper to update
               coordinates and velocities in a system, without first requiring
               that it is modified to have unique atom and residue numbers.
               Use -Oz compiler flag rather than -Os for compiling Python
               wrappers to avoid "illegal hardware instruction" error with
               Clang 14 on macOS x86_64. Fixed issue reconstructing triclinic
               box objects from a binary data stream. Added missing streaming
               operators to Sire.Unit.GeneralUnit.

    [2022.2.0] March 2022 - Fixed formatting of SOLVENT_POINTERS flag in
               AmberPrm7 parser. Removed duplicate definition of sigma_av
               in OpenMMFreEnergySt.cpp. Fixed SOMD issues related to
               assumption that perturbable molecule always has MolIdx(1)
               (@fjclark). Fixed wrappers and added significant performance
               enhancements to the SireIO::updateCoordinatesAndVelocities
               function. This significantly (200x) speeds up the remapping
               of coordinates/velocities from SOMD trajectory frames, which
               was a bottleneck for large protein-ligand simulations within
               BioSimSpace. Disabled GSL error handling to avoid a potential
               segmentation fault within a singular value decomposition
               routine called by SireMaths::align.

    [2022.1.0] Jan 2022 - Fixed counting of protons to account for dummy atoms
               when swapping water topology and ensure that original molecular
               properties are preserved. Added a fallback to the BGFS solver
               to improve robustness of FEP analysis (@kexul). Fixed a bug
               that caused distance restraints to be skipped if the ligand
               wasn't the first molecule in the FEP topology (@jmichel80,
               @fjclark). Improved atomic element inference in AMBER parsers.
               Update Sire build to latest versions of dependencies on macOS
               and Linux. This required substantial mini-changes across the
               entire codebase due to changes in APIs and deprecations. This
               includes moving away from qAlgorithm, using the new Qt
               container constructors, moving to OneAPI, switching to the
               conda-forge OpenMM and switching to the new C++ ABI (@chryswoods).
               Simplified Sire wrapper generation using a minimal Docker
               container with the latest Py++ (@chryswoods). Add support for
               native Python pickling of Sire objects (@chryswoods.) Switch
               to GitHub actions for CI. This uses a conda-forge compliant
               conda build, with packages then uploaded to the Anaconda cloud.

    [2021.1.0] Aug 2021 - Added support for multiple combining rules in SOMD
               (@SofiaBariami). Added support for triclinic simulation boxes.
               Convert Ryckaert-Bellememans form dihedral functions from
               GROMACS to Fourier series to allow conversion to AMBER format.
               Updated search functionality to enable searching for objects
               within an arbitrary distance of a point. Fixed PDB2 parser bug
               to ensure that residue names are fixed width. Ensure that
               NUMEXTRA pointer is written so that AMBER topology files can
               be read by tools such as ParmEd. Write NATYP pointer and
               correct number of SOLTY flags. Even though these aren't used,
               incorrect values break external tools, e.g. ParmEd. Added
               support for AMBER TIP5P water topology conversion. Correctly
               flag OPLS style force fields when creating MMDetail object so
               that users can reconstruct OPLS systems written to AMBER format.
               Made build Python 3.8+ compliant. (Python libraries are now
               ABI compatible.) Switched to using std::atomic since
               tbb/atomic.h is now deprecated. Switched to using HTTPS for
               sending analytics. Updated build to be able to link against
               conda version of libcpuid. Added support for generating PDB
               CONECT records from a Sire.Mol.Connectivity object. Fixed
               issue with PMEMD skipping torsions with zero periodicity.
               Fixed random number seeding bug in somd-freenrg, which
               resulted in the OpenMM generator being seeded with the same
               seed for each cycle of the simulation.

    [2020.1.0] July 2020 - Fixed bug in WaterView program to ensure that a
               molecule is extracted from the returned list. Stable sorting
               of dihedrals and other potential terms to allow reproducible
               writing of input files for SOMD (@ptosco). Updated the
               FreeEnergyAnalysis script to support different versions of the
               pymbar API. Significant performance improvement to the GroTop
               parser by looping over cut-groups during non-bonded matrix
               evaluation. Updated Miniconda and conda dependencies to latest
               cross-compatible versions. Fixed minor copiler and runtime
               issues (@nigel-cresset).

    [2019.3.0] November 2019 - Added functionality to restrict the search space
               when finding paths between atoms or searching for rings. Fixed
               performance issue in GroTop parser caused by an N^2 loop over
               atoms when searching the intrascale matrix. We now loop over
               cut-groups, which is far more efficient. Fixed issues with
               Python wrapper generation caused by issues with missing define
               symbols and a bug in the scanheaders.py script.

    [2019.2.1] October 2019 - Updated the Conda recipe to pin the dependencies
               of dependencies that are used at run time since Conda doesn't
               automatically do this for you. Added instructions detailing the
               Azure Pipeline build process and how to create a new release.

    [2019.2.0] September 2019 - Updated the Gromacs topology writer to support
               perturbable molecules containing a variable number of bonds.
               Created a Docker container for building wrappers and updated
               to using CastXML. Added support for running background
               processes on Windows (@ptosco). Updated SOMD Python wrapper
               to write restart files every cycle to simplify system monitoring
               in BioSimSpace. Fixed macOS build issue by not linking against
               libpython. Made sure that Conda dependencies are pinned
               correctly to avoid compatibility issues. Fixed bug that
               prevented upload statistics being sent and added support for
               tracking BioSimSpace usage.

    [2019.1.0] May 2019 - Updates to the Gromacs topology writer to support
               free energy perturbation simulations. The MCS matching
               functionality has been extended to allow matches between heavy
               and light atoms, and the ability to return all current matches,
               rather than just the most recent. Temporarily disabled
               parallelisation in the AmberPrm parser to avoid a threading
               issue. Switched to using Azure Pipelines for continuous
               integration to enable a fully automated build, testing, and
               deployment pipeline. In addition, we finally have created
               a Sire Conda package to simplify the installation and update
               process.

    [2018.2.0] July 2018 - Improvements to the Gromacs topology reader/writer,
               addition of code to improve matching of atoms in proteins,
               fixing compile issues on modern Ubuntu, bugfixes for crashes
               in the AmberPrm reader, added in text-based searching for
               atoms, residues etc. from Systems, MoleculeGroups, Molecules,
               etc. based on boost::spirit, updated boost to latest version,
               bugfixes for quantomm infinite rotation bug for ions,
               general bugfixes.

    [2018.1.1] May 2018 - Small bug fixes to allow single-atom solutes
               and also to fix small issues with some parsers for BioSimSpace

    [2018.1.0] March 2018 - Signficantly improved the Gromacs and Charmm
               parser and  fixed bugs. Can now write with both :-). Fixed
               compilation on Windows 7 and above. Small changes to
               the API of AtomProperties to make them easier to work
               with from Python. Added a script to automatically color
               code swap-based free energies. Fixed localisation
               problems for the PDB writer. Improved mbar analysis
               code. Added code to track forcefield data of a molecue.
               Added code to better manage processes, including redirection
               of standard output and error.

    [2017.3.0] December 2017 - Added new PDB (PDB2), Mol2, Gro87, and CharmmPSF
               parsers, as well as a GroTop parser, all part of the new
               MoleculeParser framework. Updated all of the swap-based
               methods to use this.

               Removed the ViewsOfMol Python wrapper and now have the
               code automatically return the correct python object for
               the molecule (or part of molecule) that is returned from
               the system. This makes simple scripts easier to write.

               General bugfixes and optimisations, including fixing
               bugs with the way that PropertyMap worked, cleaning
               up to/from converters from python objects to automatic
               Property wrappers, and fixing Process so that it can
               redirect to stdout and that isRunning works without the
               user having to call "wait" first!

    [2017.2.0] September 2017 - The MoleculeParser framework has
               been created to support reading and writing of molecules
               in lots of different formats. The first set of formats
               that have been completed are Amber Prmtop, Amber Rst7
               and Amber Rst/Trj. The parsers work in parallel, with
               file formats automatically detected by the parser,
               e.g. system = MoleculeParser.read( "file.prm", "file.rst" )
               will automatically do the right thing.

               Improved automatic compilation on Arch linux.

               Fixed temperature checking and general bugfixing for mbar code.

    [2017.1.0] April 20 2017 - First 2017 release. Included new
            parallel MoleculeParser code for reading molecules,
            and moved fully over to C++-14 style coding for new code.
            Included AVX-512 vectorisation for Intel KNL and can
            now successfully compile and run using GCC 5 and GCC 6,
            as well as Intel 2017 compilers and Clang.

    [2016.3.1] January 9 2017 - Minor patch release that fixes bugs:
        (1) Correctly sets MACOSX_DEPLOYMENT_TARGET to 10.8 so executable works
            on OS X 10.8 (Mountain Lion) and above
        (2) Fixed a parsing bug in Parameter that prevented integer or float
            parameters from being passed to ligandswap, waterswap etc.
        (3) Fixed a small bug in MultiDouble that meant it lost precision when
            swapping individual values
        (4) Fixed a parsing bug in Parameter that meant that windows path names
            were not interpreted correctly
        (5) Fixed the build scripts so that they placed bundled libraries into
            bundled/lib rather than bundled/lib64 (affected SUSE-based distributions)

    [2016.3.0] December 22 2016 - Public release containing full LigandSwap. Uses
     new optimised forcefields for energy calculations, built on top of Intel Threaded
     Building Blocks for parallelisation. New code is significantly faster with better
     scaling.

    [2016.2.0] June 3 2016 - Semi-private release for the CCP-BioSim workshops. Included
     the first version of LigandSwap and general bug fixes

    [2016.1.0] April 29 2016 - Merge of Bristol and Edinburgh codes, moved to miniconda
     and clean packaging system, including OpenMM fully, added in nautilus, somd etc.,
     added optimised forcefields, added a proper unit testing suite.

    [OLD] Updated gradient compuation in openmmfreenergst to finite differece gradients
    Allow the computation of reduced perturbed energies in openmmfreenergst of all computed lambda values
    Separated minimization and equilibration from production run.
    Implemented mass repartitioning for hydrogens atoms to allow for larger integration timesteps
    Added nautilus scripts

Google Code: August 7th 2006 - April 1st 2015
---------------------------------------------

Sire was developed against the subversion repository provided
by Google Code. Here is an
`archive of the repository <https://code.google.com/p/sire>`__.

This comprised 2775 commits, from developers
`@chryswoods <https://github.com/chryswoods>`__,
`@jmichel80 <https://github.com/jmichel80>`__ and
`@nividic73 <mailto:nividic73@googlemail.com>`__.

Local Subversion: February 5th 2005 - July 25th 2006
----------------------------------------------------

Sire was developed against a local subversion repository.
Here is a
`svndump of the original repository <https://sire.openbiosim.org/f/orig_sire_repository.dump.bz2>`__,
and all of the `commit history <https://sire.openbiosim.org/f/original_repository_comments.txt>`__.

This comprised 831 commits from developer `@chryswoods <https://github.com/chryswoods>`__.

Sire started as ``ProtoMS 3``, a complete C++ rewrite of
`ProtoMS 2 <https://code.google.com/archive/p/protoms/source/default/commits>`__,
developed originally as a Fortran program
by `@chryswoods <https://github.com/chryswoods>`__ and
`@jmichel80 <https://github.com/jmichel80>`__. ProtoMS has since continued
to be developed by the
`Essex Group <https://www.essexgroup.soton.ac.uk>`__ and is
itself now available as `ProtoMS 3.4 <https://protoms.org>`__.

More detail about the history and parallel development of Sire and
ProtoMS can be `found here <https://www.essexgroup.soton.ac.uk/ProtoMS/FAQ/index.html>`__.
