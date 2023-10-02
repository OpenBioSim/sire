======================================
Running Faster Free Energy Simulations
======================================

The previous section showed how to calculate the relative hydration free
energy of ethane and methanol using alchemical dynamics. Short dynamics
simulation were run for each λ-value, with the energy differences
between neighbouring λ-values used to calculate the free energy differences.

The simulations used a timestep of 1 fs, and calculated energy differences
every 0.1 ps (100 steps). This calculated a free energy that should be
accurate, but the simulation took a long time to run. In this section, we
will show how to calculate the same result, but with a much faster
simulation.

Timesteps and Constraints
-------------------------

The easiest way to speed up a simulation is to increase the dynamic
timestep. This is the amount of time between each step of the simulation,
i.e. the amount of time between each calculation of the forces. The longer
the timestep, the faster the simulation will run, but at the cost of
more unstable dynamics. At an extreme, the simulation will become
totally unstable and an
`OpenMM Particle Exception <https://github.com/openmm/openmm/wiki/Frequently-Asked-Questions#nan>`__
will be raised.

As a rule of thumb, the timestep should be less than the fastest vibrational
motion in the simulation. Since the fastest vibrations will likely involve
bonds with hydrogen atoms, we can make the simulation more stable by
contraining the bonds that involve hydrogen atoms.

We can choose the constraint to use via the ``constraint`` keyword, e.g.
after loading the molecules and minimising,

>>> import sire as sr
>>> mols = sr.load(sr.expand(sr.tutorial_url, "merged_molecule.s3"))
>>> for mol in mols.molecules("molecule property is_perturbable"):
...     mols.update(mol.perturbation().link_to_reference().commit())
>>> mols = mols.minimisation().run().commit()

...we can turn on constraints of the bonds involving hydrogen atoms by
setting ``constraint`` to ``h-bonds``.

>>> d = mols.dynamics(timestep="2fs", constraint="h-bonds")
>>> d.run("5ps")
>>> print(d)
Dynamics(completed=5 ps, energy=-50251 kcal mol-1, speed=77.8 ns day-1)

has let us run dynamics faster than if we used a 1 fs timestep...

>>> d = mols.dynamics(timestep="1fs", constraint="none")
>>> d.run("5ps")
>>> print(d)
Dynamics(completed=5 ps, energy=-50250.5 kcal mol-1, speed=69.6 ns day-1)

.. note::

   The constraints do incur a cost, hence why the simulation is not
   twice as fast

And this can go a lot faster with a 4 fs timestep...

>>> d = mols.dynamics(timestep="2fs", constraint="h-bonds")
>>> d.run("4ps")
>>> print(d)
Dynamics(completed=4 ps, energy=-50250.8 kcal mol-1, speed=72.1 ns day-1)

However, turning on ``h-bonds`` constraints would not be enough to keep the
simulation stable for larger timesteps. For example, if we use a timestep
of 5 fs...

>>> d = mols.dynamics(timestep="5fs", constraint="h-bonds")
>>> d.run("5ps")
>>> print(d)
OpenMMException: Particle coordinate is NaN.  For more information, see
https://github.com/openmm/openmm/wiki/Frequently-Asked-Questions#nan

For such timesteps, we would need to constrain all bonds and angles
involving hydrogen, using the ``h-bonds-h-angles`` constraint.

>>> d = mols.dynamics(timestep="5fs", constraint="h-bonds-h-angles")
>>> d.run("5ps")
>>> print(d)
Dynamics(completed=5 ps, energy=-50251.5 kcal mol-1, speed=181.5 ns day-1)

You can go even further by constraining all bonds, and all angles involving
hydrogen using the ``bonds-h-angles`` constraint...

>>> d = mols.dynamics(timestep="8fs", constraint="bonds-h-angles")
>>> d.run("5ps")
>>> print(d)
Dynamics(completed=5 ps, energy=-50252.1 kcal mol-1, speed=192.3 ns day-1)

.. note::

   There is also the ``bonds`` constraint which constrains
   all chemical bonds, but this doesn't really help unless
   angles involving hydrogen aren't also constrained.

.. note::

   Different molecular systems will behave differently, and may produce
   more, or less stable dynamics than that shown above. As a rule of thumb,
   start using a 4 fs timestep with ``h-bonds-h-angles`` constraints, and
   then either increase the timestep and dial up the constraints if you can,
   or reduce the timestep and reduce the timestep if the simulation becomes
   unstable.

To make things easy, :mod:`sire` automatically chooses a suitable constraint
based on the timestep. You can see the constraint chosen using the
:meth:`sire.mol.Dynamics.constraint` function.

>>> d = mols.dynamics(timestep="8fs")
>>> print(d.constraint())
bonds-h-angles

You can disable all constraints by setting ``constraint`` equal to ``none``,
e.g.

>>> d = mols.dynamics(timestep="4fs", constraint="none")
>>> print(d.constraint())
none
>>> d.run("5ps")
>>> print(d)
RuntimeError: The kinetic energy has exceeded 1000 kcal mol-1 per atom (it is 2.2202087996265908e+16 kcal mol-1 atom-1, and
2.701328046505673e+20 kcal mol-1 total). This suggests that the simulation has become unstable. Try reducing the timestep and/or minimising
the system and run again.

but do expect to see a ``ParticleException`` or other ``RuntimeError``
exceptions raised at some point!

Constraints and Perturbable Molecules
-------------------------------------

Talk about how to avoid constraining a perturbable bond to the wrong
value. Need to minimise when changing λ (or equilibrating using a
short timestep with no constraints).

Show how to use ``perturbable_constraint`` to choose a different constraint
method for perturbable molecules.

Hydrogen Mass Repartitioning
----------------------------

Bonds involving hydrogen atoms vibrate quickly because vibrational frequency
is related to atomic mass - the lighter the atom, the faster it will
vibrate. We can reduce the frequency of these vibrations by increasing the
mass of the hydrogens. Fortunately, free energy is derived from the
potential energy of the molecules, which is independent of their mass.
So, we are free to magically move mass from heavy atoms such as carbon to
their bonded hydrogen atoms without affecting the free energy.

This method, called "hydrogen mass repartitioning", is implemented in
the :func:`sire.morph.repartition_hydrogen_masses` function. It takes a
molecule as argument, and returns that same molecule with its hydrogen
masses repartitioned.

>>> mol = mols.molecule("molecule property is_perturbable")
>>> repartioned_mol = sr.morph.repartition_hydrogen_masses(mol)
>>> for atom0, atom1 in zip(mol.atoms(), repartioned_mol.atoms()):
...    print(atom0, atom0.property("mass"), atom1.property("mass"))
Atom( C1:1    [  25.71,   24.94,   25.25] ) 12.01 g mol-1 2.938 g mol-1
Atom( C2:2    [  24.29,   25.06,   24.75] ) 12.01 g mol-1 2.938 g mol-1
Atom( H3:3    [  25.91,   23.89,   25.56] ) 1.008 g mol-1 4.032 g mol-1
Atom( H4:4    [  26.43,   25.22,   24.45] ) 1.008 g mol-1 4.032 g mol-1
Atom( H5:5    [  25.86,   25.61,   26.13] ) 1.008 g mol-1 4.032 g mol-1
Atom( H6:6    [  24.14,   24.39,   23.87] ) 1.008 g mol-1 4.032 g mol-1
Atom( H7:7    [  24.09,   26.11,   24.44] ) 1.008 g mol-1 4.032 g mol-1
Atom( H8:8    [  23.57,   24.78,   25.55] ) 1.008 g mol-1 4.032 g mol-1

The repartitioned molecule has the same mass as the original molecule, but
the hydrogens have been made heavier. The mass of the carbon atoms has been
reduced to compensate.

>>> print(mol.mass(), repartioned_mol.mass())
30.068 g mol-1 30.068 g mol-1

It is normal to only repartition the hydrogen masses of perturbable molecules.
This lets us use ``h-bond-h-angles`` constraints for all molecules, with
no constraints on the perturbable molecules. But, we don't need constraints
on the perturbable molecules because their hydrogens are heavier, and so the
vibrations of their atoms should be slower.

>>> mols.update(mol)
>>> d = mols.dynamics(timestep="8fs", constraint="bonds-h-angles")
>>> d.run("5ps")
>>> print(d)
Dynamics(completed=5 ps, energy=-50252.1 kcal mol-1, speed=180.7 ns day-1)
