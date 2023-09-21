===================
Alchemical Dynamics
===================

You can create an alchemical molecular dynamics simulation in exactly the
same way as you would a normal molecular dynamics simulation. There
are two options:

* Use the high-level interface based on the :func:`~sire.system.System.minimisation` and :func:`~sire.system.System.dynamics` functions.
* Use the low-level interface that works directly with native OpenMM objects.

High level interface
--------------------

The simplest route is to use the high-level interface. Calling
:func:`~sire.system.System.minimisation` or
:func:`~sire.system.System.dynamics` on any collection of molecules (or
view) that contains one or more merged molecules will automatically create
an alchemical simulation. There are extra options that you can pass to
the simulation that will control how it is run:

* ``lambda_value`` - this sets the global λ-value for the simulation.
  λ is a parameter that controls the morph from the reference state
  (at λ=0) to the perturbed state (at λ=1).

* ``swap_end_states`` - if set to ``True``, this will swap the end states
  of the perturbation. The morph will run from the perturbed state
  (at λ=0) to the reference state (at λ=1). Note that the coordinates
  of the perturbed molecule will be used in this case to start the
  simulation. This can be useful to calculate the reverse of the free
  energy potential of mean force (PMF), to check that the reverse
  free energy equals to forwards free energy.

* ``ignore_perturbations`` - if set to ``True``, this will ignore any
  perturbations, and will run the dynamics without a λ-coordinate, using
  the properties from the reference state (or from the perturbed state
  if ``swap_end_states`` is ``True``). This is useful if you want to
  run a standard dynamics simulation of the reference or perturbed states,
  without the machinery of alchemical dynamics.

For example, we could minimise our alchemical system at λ=0.5 using

>>> import sire as sr
>>> mols = sr.load(sr.expand(sr.tutorial_url, "merged_molecule.s3"))
>>> mols = mols.minimisation(lambda_value=0.5).run().commit()

We can then run some dynamics at this λ-value using

>>> d = mols.dynamics(lambda_value=0.5, timestep="4fs", temperature="25oC")
>>> d.run("10ps", save_frequency="1ps")
>>> mols = d.commit()
>>> print(d)
Dynamics(completed=10 ps, energy=-31213.8 kcal mol-1, speed=194.7 ns day-1)

The result of dynamics is a trajectory run at this λ-value. You can view the
trajectory as you would any other, using ``mols.view()``.

Energy Trajectories
-------------------

In addition to a coordinates trajectory, dynamics also produces an
energy trajectory. This is the history of kinetic and potential energies
sampled by the molecules during the trajectory. You can access this
energy trajectory via the :func:`~sire.system.System.energy_trajectory`
function, e.g.

>>> df = mols.energy_trajectory()
>>> print(df)
          kinetic     potential
time
1.0   4461.033219 -43523.347617
2.0   6095.724385 -41312.514682
3.0   6588.720597 -40003.958276
4.0   6853.529444 -39283.565112
5.0   7080.551573 -38722.021135
6.0   7073.422652 -38677.147808
7.0   7086.066481 -38527.530217
8.0   7264.413320 -38408.505361
9.0   7271.934837 -38589.671709
10.0  7214.294605 -38444.356221

.. note::

   Dynamics uses a random number generator for the initial velocities
   and temperature control. The exact energies you get from dynamics will
   be different to what is shown here.

By default, this trajectory is returned as a
`pandas DataFrame <https://pandas.pydata.org/pandas-docs/stable/reference/frame.html>`__.
You can get the underlying :class:`sire.maths.EnergyTrajectory` function
by passing ``to_pandas=False``, e.g.

>>> t = mols.energy_trajectory(to_pandas=False)
>>> print(t)
EnergyTrajectory( size=10
time	kinetic	potential
1	4461.03	-43523.3
2	6095.72	-41312.5
3	6588.72	-40004
4	6853.53	-39283.6
5	7080.55	-38722
6	7073.42	-38677.1
7	7086.07	-38527.5
8	7264.41	-38408.5
9	7271.93	-38589.7
10	7214.29	-38444.4
)

You calculate free energies by evaluating the potential energy for different
values of λ during dynamics. You can control which values of λ are used
(the so-called "λ-windows") by setting the ``lambda_windows`` argument, e.g.

>>> d = mols.dynamics(lambda_value=0.5, timestep="4fs", temperature="25oC")
>>> d.run("10ps", save_frequency="1ps", lambda_windows=[0.4, 0.6])
>>> mols = d.commit()
>>> print(mols.energy_trajectory())
               0.4           0.5           0.6      kinetic     potential
time
1.0            NaN           NaN           NaN  4461.033219 -43523.347617
2.0            NaN           NaN           NaN  6095.724385 -41312.514682
3.0            NaN           NaN           NaN  6588.720597 -40003.958276
4.0            NaN           NaN           NaN  6853.529444 -39283.565112
5.0            NaN           NaN           NaN  7080.551573 -38722.021135
6.0            NaN           NaN           NaN  7073.422652 -38677.147808
7.0            NaN           NaN           NaN  7086.066481 -38527.530217
8.0            NaN           NaN           NaN  7264.413320 -38408.505361
9.0            NaN           NaN           NaN  7271.934837 -38589.671709
10.0           NaN           NaN           NaN  7214.294605 -38444.356221
11.0 -40712.392040 -40712.520657 -40712.684409  6154.441842 -40712.520657
12.0 -39714.124832 -39714.253449 -39714.327574  6611.854460 -39714.253449
13.0 -39152.132719 -39152.769223 -39153.291484  6924.754530 -39152.769223
14.0 -38791.263933 -38791.601680 -38791.944687  7027.605561 -38791.601680
15.0 -38568.331333 -38568.519701 -38568.683453  7040.809156 -38568.519701
16.0 -38558.203465 -38558.302206 -38558.376331  7166.843369 -38558.302206
17.0 -38605.586352 -38606.133229 -38606.595738  7195.020226 -38606.133229
18.0 -38618.731667 -38618.412148 -38618.068013  7170.302626 -38618.412148
19.0 -38611.651122 -38612.317502 -38612.989142  7182.904255 -38612.317502

This has run a new trajectory, evaluating the potential energy at the
simulation λ-value (0.5) as well as at λ-windows 0.4 and 0.6. You can pass in
as many or as few λ-windows as you want.

.. note::

   Notice how the potential energy is evaluated at λ=0.4, λ=0.5 and λ=0.6
   only from 11ps onwards. The first 10ps was the first block of dynamics
   where we only evaluted the energy at the simulated λ-value. The second
   block of 10ps also evaluated the energy at λ=0.4 and λ=0.6.

Controlling the trajectory frequency
------------------------------------

The ``save_frequency`` parameter controls the frequency at which both
coordinate frames and potential energies are saved to the trajectory.

Typically you want to evaluate the energies at a much higher frequency than
you want to save frames to the coordinate trajectory. You can choose
a different frequency by either using the ``frame_frequency`` option to
choose a different coordinate frame frequency, and/or using the
``energy_frequency`` option to choose a different energy frequency.

For example, here we will run dynamics saving coordinates every picosecond,
but saving energies every 20 femtoseconds.

>>> d = mols.dynamics(lambda_value=0.5, timestep="4fs", temperature="25oC")
>>> d.run("10ps", frame_frequency="1ps", energy_frequency="20fs",
...       lambda_windows=[0.4, 0.6])
>>> mols = d.commit()
>>> print(mols.energy_trajectory())
                0.4           0.5           0.6      kinetic     potential
time
1.00            NaN           NaN           NaN  4461.033219 -43523.347617
2.00            NaN           NaN           NaN  6095.724385 -41312.514682
3.00            NaN           NaN           NaN  6588.720597 -40003.958276
4.00            NaN           NaN           NaN  6853.529444 -39283.565112
5.00            NaN           NaN           NaN  7080.551573 -38722.021135
...             ...           ...           ...          ...           ...
29.92 -38389.435539 -38389.414777 -38389.279773  7332.296186 -38389.414777
29.94 -38353.913312 -38354.131556 -38354.235557  7301.874742 -38354.131556
29.96 -38371.599736 -38371.549099 -38371.384218  7293.633958 -38371.549099
29.98 -38253.351648 -38253.271135 -38253.195882  7186.511174 -38253.271135
30.00 -38305.305520 -38305.195131 -38305.090002  7234.784841 -38305.195131

Setting up a λ-schedule
-----------------------

* ``schedule`` - set the λ-schedule which specifies how λ morphs between
  the reference and perturbed molecules.

A λ-schedule (represented using the :class:`sire.cas.LambdaSchedule` class)
specifies how the λ-parameter morphs from the reference to the perturbed
molecules. The λ-schedule achieves this...

WRITE MORE ABOUT THE λ-schedule

Softening potentials
--------------------

* ``shift_delta`` - set the ``shift_delta`` parameter which is used to
  control the electrostatic and van der Waals softening potential that
  smooths the creation or deletion of ghost atoms. This is a floating
  point number that defaults to ``1.0``, which should be good for
  most perturbations.

* ``coulomb_power`` - set the ``coulomb_power`` parameter which is used
  to control the electrostatic softening potential that smooths the
  creation and deletion of ghost atoms. This is an integer that defaults
  to ``0``, which should be good for most perturbations.

Low level interface
-------------------

The high-level interface is just a set of convienient wrappers around the
OpenMM objects which are used to run the simulation. If you convert
any set of views (or view) that contains merged molecules, then an
alchemical OpenMM context will be returned.

>>> context = sr.convert.to(mols, "openmm")
>>> print(context)
OUTPUT

The context is held in a low-level class,
:class:`~sire.Convert.SireOpenMM.SOMMContext`, inherits from the
standard `OpenMM Context <https://docs.openmm.org/latest/api-python/generated/openmm.openmm.Context.html#openmm.openmm.Context>`__
class.

The class adds some additional metadata and control functions that are needed
to update the atomic parameters in the OpenMM Context to represent the
molecular system at different values of λ.

The key additional functions provided by :class:`~sire.Convert.SireOpenMM.SOMMContext`
are;

* :func:`~sire.Convert.SireOpenMM.SOMMContext.get_lambda` - return the
  current value of λ for the context.
* :func:`~sire.Convert.SireOpenMM.SOMMContext.set_lambda` - set the
  new value of λ for the context. Note that this should only really
  be used to change λ to evaluate energies at different λ-windows.
  It is better to re-create the context if you want to simulate
  at a different λ-value.
* :func:`~sire.Convert.SireOpenMM.SOMMContext.get_lambda_schedule` - return the
  λ-schedule used to control the morph.
* :func:`~sire.Convert.SireOpenMM.SOMMContext.set_lambda_schedule` - set the
  λ-schedule used to control the morph.
* :func:`~sire.Convert.SireOpenMM.SOMMContext.get_energy` - return the
  current potential energy of the context. This will be in :mod:`sire`
  units if ``to_sire_units`` is ``True`` (the default).
