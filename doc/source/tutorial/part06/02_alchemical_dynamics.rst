===================
Alchemical Dynamics
===================

You can create an alchemical molecular dynamics simulation in exactly the
same way as you would a normal molecular dynamics simulation. There
are two options:

* Use the high-level interface based on the
:func:`~sire.system.System.minimisation` and
:func:`~sire.system.System.dynamics` functions.
* Use the low-level interface that works directly with native OpenMM objects.

High level interface
--------------------

The simplest route is to use the high-level interface. Calling
:func:`~sire.system.System.minimisation` or
:func:`~sire.system.System.dynamics` on any collection of molecules (or
view) that contains any perturbable molecule will automatically create
an alchemical simulation. There are extra options that you can pass to
the simulation that will control how it is run:

* ``lambda_value`` - this sets the global λ-value for the simulation.
  λ is a parameter that controls the morph from the reference molecule
  (at λ=0) to the perturbed molecule (at λ=1).

* ``swap_end_states`` - if set to ``True``, this will swap the end states
  of the perturbation. The morph will run from the perturbed molecule
  (at λ=0) to the reference molecule (at λ=1). Note that the coordinates
  of the perturbed molecule will be used in this case to start the
  simulation.

* ``schedule`` - set the λ-schedule which specifies how λ morphs between
  the reference and perturbed molecules.

* ``shift_delta`` - set the ``shift_delta`` parameter which is used to
  control the electrostatic and van der Waals softening potential that
  smooths the creation or deletion of ghost atoms. This is a floating
  point number that defaults to ``1.0``, which should be good for
  most perturbations.

* ``coulomb_power`` - set the ``coulomb_power`` parameter which is used
  to control the electrostatic softening potential that smooths the
  creation and deletion of ghost atoms. This is an integer that defaults
  to ``0``, which should be good for most perturbations.

For example, we could minimise our alchemical system at λ=0.5 using

>>> mols = mols.minimisation(lambda_value=0.5).run().commit()

We can then run some dynamics at this λ-value using

>>> d = mols.dynamics(lambda_value=0.5, timestep="4fs", temperature="25oC")
>>> d.run("10ps", save_frequency="1ps")
>>> mols = d.commit()

The result of dynamics is a trajectory run at this λ-value. You can view the
trajectory as you would any other, e.g.

>>> mols.view()

You could view specifically the perturbable molecules using

>>> mols["molecule property is_perturbable"].view()

Energy Trajectories
-------------------

In addition to a coordinates trajectory, dynamics also produces an
energy trajectory. This is the history of kinetic and potential energies
sampled by the molecules during the trajectory. You can access this
energy trajectory via the :func:`~sire.system.System.energy_trajectory`
function, e.g.

>>> df = mols.energy_trajectory()
>>> print(df)
PANDAS DATAFRAME

By default, this trajectory is returned as a
`pandas DataFrame <https://pandas.pydata.org/pandas-docs/stable/reference/frame.html>`__.
You can get the underlying :class:`sire.maths.EnergyTrajectory` function
by passing ``to_pandas=False``, e.g.

>>> t = mols.energy_trajectory(to_pandas=False)
>>> print(t)
EnergyTrajectory()

You calculate free energies by evaluating the potential energy for different
values of λ during dynamics. You can control which values of λ are used
(the so-called "λ-windows") by setting the ``lambda_windows`` argument, e.g.

>>> d = mols.dynamics(lambda_value=0.5, timestep="4fs", temperature="25oC")
>>> d.run("10ps", save_frequency="1ps", lambda_windows=[0.4, 0.6])
>>> mols = d.commit()
>>> print(mols.get_energy_trajectory())

This has run a new trajectory, evaluating the potential energy at the
simulation λ-value (0.5) as well as at λ-windows 0.4 and 0.6. You can pass in
as many or as few λ-windows as you want.

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
...       lambda_windows=[0.4, 0.6], save_velocities=False)
>>> mols = d.commit()
>>> print(mols.get_energy_trajectory())

.. note::

   The ``save_velocities`` option turns on or off the saving of atomic
   velocities to the frame trajectory. Typically you don't need to
   record velocities, so it is safe to switch them off. This can
   reduce memory consumption and also slightly speed up your simulation.

Setting up a λ-schedule
-----------------------

A λ-schedule (represented using the :class:`sire.cas.LambdaSchedule` class)
specifies how the λ-parameter morphs from the reference to the perturbed
molecules. The λ-schedule achieves this...

WRITE MORE ABOUT THE λ-schedule

Low level interface
-------------------

The high-level interface is just a set of convienient wrappers around the
OpenMM objects which are used to run the simulation. If you convert
any set of views (or view) that contains perturbable molecules, then an
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
