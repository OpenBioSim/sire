==================
OpenMM Integration
==================

The :doc:`dynamics functionality in sire <../tutorial/part05/05_dynamics>`
is based on tight integration with `OpenMM <https://openmm.org>`__.

This is achieved via code in :mod:`sire.convert` that can convert
:mod:`sire` objects to their `OpenMM <https://openmm.org>`__
equivalents.

You can convert any molecule, molecule view, collection of molecules or
system objects into an OpenMM object using :func:`sire.convert.to`
and passing in ``openmm`` as the desired format.

>>> import sire as sr
>>> mols = sr.load(sr.expand(sr.tutorial_url, "kigaki.gro", "kigaki.top"),
...                show_warnings=False)
>>> omm = sr.convert.to(mols, "openmm")
>>> print(omm)
<openmm.openmm.Context; proxy of <Swig Object of type 'OpenMM::Context *' at 0x2846e5cb0> >

The result is an `OpenMM Context <https://docs.openmm.org/latest/api-python/generated/openmm.openmm.Context.html#openmm.openmm.Context>`__
which can be used just like any other ``Context`` object.

For example, you can extract the
`OpenMM System <https://docs.openmm.org/latest/api-python/generated/openmm.openmm.System.html#openmm.openmm.System>`__
via

>>> print(omm.getSystem())
<openmm.openmm.System; proxy of <Swig Object of type 'OpenMM::System *' at 0x2846e5d40> >

or the
`OpenMM Integrator <https://docs.openmm.org/latest/api-python/library.html#integrators>`__

>>> print(omm.getIntegrator())
<openmm.openmm.VerletIntegrator; proxy of <Swig Object of type 'OpenMM::VerletIntegrator *' at 0x2846e5a70> >

or the
`OpenMM Platform <https://docs.openmm.org/latest/api-python/generated/openmm.openmm.Platform.html#openmm.openmm.Platform>`__

>>> print(omm.getPlatform())
<openmm.openmm.Platform; proxy of <Swig Object of type 'OpenMM::Platform *' at 0x2846e5ce0> >

as you could for any other
`OpenMM Context <https://docs.openmm.org/latest/api-python/generated/openmm.openmm.Context.html#openmm.openmm.Context>`__.

You could then recombine these with your own choice of ``Integrator`` or
``Platform`` if you want to use something different to what was created
by :mod:`sire`.

Choosing options
----------------

:mod:`sire` made a series of choices for various OpenMM options and parameters
when it created the OpenMM Context. These choices were based on sensible
defaults derived from what it could find within the molecules/system being
converted.

You can override these choices by passing in a dictionary of key-value pairs
as the ``map`` option to :func:`sire.convert.to`.

For example, you can choose the integrator by setting a value
for the ``integrator`` key, and set the temperature and pressure
by setting values for the ``temperature`` and ``pressure`` keys.

>>> m = {"integrator": "langevin_middle",
...      "temperature": 25*sr.units.celsius,
...      "pressure": 1*sr.units.atm}
>>> omm = sr.convert.to(mols, "openmm", map=m)
>>> print(omm.getIntegrator())
<openmm.openmm.LangevinMiddleIntegrator; proxy of <Swig Object of type 'OpenMM::LangevinMiddleIntegrator *' at 0x295a07cc0> >

Available keys and allowable values are listed below.

+---------------------------+----------------------------------------------------------+
| Key                       | Valid values                                             |
+===========================+==========================================================+
| com_reset_frequency       | The frequency at which the ``CMMotionRemover`` acts to   |
|                           | remove center of mass relative motion. If this is not    |
|                           | set (the default) then center of mass motion is not      |
|                           | removed.                                                 |
+---------------------------+----------------------------------------------------------+
| constraint                | Type of constraint to use for bonds and/or angles.       |
|                           | Valid strings are ``none``, ``h-bonds``, ``bonds``,      |
|                           | ``h-bonds-h-angles`` and ``bonds-h-angles``.             |
+---------------------------+----------------------------------------------------------+
| coulomb_power             | The coulomb power parameter used by the softening        |
|                           | potential used to soften interactions involving          |
|                           | ghost atoms.                                             |
+---------------------------+----------------------------------------------------------+
| cutoff                    | Size of the non-bonded cutoff, e.g.                      |
|                           | ``7.5*sr.units.angstrom``                                |
+---------------------------+----------------------------------------------------------+
| cutoff_type               | Type of cutoff, e.g. ``PARTICLE_MESH_EWALD``, ``PME``,   |
|                           | ``NO_CUTOFF``, ``REACTION_FIELD``, ``RF``, ``EWALD``     |
+---------------------------+----------------------------------------------------------+
| cpu_pme                   | Boolean value, e.g. ``True`` or ``False`` as to whether  |
|                           | or not PME should be evaluated on the CPU rather than    |
|                           | on the GPU.                                              |
+---------------------------+----------------------------------------------------------+
| device                    | Any valid OpenMM device number or device string, e.g.    |
|                           | ``0`` would select device 0                              |
+---------------------------+----------------------------------------------------------+
| dielectric                | Dielectric value if a reaction field cutoff is used,     |
|                           | e.g. ``78.0``                                            |
+---------------------------+----------------------------------------------------------+
| fixed                     | The atoms in the system that should be fixed (not moved) |
+---------------------------+----------------------------------------------------------+
| ignore_perturbations      | Whether or not to ignore any perturbations and only set  |
|                           | up a perturbable molecule as a non-perurbable molecule   |
|                           | from only the reference state.                           |
+---------------------------+----------------------------------------------------------+
| integrator                | The MD integrator to use, e.g.                           |
|                           | ``verlet``, ``leapfrog``, ``langevin``,                  |
|                           | ``langevin_middle``, ``nose_hoover``,                    |
|                           | ``brownian``, ``andersen``                               |
+---------------------------+----------------------------------------------------------+
| friction                  | Friction value for the integrator, in inverse time, e.g. |
|                           | ``5.0 / sr.units.picosecond``                            |
+---------------------------+----------------------------------------------------------+
| lambda                    | The λ-value at which to set up the system (assuming this |
|                           | contains any perturbable molecules or restraints)        |
+---------------------------+----------------------------------------------------------+
| platform                  | Any valid OpenMM platform string, e.g. ``CUDA``,         |
|                           | ``OpenCL``, ``Metal``, ```CPU``, ``Reference``           |
+---------------------------+----------------------------------------------------------+
| precision                 | Any valid OpenMM platform precision value, e.g.          |
|                           | ``single``, ``mixed`` or ``double``.                     |
+---------------------------+----------------------------------------------------------+
| pressure                  | Any pressure value, e.g. ``1*sr.units.atm``              |
+---------------------------+----------------------------------------------------------+
| restraints                | The :class:`~sire.mm.Restraints` object (or list of      |
|                           | objects) that describe the restraints that should be     |
|                           | added to the system.                                     |
+---------------------------+----------------------------------------------------------+
| schedule                  | The :class:`~sire.cas.LambdaSchedule` to use that        |
|                           | controls how parameters are modified with λ              |
+---------------------------+----------------------------------------------------------+
| shift_delta               | The shift_delta parameter to use for the softening       |
|                           | potential used to soften interactions involving          |
|                           | ghost atoms.                                             |
+---------------------------+----------------------------------------------------------+
| space                     | Space in which the simulation should be conducted, e.g.  |
|                           | `sr.vol.Cartesian`                                       |
+---------------------------+----------------------------------------------------------+
| swap_end_states           | Whether to swap the end states of a perturbable molecule |
|                           | (i.e. treat the perturbed state as the reference state   |
|                           | and vice versa).                                         |
+---------------------------+----------------------------------------------------------+
| temperature               | Any temperature value, e.g. ``25*sr.units.celsius``      |
+---------------------------+----------------------------------------------------------+
| threads                   | The number of threads to use in the CPU platform         |
+---------------------------+----------------------------------------------------------+
| timestep                  | Time between integration steps, e.g.                     |
|                           | ``2 * sr.units.femtosecond``                             |
+---------------------------+----------------------------------------------------------+
| use_dispersion_correction | Whether or not to use the dispersion correction to       |
|                           | deal with cutoff issues. This is very expensive.         |
+---------------------------+----------------------------------------------------------+

Higher level API
----------------

The :class:`~sire.mol.Dynamics` object and :func:`~sire.mol.SelectorMol.dynamics`
function provides a higher level API for running molecular dynamics using the
`OpenMM Context <https://docs.openmm.org/latest/api-python/generated/openmm.openmm.Context.html#openmm.openmm.Context>`__
created by :mod:`sire`.

You create a :class:`~sire.mol.Dynamics` object by calling the
:func:`~sire.mol.SelectorMol.dynamics` function on the molecule,
molecule view, collection or system that you want to simulate.
For example

>>> d = mols.dynamics()

You can use this object to query the options that were passed into OpenMM.

>>> print(d.ensemble())
microcanonical (NVE) ensemble

You can set most of the OpenMM options via arguments to the :func:`~sire.mol.SelectorMol.dynamics`
function, e.g.

>>> d = mols.dynamics(temperature="25oC")
>>> print(d.ensemble())
canonical (NVT) ensemble { temperature = 298.15 C }

... note::

    The function will automatically convert strings to units if these are
    needed, e.g. ``25oC`` will automatically be converted to 25 Celsius.

or

>>> d = mols.dynamics(timestep="4fs", lambda_value=0.5)

You can also set OpenMM options by passing the dictionary of key-value pairs
as the ``map`` option.

>>> d = mols.dynamics(map={"temperature": "25oC"})
>>> print(d.ensemble())
canonical (NVT) ensemble { temperature = 298.15 C }

.. note::

   :mod:`sire` automatically chooses the right OpenMM Integrator and
   barostat options based on the ensemble parameters.

It is a mistake to use an OpenMM Integrator that is not suited
for the chosen ensemble.

>>> d = mols.dynamics(temperature="25oC", integrator="verlet")
ValueError: You cannot use a verlet integrator with the ensemble canonical (NVT) ensemble { temperature = 298.15 C }

You can also query other parameters.

>>> d = mols.dynamics(timestep="1fs")
>>> print(d.constraint())
none
>>> d = mols.dynamics(timestep="5fs")
>>> print(d.constraint())
bonds-h-angles
>>> print(d.timestep())
0.005 ps

Forcefield properties are automatically set based on the properties
contained by the molecules. You can get a summary of these properties
using the :func:`~sire.mol.Dynamics.info` function.

>>> print(d.info())
ForceFieldInfo(
  space=PeriodicBox( ( 48.3263, 48.3263, 48.3263 ) ),
  cutoff_type=PME,
  cutoff=7.5 Å,
  params=Properties( tolerance => 0.0001 ),
  detail=MM ForceField{ amber::ff,
               combining_rules = arithmetic,
               1-4 scaling = 0.833333, 0.5,
               nonbonded = coulomb, lj,
               bond = harmonic, angle = harmonic,
               dihedral = cosine }
)

Some of these properties, such as ``detail``, come from the forcefield
parameters of the converted molecules. Others, such as the
``cutoff_type`` and ``cutoff`` are passed from the options given
by the user (or derived as defaults). The ``space`` property is
extracted from the :class:`~sire.system.System` if that is passed,
or is found from the ``space`` property from the first molecule that
contains such a property. Sometimes, particularly if you aren't using
a :class:`~sire.system.System`, it can be a good idea to manually
set the ``space``, e.g. to :class:`~sire.vol.Cartesian` if you are
running a gas-phase simulation. In this case setting the
``cutoff_type`` to ``NO_CUTOFF`` will set the cutoff to a sufficiently
large value so that the effect is that there is no cutoff. Setting
the ``space`` to :class:`~sire.vol.Cartesian` will require disabling
``PME``, as this cutoff type requires a periodic space. Instead, choose
a cutoff type like reaction field.

>>> d = mols.dynamics(map={"space": sr.vol.Cartesian(),
...                        "cutoff_type": "NO_CUTOFF"})
>>> print(d.info())
ForceFieldInfo(
  space=Infinite cartesian space,
  cutoff_type=NO_CUTOFF,
  detail=MM ForceField{ amber::ff,
               combining_rules = arithmetic,
               1-4 scaling = 0.833333, 0.5,
               nonbonded = coulomb, lj,
               bond = harmonic, angle = harmonic,
               dihedral = cosine }
)
>>> d = mols.dynamics(map={"space": sr.vol.Cartesian(),
...                        "cutoff_type": "RF"})
>>> print(d.info())
ForceFieldInfo(
  space=Infinite cartesian space,
  cutoff_type=REACTION_FIELD,
  cutoff=7.5 Å,
  params=Properties( dielectric => 78.3 ),
  detail=MM ForceField{ amber::ff,
               combining_rules = arithmetic,
               1-4 scaling = 0.833333, 0.5,
               nonbonded = coulomb, lj,
               bond = harmonic, angle = harmonic,
               dihedral = cosine }
)

Running dynamics and saving frames and energies
-----------------------------------------------

You can run dynamics via the :func:`~sire.mol.Dynamics.run` function, e.g.

>>> d = mols.dynamics(timestep="4fs", temperature="25oC")
>>> d.run("100ps")

would run 100 picoseconds of dynamics.

At the end, you can extract the final system using the
:func:`~sire.mol.Dynamics.commit` function, e.g.

>>> mols = d.commit()

You can set the frequency at which trajectory frames and energies are saved
via the ``save_frequency`` argument, e.g.

>>> d.run("100ps", save_frequency="10ps")

would save energies and trajectory frames every 10 picoseconds. You can
specifiy different frequencies for these via the
``energy_frequency`` and/or ``frame_frequency`` arguments, e.g.

>>> d.run("1ns", energy_frequency="1ps", frame_frequency="100ps")

would save energies every picosecond and frames every 100 picoseconds.

By default, only coordinates are saved. You can choose to save velocities
as well by setting ``save_velocities=True``, e.g.

>>> d.run("10ps", save_frequency="1ps", save_velocities=True)

By default, energies are saved only for the simulated λ-value of the
system. You can request energies to be saved for other λ-values using
the ``lambda_windows`` argument, e.g.

>>> d.run("100ps", energy_frequency="1ps", lambda_windows=[0.0, 0.5, 1.0])

would save the energies at λ-values 0.0, 0.5 and 1.0 for every picosecond
of the trajectory. You can pass in as many or few λ-windows as you wish.

The coordinate/velocity frames are saved to the ``trajectory`` property of
the molecules, and are accessible identically to trajectories loaded
from files (e.g. via that property of the ``.trajectory()`` function).

The energies are saved to the ``energy_trajectory`` property of the
returned molecules, and accessible via that property or the
:func:`~sire.system.System.energy_trajectory` function.

