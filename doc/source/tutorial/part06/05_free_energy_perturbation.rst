=============================================
Running a Free Energy Perturbation Simulation
=============================================

Now that you have your (possibly restrained) alchemical system, you can
run a free energy perturbation simulation. There are many ways to do this,
and the exact protocol will depend on what type of free energy you are
calculating, what level of accuracy you require and what compute resource
you have available. :mod:`sire` provides the building blocks to run a
simulation in any way you like, but in this tutorial we will use the
a very simple protocol. Take a look at
`BioSimSpace <https://biosimspace.openbiosim.org>`__ or the upcoming
`somd2 software <https://github.com/openbiosim/somd2>`__ if you are
interested in higher-level interfaces to this functionality that
automatically run more complex protocols.

At its simplest, a free energy perturbation simulation consists of a series
of alchemical molecular dynamics simulations run at different values
of the alchemical parameter, λ. Let's start with the ethane to methanol
merged molecule from the previous section.

>>> import sire as sr
>>> mols = sr.load(sr.expand(sr.tutorial_url, "merged_molecule.s3"))
>>> print(mols)
System( name=BioSimSpace_System num_molecules=4054 num_residues=4054 num_atoms=12167 )

And lets link the properties to the reference state.

>>> for mol in mols.molecules("molecule property is_perturbable"):
...     mols.update(mol.perturbation().link_to_reference().commit())

Next we will run through 20 evenly-space λ values from 0 to 1. We've picked
20 because it is a reasonable number that should over-sample the λ-coordinate.

For each λ-value, we will minimise the system, equilibrate it for 2 ps, then
run dynamics for 25 ps. We will calculate the energy of each simulation at
each λ-value, plus at the neighbouring λ-values. This will allow us to
calculate energy differences which we will later use to calculate the
free energy difference.

We will calculate these energies every 0.1 ps, and will write them at the
end of each block of dynamics via the :class:`~sire.maths.EnergyTrajectory`
output which we will stream as a :mod:`sire` object.

This will let us subsequently calculate the free energy across λ using
`alchemlyb <https://alchemlyb.readthedocs.io/en/latest/>`__.

>>> for l in range(0, 105, 5):
...     # turn l into the lambda value by dividing by 100
...     lambda_value = l / 100.0
...     print(f"Simulating lambda={lambda_value:.2f}")
...     # create a dynamics object for the system
...     d = mols.dynamics(timestep="1fs", temperature="25oC",
...                       lambda_value=lambda_value)
...     # minimise
...     d.minimise()
...     # equilibrate, not saving anything
...     d.run("2ps", save_frequency=0)
...     print("Equilibration complete")
...     print(d)
...     # get the values of lambda for neighbouring windows
...     lambda_windows = [lambda_value]
...     if lambda_value > 0:
...         lambda_windows.insert(0, (l-5)/100.0)
...     if lambda_value < 1:
...         lambda_windows.append((l+5)/100.0)
...     # run the dynamics, saving the energy every 0.1 ps
...     d.run("25ps", energy_frequency="0.1ps", frame_frequency=0,
...           lambda_windows=lambda_windows)
...     print("Dynamics complete")
...     print(d)
...     # stream the EnergyTable to a sire save stream object
...     sr.stream.save(d.commit().energy_trajectory(to_pandas=False),
...                    f"energy_{lambda_value:.2f}.s3")

This is a very simple protocol for a free energy simulation, with simulation
times chosen so that it will run relatively quickly. For example, with a
modern GPU this should run at about 60 nanoseconds of sampling per day,
so take only 45 seconds or so for each of the 21 λ-windows.

.. note::

   This is not enough sampling to calculate a properly converged free energy
   for most molecular systems. You would need to experiment with different
   equilibration and simulation times, or use a more sophisticated algorithm
   implemented via `BioSimSpace <https://biosimspace.openbiosim.org>`__ or
   `somd2 <https://github.com/openbiosim/somd2>`__.

.. note::

   We set the forwards and backwards λ-values to ``(l-5)/100.0`` and
   ``(l+5)/100.0`` so that they exactly match the way we set the reference
   λ-value. This is because alchemlyb matches energies using a
   floating point representation of λ, so small numerical differences
   between, e.g. ``(l+5)/100.0`` and ``lambda_value + 0.05`` will lead
   to errors (as in this case, ``(l+5)/100.0`` is exactly every 0.05, while
   ``lambda_value + 0.05`` has numerical error, meaning some values are
   at, e.g. ``0.15000000000000002``).

At the end of the simulation, you should have 21 energy files, one for each
λ-window. These are called ``energy_0.00.pq``, ``energy_0.05.pq``, ...,
``energy_1.00.pq``. They are in
`parquet format <https://parquet.apache.org/>`__ and are ready for processing
directly in the `alchemlyb <https://alchemlyb.readthedocs.io/en/latest/>`__ package.

Processing Free Energy Data using alchemlyb
--------------------------------------------

`alchemlyb <https://alchemlyb.readthedocs.io/en/latest/>`__ is a simple
alchemistry python package for easily analysing the results of free energy
simulations. They have excellent documentation on their website that you
can use if you want to go into the details of how to calculate free
energies.

Here, we will show a simple BAR analysis of the data that we have just
generated. We can do this because we have calculated data which
alchemlyb can convert into reduced potentials for each λ-window.

First, we need to import alchemlyb

>>> import alchemlyb

.. note::

   If you see an error then you may need to install (or reinstall)
   alchemlyb. You can do this using conda or mamba, e.g.
   ``mamba install -c conda-forge alchemlyb``.

Next, we will load all of the DataFrames up into alchemlyb dataframes.

>>> import sire as sr
>>> from glob import glob
>>> dfs = []
>>> for energy_file in glob("energy*.s3"):
...     dfs.append(sr.stream.load(energy_file).to_pandas(to_alchemlyb=True, temperature="25oC"))

.. note::

   We have to manually set the temperature to 25°C here because
   the EnergyTable doesn't (yet) contain the system temperature.
   This is something that we are looking to fix in a later release.

.. note::

   If you wanted, you could have put the dataframes generated above
   directly into the ``dfs`` list here, and not saved them to disk.
   However, this would risk you having to re-run
   all of the simulation if you wanted to change the analysis below.

Next, we will join together all of these DataFrames into a single
DataFrame.

>>> import pandas as pd
>>> df = pd.concat(dfs)
>>> print(df)
                         0.75          0.80          0.85  0.90          0.95          1.00  0.05  ...  0.25  0.40  0.45  0.50  0.55  0.30  0.35
time fep-lambda                                                                                    ...
2.1  0.8        -40294.138361 -40294.528643 -40294.867959   NaN           NaN           NaN   NaN  ...   NaN   NaN   NaN   NaN   NaN   NaN   NaN
2.2  0.8        -40132.809489 -40133.199771 -40133.479335   NaN           NaN           NaN   NaN  ...   NaN   NaN   NaN   NaN   NaN   NaN   NaN
2.3  0.8        -39940.111114 -39940.292266 -39940.422452   NaN           NaN           NaN   NaN  ...   NaN   NaN   NaN   NaN   NaN   NaN   NaN
2.4  0.8        -39709.709585 -39710.368748 -39711.006821   NaN           NaN           NaN   NaN  ...   NaN   NaN   NaN   NaN   NaN   NaN   NaN
2.5  0.8        -39786.131669 -39786.342696 -39786.532634   NaN           NaN           NaN   NaN  ...   NaN   NaN   NaN   NaN   NaN   NaN   NaN
...                       ...           ...           ...   ...           ...           ...   ...  ...   ...   ...   ...   ...   ...   ...   ...
4.6  1.0                  NaN           NaN           NaN   NaN -37920.811949 -37921.595875   NaN  ...   NaN   NaN   NaN   NaN   NaN   NaN   NaN
4.7  1.0                  NaN           NaN           NaN   NaN -37794.676672 -37794.833207   NaN  ...   NaN   NaN   NaN   NaN   NaN   NaN   NaN
4.8  1.0                  NaN           NaN           NaN   NaN -37823.895123 -37824.201037   NaN  ...   NaN   NaN   NaN   NaN   NaN   NaN   NaN
4.9  1.0                  NaN           NaN           NaN   NaN -37739.077963 -37739.503380   NaN  ...   NaN   NaN   NaN   NaN   NaN   NaN   NaN
5.0  1.0                  NaN           NaN           NaN   NaN -37724.349234 -37724.894154   NaN  ...   NaN   NaN   NaN   NaN   NaN   NaN   NaN

Now we can tell alchemlyb to calculate the free energy using the BAR method.

>>> from alchemlyb.estimators import BAR
>>> b = BAR()
>>> b.fit(df)
>>> print(b.delta_f_.loc[0.00, 1.00])
0.3736
