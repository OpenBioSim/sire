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
>>> energy_files = glob("energy*.s3")
>>> energy_files.sort()
>>> for energy_file in energy_files:
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

.. note::

   Be careful to load the DataFrames in λ-order. The ``glob`` function
   can return the files in a random order, hence why we need to sort
   this list. This sort only works because we have used a naming convention
   for the files that puts them in λ-order. They must be in the right
   order or else alchemlyb will calculate the free energy incorrectly
   (it uses the column-order rather than the λ-order when calculating
   free energies).

Next, we will join together all of these DataFrames into a single
DataFrame.

>>> import pandas as pd
>>> df = pd.concat(dfs)
>>> print(df)
                         0.00          0.05  0.10  0.15  0.20  0.25  0.30  0.35  ...  0.65  0.70  0.75  0.80  0.85  0.90          0.95          1.00
time fep-lambda                                                                  ...
2.1  0.0        -40299.914145 -40300.859816   NaN   NaN   NaN   NaN   NaN   NaN  ...   NaN   NaN   NaN   NaN   NaN   NaN           NaN           NaN
2.2  0.0        -40028.164623 -40028.811537   NaN   NaN   NaN   NaN   NaN   NaN  ...   NaN   NaN   NaN   NaN   NaN   NaN           NaN           NaN
2.3  0.0        -39862.145264 -39862.881804   NaN   NaN   NaN   NaN   NaN   NaN  ...   NaN   NaN   NaN   NaN   NaN   NaN           NaN           NaN
2.4  0.0        -39795.910799 -39796.438210   NaN   NaN   NaN   NaN   NaN   NaN  ...   NaN   NaN   NaN   NaN   NaN   NaN           NaN           NaN
2.5  0.0        -39597.177529 -39597.256804   NaN   NaN   NaN   NaN   NaN   NaN  ...   NaN   NaN   NaN   NaN   NaN   NaN           NaN           NaN
...                       ...           ...   ...   ...   ...   ...   ...   ...  ...   ...   ...   ...   ...   ...   ...           ...           ...
26.6 1.0                  NaN           NaN   NaN   NaN   NaN   NaN   NaN   NaN  ...   NaN   NaN   NaN   NaN   NaN   NaN -37356.997418 -37357.661840
26.7 1.0                  NaN           NaN   NaN   NaN   NaN   NaN   NaN   NaN  ...   NaN   NaN   NaN   NaN   NaN   NaN -37368.648948 -37368.924986
26.8 1.0                  NaN           NaN   NaN   NaN   NaN   NaN   NaN   NaN  ...   NaN   NaN   NaN   NaN   NaN   NaN -37380.838240 -37381.920923
26.9 1.0                  NaN           NaN   NaN   NaN   NaN   NaN   NaN   NaN  ...   NaN   NaN   NaN   NaN   NaN   NaN -37409.937188 -37410.780865
27.0 1.0                  NaN           NaN   NaN   NaN   NaN   NaN   NaN   NaN  ...   NaN   NaN   NaN   NaN   NaN   NaN -37412.237619 -37412.812414

.. note::

   Do not worry about the large number of ``NaN`` values. These just show that
   we have only calculated free energy differences along the diagonal of this
   DataFrame, i.e. only between the simulated and neighbouring λ-windows.

Now we can tell alchemlyb to calculate the free energy using the BAR method.

>>> from alchemlyb.estimators import BAR
>>> b = BAR()
>>> b.fit(df)
>>> print(b.delta_f_.loc[0.00, 1.00])
-2.9285302096352157

You can get a convergence plot, showing how the free energy changes as
a function of the simulation length using the ``convergence_plot`` function.

>>> from alchemlyb.convergence import forward_backward_convergence
>>> from alchemlyb.visualisation import plot_convergence
>>> f = forward_backward_convergence(dfs, "bar")
>>> plot_convergence(f)

.. image:: images/06_05_01.jpg
   :alt: Convergence of the free energy estimate as a function of the fraction of simulation length

All of this shows that the relative free energy for the perturbation of
ethane to methanol in water is about -2.9 kcal mol-1.

To get the relative hydration free energy, we would need to complete the
cycle by calculating the relative free energy for the perturbation in the
gas phase. We could do this using this code (which is almost identical to
above, except we only simulate the perturbable molecule, and save
the :class:`~sire.maths.EnergyTable` objects to ``energy_gas_{lambda}.s3``
instead of ``energy_{lambda}.s3``).

>>> import sire as sr
>>> mols = sr.load(sr.expand(sr.tutorial_url, "merged_molecule.s3"))
>>> mol = mols.molecule("molecule property is_perturbable")
>>> for l in range(0, 105, 5):
...     # turn l into the lambda value by dividing by 100
...     lambda_value = l / 100.0
...     print(f"Simulating lambda={lambda_value:.2f}")
...     # create a dynamics object for the system
...     d = mol.dynamics(timestep="1fs", temperature="25oC",
...                      lambda_value=lambda_value)
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
...                    f"energy_gas_{lambda_value:.2f}.s3")

This should run more quickly than the simulation in water, e.g. about
15 seconds per window (at about 150 nanoseconds per day of sampling).

We can then analyse the results using the same analysis code, except we
switch to analysing the ``energy_gas_{lambda}.s3`` files instead.

>>> import sire as sr
>>> from glob import glob
>>> dfs = []
>>> energy_files = glob("energy_gas_*.s3")
>>> energy_files.sort()
>>> for energy_file in energy_files:
...     dfs.append(sr.stream.load(energy_file).to_pandas(to_alchemlyb=True, temperature="25oC"))
>>> import pandas as pd
>>> df = pd.concat(dfs)
>>> print(df)
                     0.00      0.05  0.10  0.15  0.20  0.25  0.30  0.35  0.40  ...  0.60  0.65  0.70  0.75  0.80  0.85  0.90       0.95       1.00
time fep-lambda                                                                ...
2.1  0.0         7.707173  7.686233   NaN   NaN   NaN   NaN   NaN   NaN   NaN  ...   NaN   NaN   NaN   NaN   NaN   NaN   NaN        NaN        NaN
2.2  0.0         5.048624  5.303638   NaN   NaN   NaN   NaN   NaN   NaN   NaN  ...   NaN   NaN   NaN   NaN   NaN   NaN   NaN        NaN        NaN
2.3  0.0         5.797784  5.446916   NaN   NaN   NaN   NaN   NaN   NaN   NaN  ...   NaN   NaN   NaN   NaN   NaN   NaN   NaN        NaN        NaN
2.4  0.0         4.080923  4.189452   NaN   NaN   NaN   NaN   NaN   NaN   NaN  ...   NaN   NaN   NaN   NaN   NaN   NaN   NaN        NaN        NaN
2.5  0.0         5.386196  5.271694   NaN   NaN   NaN   NaN   NaN   NaN   NaN  ...   NaN   NaN   NaN   NaN   NaN   NaN   NaN        NaN        NaN
...                   ...       ...   ...   ...   ...   ...   ...   ...   ...  ...   ...   ...   ...   ...   ...   ...   ...        ...        ...
26.6 1.0              NaN       NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN  ...   NaN   NaN   NaN   NaN   NaN   NaN   NaN  10.989612  11.139099
26.7 1.0              NaN       NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN  ...   NaN   NaN   NaN   NaN   NaN   NaN   NaN   9.000112   9.071128
26.8 1.0              NaN       NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN  ...   NaN   NaN   NaN   NaN   NaN   NaN   NaN  11.481547  11.700362
26.9 1.0              NaN       NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN  ...   NaN   NaN   NaN   NaN   NaN   NaN   NaN  11.801311  11.934601
27.0 1.0              NaN       NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN  ...   NaN   NaN   NaN   NaN   NaN   NaN   NaN  10.316325  10.489984
>>> from alchemlyb.estimators import BAR
>>> b = BAR()
>>> b.fit(df)
>>> print(b.delta_f_.loc[0.00, 1.00])
3.177901359408199

This shows that the relative free energy for the perturbation of ethane
to methanol in the gas phase is about 3.2 kcal mol-1. Subtracting this
from the free energy in water gives a relative hydration free energy of
about -6.1 kcal mol-1, which is in reasonable agreement with
`published results from other codes <https://www.pure.ed.ac.uk/ws/portalfiles/portal/75900057/20181010_Michel_reprod.pdf>`__
which are typically -6.1 kcal mol-1 to -6.2 kcal mol-1.
