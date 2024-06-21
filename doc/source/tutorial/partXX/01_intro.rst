============
Introduction
============

The ``sire`` QM/MM implementation takes advantage of the new means of writing
`platform independent force calculations <http://docs.openmm.org/development/developerguide/09_customcppforceimpl.html>`_
introduced in `OpenMM <http://openmm.org/>`_ 8.1. This allows us to interface
with any external package to modify atomic forces within the ``OpenMM`` context.
While OpenMM already directly supports ML/MM simulations via the `OpenMM-ML <https://github.com/openmm/openmm-ml>`_
package, it is currently limited to specific backends and only supports mechanical
embedding. The ``sire`` QM/MM implementation provides a simple way to interface
with any external QM package using a simple Python callback. This approach is
designed with generarilty and flexibility in mind, rather than performance,
allowing a user to quickly prototype new ideas.

Creating a QM engine
--------------------

In order to run QM/MM with ``sire``, we first need to create a QM engine. This
is passed as a keword argument to the ``dynamics`` function and is used to
perform the QM part of the calculation at each timestep.

As an example, we will consider the case of running a QM/MM simulation of alanine
dipeptide in water. First, let us load the molecular system:

>>> import sire as sr
>>> mols = sr.load_test_files("ala.crd", "ala.top")

We now need to set up the molecular system for the QM/MM simulation and create
an engine to perform the calculation:

>>> qm_mols, engine = sr.qm.create_engine(
>>> ...     mols,
>>> ...     mols[0],
>>> ...     py_object,
>>> ...     callback="callback",
>>> ...     cutoff="7.5A",
>>> ...     neighbour_list_update_frequency=20,
>>> ...     mechanical_embedding=False,
>>> ... )

Here the first argument is the molecules that we are simulating, the second
selection coresponding to the QM region (here this is the first molecule).
The third argument is the Python object that will be used to perform the QM
calculation. The fourth argument is the name of the callback function that will
be used. If ``None``, then it assumed that the ``py_object`` itself is a callable,
i.e.  it is the callback function. The callback function should have the following
signature:

.. code-block:: python

    def callback(
        numbers_qm: List[int],
        charges_mm: List[float],
        xyz_qm: List[List[float]],
        xyz_mm: List[List[float]],
    ) -> Tuple[float, List[List[float]], List[List[float]]]:

The function takes the atomic numbers of the QM atoms, the charges of the MM
atoms in mod electron charge, the coordinates of the QM atoms in Angstrom, and
the coordinates of the MM atoms in Angstrom. It should return the calculated
energy in kJ/mol, the forces on the QM atoms in kJ/mol/nm, and the forces
on the MM atoms in kJ/mol/nm. The remaining arguments are optional and specify
the QM cutoff distance, the neigbour list update frequency, and whether the
electrostatics should be treated with mechanical embedding. When mechanical
embedding is used, the electrostatics are treated at the MM level by ``OpenMM``.
Note that this doesn't change the signature of the callback function, i.e. it
will be passed empty lists for the MM specific arguments and should return an
empty list for the MM forces. Atomic positions passed to the callback function
will already be unwrapped with the QM region in the center.

The ``create_engine`` function returns a modified version of the molecules
containing a "merged" dipeptide that can be interpolated between MM and QM
levels of theory, along with the QM engine. This approach is extremely flexible
and allows the user to easily create a QM engine for a wide variety of QM packages.

Running a QM/MM simulation
--------------------------

In order to run a QM/MM simulation with ``sire`` we just need to specify our
QM engine when creating a dynamics object, for example:

>>> d = qm_mols.dynamics(
...     timestep="1fs",
...     constraint="none",
...     qm_engine=engine,
...     platform="cpu",
... )

For QM/MM simulations it is recommended to use a 1 femtosecond timestep and no
constraints. The simulation can then be run as usual:

>>> d.run("100ps", energy_frequency="1ps", frame_frequency="1ps")

This will run 100 picoseconds of dynamics, recording the energy and coordinates
every picosecond.

In next section we will show how to use `emle-engine <https://github.com/chemle/emle-engine>`_
package as QM engine via a simple specialisation of the interface shown above.
