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
...     mols,
...     mols[0],
...     py_object,
...     callback="callback",
...     cutoff="7.5A",
...     neighbour_list_frequency=20,
...     mechanical_embedding=False,
... )

Here the first argument is the molecules that we are simulating, the second
selection coresponding to the QM region (here this is the first molecule).
The selection syntax for QM atoms is extremely flexible. Any valid search string,
atom index, list of atom indicies, or molecule view/container that can be used.
Support for modelling partial molecules at the QM level is provided via the link
atom approach, via the charge shifting method. For details of this implementation,
see, e.g., the NAMD user guide `here <https://www.ks.uiuc.edu/Research/qmmm/>`_.
While we support multiple QM fragments, we do not currently support multiple
*independent* QM regions. We plan on adding support for this in the near future.
The third argument is the Python object that will be used to perform the QM
calculation. The fourth argument is the name of the callback function that will
be used. If ``None``, then it assumed that the ``py_object`` itself is a callable,
i.e.  it is the callback function. The callback function should have the following
signature:

.. code-block:: python

    from typing import List, Optional, Tuple

    def callback(
        numbers_qm: List[int],
        charges_mm: List[float],
        xyz_qm: List[List[float]],
        xyz_mm: List[List[float]],
        idx_mm: Optional[List[int]] = None,
    ) -> Tuple[float, List[List[float]], List[List[float]]]:

The function takes the atomic numbers of the QM atoms, the charges of the MM
atoms in mod electron charge, the coordinates of the QM atoms in Angstrom, and
the coordinates of the MM atoms in Angstrom. Optionally, it should also take the
indices of the true MM atoms (not link atoms or virtual charges) within the
QM/MM region. This is useful for obtaining any additional atomic properties
that may be required by the callback. The function should return the calculated
energy in kJ/mol, the forces on the QM atoms in kJ/mol/nm, and the forces
on the MM atoms in kJ/mol/nm. The remaining arguments are optional and specify
the QM cutoff distance, the neighbour list update frequency, and whether the
electrostatics should be treated with mechanical embedding. When mechanical
embedding is used, the electrostatics are treated at the MM level by ``OpenMM``.
Note that this doesn't change the signature of the callback function, i.e. it
will be passed empty lists for the MM specific arguments and should return an
empty list for the MM forces. Atomic positions passed to the callback function
will already be unwrapped with the QM region in the center. By default, no
neighbour list will be used. (The same thing can be achieved by passing
``neighbour_list_frequency=0``.) This is useful when using the engine as
a calculator for different input structures, where there may be no correlation
between coordinates. For regular molecular dynamics simulations, setting a
non-zero neighbour list frequency can improve performance.

The ``create_engine`` function returns a modified version of the molecules
containing a "merged" dipeptide that can be interpolated between MM and QM
levels of theory, along with the QM engine. This approach is extremely flexible
and allows the user to easily create a QM engine for a wide variety of QM packages.

Note that while the callback interface described above is designed to be used
for QM/MM, it is completely general so could be used to apply *any* external
force based on the local environment around a subset of atoms. For example, you
could apply a biasing potential on top of the regular MM force field.

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

If you are using the callback interface and wish to apply a force on top of the
existing MM force field, rather than perform QM/MM, then you can pass
``swap_end_states=True`` to the ``dynamics`` function. This will swap the QM and
MM end states of all *perturbable* molecules within ``qm_mols``, so that the MM
state corresponds to λ = 1. More details on on λ interpolation can be found in
the `next section <https://github.com/chemle/emle-engine>`_.

In next section we will show how to use `emle-engine <https://github.com/chemle/emle-engine>`_
package as QM engine via a simple specialisation of the interface shown above.
