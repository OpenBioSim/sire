=========
Sire-EMLE
=========

The ``sire`` QM/MM implementation takes advantage of the new means of writing
`platform independent force calculations <http://docs.openmm.org/development/developerguide/09_customcppforceimpl.html>`_
introduced in `OpenMM <http://openmm.org/>`_ 8.1. This allows us to interface
with any external package to modify atomic forces within the ``OpenMM`` context.
While OpenMM already directly supports ML/MM simulations via the `OpenMM-ML <https://github.com/openmm/openmm-ml>`_
package, it is currently limited to specific backends and only supports mechanical
embedding. The ``sire`` QM/MM implementation performs the QM calculation using
the `emle-engine <https://github.com/chemle/emle-engine>`_ package, which has
support for a wide range of backends and embedding models, importantly providing
a simple and efficient ML model for electrostatic embedding.

In order to use QM/MM functionality within ``sire`` you will first need to
create the following ``conda`` environment:

.. code-block:: bash

   $ git clone https://github.com/chemle/emle-engine.git
   $ cd emle-engine
   $ conda env create -f environment_sire.yaml
   $ conda activate emle-sire
   $ pip install -e .

In this tutorial, we will perform a short ML/MM simulation of alanine dipeptide
in water. First, let us load the molecular system:

>>> import sire as sr
>>> mols = sr.load_test_files("ala.crd", "ala.top")

Creating an EMLE calculator
---------------------------

Next we will create an ``emle-engine`` calculator to perform the QM (or ML) calculation
for the dipeptide along with the ML electrostatic embedding. Since this is a small molecule
it isn't beneficial to perform the calculation on a GPU, so we will use the CPU instead.

>>> from emle.calculator import EMLECalculator
>>> calculator = EMLECalculator(device="cpu")

By default, ``emle-engine`` will use `TorchANI <https://aiqm.github.io/torchani/>`_
as the backend for in vacuo calculation of energies and gradients. However,
it is possible to use a wide variety of other backends, including your own
as long as  it supports the standand `Atomic Simulation Environment (ASE) <https://wiki.fysik.dtu.dk/ase/>`_
`calculator <https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html>`_ interface.
For details, see the `backends <https://github.com/chemle/emle-engine#backends>`_
section of the ``emle-engine`` documentation. At present, the default embedding
model provided with ``emle-engine`` supports only the elements H, C, N, O, and S.
We plan on adding support for other elements in the near future.

Creating a QM engine
--------------------

We now need to set up the molecular system for the QM/MM simulation and create
an engine to perform the calculation:

>>> qm_mols, engine = sr.qm.emle(mols, mols[0], calculator, "7.5A", 20)

Here the first argument is the molecules that we are simulating, the second
selection coresponding to the QM region (here this is the first molecule), and
the third is calculator that was created above. The fourth and fifth arguments
are optional, and specify the QM cutoff distance and the neigbour list update
frequency respectively. (Shown are the default values.) The function returns a
modified version of the molecules containing a "merged" dipeptide that can be
interpolated between MM and QM levels of theory, along with an engine. The
engine registers a Python callback that uses ``emle-engine`` to perform the QM
calculation.

The selection syntax for QM atoms is extremely flexible. Any valid search string,
atom index, list of atom indicies, or molecule view/container that can be used.
Support for modelling partial molecules at the QM level is provided via the link
atom approach, via the charge shifting method. For details of this implementation,
see, e.g., the NAMD user guide `here <https://www.ks.uiuc.edu/Research/qmmm/>`_.
While we support multiple QM fragments, we do not currently support multiple
*independent* QM regions. We plan on adding support for this in the near future.

Running a QM/MM simulation
--------------------------

Next we need to create a dynamics object to perform the simulation. For QM/MM
simulations it is recommended to use a 1 femtosecond timestep and no constraints.
In this example we will use the ``lambda_interpolate`` keyword to  interpolate
the dipeptide potential between pure MM (λ=0) and QM (λ=1) over the course of
the simulation, which can be used for end-state correction of binding free
energy calculations.

>>> d = qm_mols.dynamics(timestep="1fs", constraint="none", qm_engine=engine, lambda_interpolate=[0, 1])

We can now run the simulation. The options below specify the run time, the
frequency at which trajectory frames are saved, and the frequency at which
energies are recorded. The ``energy_frequency`` also specifies the frequency
at which the λ value is updated.

>>> d.run("1ps", frame_frequency="0.05ps", energy_frequency="0.05ps")

.. note::

    Updating λ requires the updating of force field parameters in the ``OpenMM``
    context. For large systems, this can be quite slow so it isn't recommended
    to set the ``energy_frequency`` to a value that is too small. We have a custom
    `fork <https://github.com/chryswoods/openmm>`_ of ``OpenMM`` that provides a
    significant speedup for this operation by only updating a subset of the parameters.
    Installation instructions can be provided on request.

.. note::

    If you don't require a trajectory file, then better performance can be achieved
    leaving the ``frame_frequency`` keyword argument unset.

.. note::

    ``emle-engine`` currently requires the use of `librascal <https://lab-cosmo.github.io/librascal/#/>`_
    for the calculation of SOAP (Smooth Overlap of Atomic Positions) descriptors.
    This is a serial code, so you may see better performance by restricting the
    number of ``OpenMP`` threads to 1, e.g. by setting the ``OMP_NUM_THREADS``
    environment variable.

Once the simulation has finished we can get back the trajectory of energy values.
This can be obtained as a `pandas <https://pandas.pydata.org/>`_ ``DataFrame``,
allowing for easy plotting and analysis. The table below shows the instantaneous
kintetic and potential energies as a function of λ, along with the pure MM and
QM potential energies. (Times are in picoseconds and energies are in kcal/mol.)

>>> nrg_traj = d.energy_trajectory(to_pandas=True)
>>> print(nrg_traj)
           lambda           KE     PE(lambda)  PE(lambda=0)   PE(lambda=1)
time
6000.05  0.000000   980.181564   -6954.938694  -6954.938694 -318014.135823
6000.10  0.052632   871.904630  -23214.139963  -6843.385099 -317910.734657
6000.15  0.105263  1074.693130  -39796.029943  -7056.370765 -318111.343285
6000.20  0.157895   979.813677  -56061.595767  -6952.183998 -318008.475588
6000.25  0.210526  1009.571276  -72462.277097  -6981.451657 -318040.986409
6000.30  0.263158  1016.026458  -88842.745858  -6991.337337 -318046.238677
6000.35  0.315789  1003.273813 -105199.347795  -6976.690749 -318031.016925
6000.40  0.368421  1021.295211 -121583.564572  -6991.838146 -318041.438719
6000.45  0.421053  1027.366329 -137961.602333  -7000.530076 -318049.949920
6000.50  0.473684  1049.387973 -154355.318394  -7023.254018 -318072.387286
6000.55  0.526316  1040.626785 -170718.777695  -7016.367279 -318066.329145
6000.60  0.578947  1047.005579 -187097.460730  -7015.987089 -318076.072803
6000.65  0.631579  1030.218148 -203453.572350  -6997.132190 -318063.875864
6000.70  0.684211  1022.362023 -219819.959312  -6994.205184 -318058.533453
6000.75  0.736842  1044.950320 -236216.451165  -7012.311296 -318084.096807
6000.80  0.789474  1024.087813 -252561.720268  -6985.090189 -318055.746705
6000.85  0.842105  1056.241205 -268962.249393  -7016.702075 -318082.555659
6000.90  0.894737  1053.591066 -285328.646842  -7013.509852 -318075.626766
6000.95  0.947368  1033.013716 -301672.026582  -6986.164439 -318045.397622
6001.00  1.000000  1045.687318 -318056.550581  -6991.865785 -318056.550599

.. note::

   In the table above, the time doesn't start from zero because the example
   molecular system was loaded from an existing trajectory restart file.

Interfacing with OpenMM-ML
--------------------------

In the example above we used a sire dynamics object ``d`` to run the simulation.
This is wrapper around a standard OpenMM context object, providing a simple
convenience functions to make it easier to run and analyse simulations. However,
if you are already familiar with OpenMM, then it is possible to use ``emle-engine``
with OpenMM directly. This allows for fully customised simulations, or the use
of `OpenMM-ML <https://github.com/openmm/openmm-ml>`_ as the backend for
calculation of the intramolecular force for the QM region.

To use ``OpenMM-ML`` as the backend for the QM calculation, you will first need
to install the package:

.. code-block:: bash

   $ conda install -c conda-forge openmm-ml

Next, you will need to create an ``MLPotential`` for desired backend. Here we
will use the ``ani2x``, as was used for the ``EMLECalculator`` above. The

>>> import openmm
>>> from openmmml import MLPotential
>>> potential = MLPotential("ani2x")

Since we are now using the ``MLPotential`` for the QM calculation, we need to
create a new ``EMLECalculator`` object with no backend, i.e. one that only
computes the electrostatic embedding:

>>> calculator = EMLECalculator(backend=None, device="cpu")

Next we create a new engine bound to the calculator:

>>> qm_mols, engine = sr.qm.emle(mols, mols[0], calculator, "7.5A", 20)

Rather than using this engine with a ``sire`` dynamics object, we can instead
extract the underlying ``OpenMM`` force object and add it to an existing
``OpenMM`` system. The forces can be extracted from the engine as follows:

>>> emle_force, interpolation_force = engine.get_forces()

The ``emle_force`` object is the ``OpenMM`` force object that calculates the
electrostatic embedding interaction. The ``interpolation_force`` is a null
``CustomBondForce`` object that contains a ``lambda_emle`` global parameter
than can be used to scale the electrostatic embedding interaction. (By default,
this is set to 1, but can be set to any value between 0 and 1.)

.. note::

    The ``interpolation_force`` has no energy contribution. It is only required
    as there is currently no way to add global parameters to the ``EMLEForce``.

Since we want to use electrostatic embedding, we will also need to zero the charges
on the atoms within the QM region before creating an ``OpenMM`` system. If not,
then we would also calculate the mechanical embedding interaction. This can be
done using the ``qm_mols`` object generated above. This system is *perturbable*
so can be converted between an MM reference state and QM perturbed state. Here
we require the perturbed state, which has zeroed charges for the QM region:

>>> qm_mol = sr.morph.link_to_perturbed(qm_mols[0])
>>> qm_mols.update(qm_mol)

We now write the modified system to an AMBER format topology and coordinate file
so that we can load them with ``OpenMM``:

>>> sr.save(qm_mols, "ala_qm.prm7")
>>> sr.save(qm_mols, "ala_qm.rst7")

We can now read them back in with ``OpenMM``:

>>> prmtop = openmm.app.AmberPrmtopFile("ala_qm.prm7")
>>> inpcrd = openmm.app.AmberInpcrdFile("ala_qm.rst7")

Next we use the ``prmtop`` to create the MM system:

>>> mm_system = prmtop.createSystem(
...     nonbondedMethod=openmm.app.PME,
...     nonbondedCutoff=1 * openmm.unit.nanometer,
...     constraints=openmm.app.HBonds,
... )

In oder to create the ML system, we first define the ML region. This is a list
of atom indices that are to be treated with the ML model.

>>> ml_atoms = list(range(qm_mols[0].num_atoms()))

We can now create the ML system:

>>> ml_system = potential.createMixedSystem(
...     topology, mm_system, ml_atoms, interpolate=True
... )

By setting ``interpolate=True`` we are telling the ``MLPotential`` to create
a *mixed* system that can be interpolated between MM and ML levels of theory
using the ``lambda_interpolate`` global parameter. (By default this is set to 1.)

.. note::

    If you choose not to add the ``emle`` interpolation force to the system, then
    the ``EMLEForce`` will also use the ``lambda_interpolate`` global parameter.
    This allows for the electrostatic embedding to be alongside or independent of
    the ML model.

We can now add the ``emle`` forces to the system:

>>> ml_system.addForce(emle_force)
>>> ml_system.addForce(interpolation_force)

In oder to run a simulation we need to create an integrator and context. First
we create the integrator:

>>> integrator = openmm.LangevinMiddleIntegrator(
...     300 * openmm.unit.kelvin,
...     1.0 / openmm.unit.picosecond,
...     0.002 * openmm.unit.picosecond,
... )

And finally the context:

>>> context = openmm.Context(ml_system, integrator)
>>> context.setPositions(inpcrd.positions)
