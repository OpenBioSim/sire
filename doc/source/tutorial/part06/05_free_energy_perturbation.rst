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


