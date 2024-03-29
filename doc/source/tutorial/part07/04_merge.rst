========================
Creating merge molecules
========================


Merged molecules are used in free energy calculations to represent the
perturbation between two molecules; the reference molecule (at λ=0)
and the perturbed molecule (at λ=1).

To start, let's load up two molecules, neopentane and methane, which we
will use to create a merged molecule to calculate the relative hydration
free energy.

>>> neopentane = sr.load_test_files("neopentane.prm7", "neopentane.rst")[0]
>>> neopentane.view()

.. image:: images/07_04_01.jpg
   :alt: A picture of neopentane

>>> methane = sr.load_test_files("methane.prm7", "methane.rst")[0]
>>> methane.view()

.. image:: images/07_04_02.jpg
   :alt: A picture of methane

Matching atoms
--------------

The first step to creating a merged molecule is to decide how atoms
should relate between the two end states. This is done by matching atoms
from the reference molecule (in this case neopentane) to the perturbed
molecule (in this case methane).

For example, let's say that we want the central carbon of neopentane to
perturb into the central carbon of methane. We could specify this by
creating a dictionary that says that the name of this atom in neopentane
should map to the name of the equivalent atom in methane.

We can get the name of these atoms using the 3D viewer, as shown above.
The central carbon of neopentane is called ``C1``, while the central
carbon of methane is also called ``C1``.

>>> matching = {"C1": "C1"}

Next, we would match each of the other carbon atoms in neopentane to
the hydrogen atoms in methane.

>>> matching["C2"] = "H2"
>>> matching["C3"] = "H3"
>>> matching["C4"] = "H4"
>>> matching["C5"] = "H5"
>>> print(matching)
{'C1': 'C1', 'C2': 'H2', 'C3': 'H3', 'C4': 'H4', 'C5': 'H5'}

Merging molecules
-----------------

We use the :func:`sire.morph.merge`



