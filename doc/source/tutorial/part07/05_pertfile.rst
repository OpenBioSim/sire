==============================
Perturbation Files (pertfiles)
==============================

Perturbation files (or pertfiles) are an older mechanism that was used
in ``somd`` to create a merged molecule from a passed input molecule.
They are a simple text file that describes the perturbation in terms
of changing forcefield parameters. You can create a merged molecule
from a single molecule plus pertfile using the
:func:`sire.morph.create_from_pertfile` function.

>>> merged_mol = sr.morph.create_from_pertfile(mol, "neopentane_methane.pert")
>>> print(merged_mol.property("charge0"))

>>> print(merged_mol.property("charge1"))

.. note::

   This is an older mechanism that has many limitations due to the
   inherent limits of the pertfile format. It is provided to aid
   compatibility with older ``somd`` workflows, but is not
   recommended for new use cases.

Updating internals involving ghost atoms
----------------------------------------

Sometimes you want to update the internals (bonds, angles, torsions) when
one or more of the atoms involved are ghosts in either the reference or
perturbed states.

The function :func:`sire.morph.zero_ghost_torsions` will automatically add
torsion perturbations that zero the force constant of any torsions that
involve ghost atoms in that state. This is useful when you want to
fade in or out torsion forces as ghost atoms appear or disappear.

>>> mols = sr.morph.zero_ghost_torsions(mols)
>>> print(mols[0].perturbation().to_openmm().changed_torsions())
                  torsion       k0   k1  periodicity0  periodicity1  phase0   phase1
0   C5:5-C2:2-C3:3-H11:11  0.66944  0.0             3             3    -0.0    -0.0
1   C1:1-C2:2-C3:3-H11:11  0.66944  0.0             3             3    -0.0    -0.0
2   C5:5-C2:2-C3:3-H10:10  0.66944  0.0             3             3    -0.0    -0.0
3   C3:3-C2:2-C4:4-H12:12  0.66944  0.0             3             3    -0.0    -0.0
4   C1:1-C2:2-C3:3-H10:10  0.66944  0.0             3             3    -0.0    -0.0
5   C1:1-C2:2-C5:5-H16:16  0.66944  0.0             3             3    -0.0    -0.0
6     C5:5-C2:2-C3:3-H9:9  0.66944  0.0             3             3    -0.0    -0.0
7     C1:1-C2:2-C3:3-H9:9  0.66944  0.0             3             3    -0.0    -0.0
8     C3:3-C2:2-C1:1-H7:7  0.66944  0.0             3             3    -0.0    -0.0
9     C3:3-C2:2-C1:1-H8:8  0.66944  0.0             3             3    -0.0    -0.0
10    C3:3-C2:2-C1:1-H6:6  0.66944  0.0             3             3    -0.0    -0.0
11  C1:1-C2:2-C5:5-H15:15  0.66944  0.0             3             3    -0.0    -0.0
12  C4:4-C2:2-C5:5-H17:17  0.66944  0.0             3             3    -0.0    -0.0
13  C4:4-C2:2-C3:3-H11:11  0.66944  0.0             3             3    -0.0    -0.0
14  C4:4-C2:2-C3:3-H10:10  0.66944  0.0             3             3    -0.0    -0.0
15  C4:4-C2:2-C5:5-H16:16  0.66944  0.0             3             3    -0.0    -0.0
16    C4:4-C2:2-C3:3-H9:9  0.66944  0.0             3             3    -0.0    -0.0
17  C5:5-C2:2-C4:4-H14:14  0.66944  0.0             3             3    -0.0    -0.0
18  C5:5-C2:2-C4:4-H13:13  0.66944  0.0             3             3    -0.0    -0.0
19  C1:1-C2:2-C4:4-H14:14  0.66944  0.0             3             3    -0.0    -0.0
20  C4:4-C2:2-C5:5-H15:15  0.66944  0.0             3             3    -0.0    -0.0
21  C3:3-C2:2-C5:5-H17:17  0.66944  0.0             3             3    -0.0    -0.0
22  C1:1-C2:2-C4:4-H13:13  0.66944  0.0             3             3    -0.0    -0.0
23  C5:5-C2:2-C4:4-H12:12  0.66944  0.0             3             3    -0.0    -0.0
24  C1:1-C2:2-C4:4-H12:12  0.66944  0.0             3             3    -0.0    -0.0
25  C3:3-C2:2-C5:5-H16:16  0.66944  0.0             3             3    -0.0    -0.0
26    C5:5-C2:2-C1:1-H7:7  0.66944  0.0             3             3    -0.0    -0.0
27    C5:5-C2:2-C1:1-H8:8  0.66944  0.0             3             3    -0.0    -0.0
28    C5:5-C2:2-C1:1-H6:6  0.66944  0.0             3             3    -0.0    -0.0
29  C3:3-C2:2-C5:5-H15:15  0.66944  0.0             3             3    -0.0    -0.0
30  C3:3-C2:2-C4:4-H14:14  0.66944  0.0             3             3    -0.0    -0.0
31  C3:3-C2:2-C4:4-H13:13  0.66944  0.0             3             3    -0.0    -0.0
32    C4:4-C2:2-C1:1-H7:7  0.66944  0.0             3             3    -0.0    -0.0
33    C4:4-C2:2-C1:1-H8:8  0.66944  0.0             3             3    -0.0    -0.0
34    C4:4-C2:2-C1:1-H6:6  0.66944  0.0             3             3    -0.0    -0.0
35  C1:1-C2:2-C5:5-H17:17  0.66944  0.0             3             3    -0.0    -0.0

Similarly, the :func:`sire.morph.shrink_ghost_atoms` function will automatically
adjust the bond lengths of bonds that involve ghost atoms, so that they will
either be pulled into, or emerge from their connected atoms.

>>> mols = sr.morph.shrink_ghost_atoms(mols)
>>> print(mols[0].perturbation().to_openmm(constraint="bonds").changed_constraints())
       atompair  length0  length1
0     C2:2-C3:3  0.15375  0.06000
1   C5:5-H17:17  0.10969  0.06000
2     C2:2-C4:4  0.15375  0.10969
3     C1:1-C2:2  0.15375  0.06000
4     C1:1-H7:7  0.10969  0.06000
5     C2:2-C5:5  0.15375  0.06000
6     C1:1-H8:8  0.10969  0.06000
7     C1:1-H6:6  0.10969  0.06000
8   C3:3-H11:11  0.10969  0.06000
9   C5:5-H15:15  0.10969  0.06000
10    C3:3-H9:9  0.10969  0.06000
11  C5:5-H16:16  0.10969  0.06000
12  C3:3-H10:10  0.10969  0.06000

.. note::

   You can control the length of the ghost bond using the ``length``
   argument, e.g. ``shrink_ghost_atoms(mols, length="0.2A")`` would
   shrink the bond to 0.2 Å. The default length is 0.6 Å.

.. note::

   In general, you don't often need to pull ghost atoms into or out
   from their connected atoms. This is because a soft-core potential
   is used to soften interactions involving ghost atoms, such that
   they fade away smoothly as they disappear.
