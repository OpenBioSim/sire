============================================
Searching for Fragments using Smarts Strings
============================================

In the :doc:`last page <03_smiles>` you saw how we can search for molecules
using their smiles string.
`Smarts strings <https://en.wikipedia.org/wiki/SMILES_arbitrary_target_specification>`__
provide a much more powerful way to search for molecules, or fragments within
molecules.

They can be vieweed as a "regular expression" for chemical structures.
`This blog post <https://russodanielp.github.io/blog/a-brief-introduction-to-smarts/>`__
gives an excellent introduction to what they are and how they can be used.

You can search for matches using a smarts string just by putting
``smarts`` before the query, e.g.

>>> mols = sr.load(sr.expand(sr.tutorial_url, "ala.crd", "ala.top"))
>>> print(mols["smarts [#6]"])
AtomMatch( size=6
0: [1] CH3:2
1: [1] C:5
2: [1] CA:9
3: [1] CB:11
4: [1] C:15
5: [1] CH3:19
)

will search for all carbon atoms. In this case there are six matches
to the query, which can be seen in the returned :class:`~sire.mol.AtomMatch`.

>>> m = mols["smarts [#6]"]
>>> print(m.num_groups())
6

You can extract individual matches using the
:func:`~sire.mol.AtomMatch.group` function, and all matches using
the :func:`~sire.mol.AtomMatch.groups` function.

>>> print(m.group(3))
Selector<SireMol::Atom>( size=1
0:  Atom( CB:11   [  16.05,    6.39,   15.26] )
)

>>> print(m.groups())
[Selector<SireMol::Atom>( size=1
0:  Atom( CH3:2   [  18.98,    3.45,   13.39] )
), Selector<SireMol::Atom>( size=1
0:  Atom( C:5     [  18.48,    4.55,   14.35] )
), Selector<SireMol::Atom>( size=1
0:  Atom( CA:9    [  16.54,    5.03,   15.81] )
), Selector<SireMol::Atom>( size=1
0:  Atom( CB:11   [  16.05,    6.39,   15.26] )
), Selector<SireMol::Atom>( size=1
0:  Atom( C:15    [  15.37,    4.19,   16.43] )
), Selector<SireMol::Atom>( size=1
0:  Atom( CH3:19  [  13.83,    3.94,   18.35] )
)]

The atoms contained in the :class:`~sire.mol.AtomMatch` object are all
the ones that are involved in any of the matches. Since this object is
also a standard molecule view, it can be used like any other
molecule view, e.g.

>>> print(m.atoms())
Selector<SireMol::Atom>( size=6
0:  Atom( CH3:2   [  18.98,    3.45,   13.39] )
1:  Atom( C:5     [  18.48,    4.55,   14.35] )
2:  Atom( CA:9    [  16.54,    5.03,   15.81] )
3:  Atom( CB:11   [  16.05,    6.39,   15.26] )
4:  Atom( C:15    [  15.37,    4.19,   16.43] )
5:  Atom( CH3:19  [  13.83,    3.94,   18.35] )
)

More complex searches
---------------------

Smarts searches can match multiple atoms in a group. For example,
``[#6]!:[#6]`` matches a carbon bonded via an aliphatic bond to another
carbon (i.e. it matches aliphatic carbons). There are two atoms
contained in this search, so each group match will contain two
atoms.

>>> m = mols["smarts [#6]!:[#6]"]
>>> print(m)
AtomMatch( size=3
0: [2] CH3:2,C:5
1: [2] CA:9,CB:11
2: [2] CA:9,C:15
)

The atoms in each group match are returned in the same order as they are
matched in the smarts string. They can be accessed via the standard index
operator of the atom container, e.g.

>>> print(m.group(0)[1])
Atom( C:5     [  18.48,    4.55,   14.35] )

prints the second matching atom in the first matched group.

Matching multiple molecules
---------------------------

A smarts match may match multiple molecules. For example,
``[#1]-[#8]-[#1]`` would match ``H-O-H``, in that order.

>>> m = mols["smarts [#1]-[#8]-[#1]"]
>>> print(m)
AtomMatchM( size=630
0: WAT:1965 [3] H1:24,O:23,H2:25
1: WAT:1974 [3] H1:27,O:26,H2:28
2: WAT:1984 [3] H1:30,O:29,H2:31
3: WAT:1993 [3] H1:33,O:32,H2:34
4: WAT:2002 [3] H1:36,O:35,H2:37
...
625: WAT:2395 [3] H1:1899,O:1898,H2:1900
626: WAT:2419 [3] H1:1902,O:1901,H2:1903
627: WAT:2430 [3] H1:1905,O:1904,H2:1906
628: WAT:2463 [3] H1:1908,O:1907,H2:1909
629: WAT:2475 [3] H1:1911,O:1910,H2:1912
)

Note here how a :class:`~sire.mol.AtomMatchM` object is returned, signifying
this contains multiple molecules. The same
:func:`~sire.mol.AtomMatchM.num_groups`,
:func:`~sire.mol.AtomMatchM.group` and
:func:`~sire.mol.AtomMatchM.groups`,  functions can be used to
extract individual match groups.

>>> print(m.group(5))
Selector<SireMol::Atom>( size=3
0:  Atom( H1:39   [  19.17,   10.39,   17.12] )
1:  Atom( O:38    [  18.83,    9.65,   16.62] )
2:  Atom( H2:40   [  19.39,    9.61,   15.84] )
)

Note again how the atoms are returned in the same order as they matched
the smarts string.

>>> print(mols[6].atoms())
Selector<SireMol::Atom>( size=3
0:  Atom( O:38    [  18.83,    9.65,   16.62] )
1:  Atom( H1:39   [  19.17,   10.39,   17.12] )
2:  Atom( H2:40   [  19.39,    9.61,   15.84] )
)

Generating smarts strings from molecules
----------------------------------------

You can generate a smarts string from a molecule by calling its
:func:`~sire.mol.SelectorMol.smarts` function, e.g.

>>> print(mols[0].smarts())
[#6H3]-[#6](=[#8])-[#7H]-[#6H](-[#6H3])-[#6](=[#8])-[#7H]-[#6H3]

You can include hydrogens by specifying ``include_hydrogens=True``

>>> print(mols[0].smarts(include_hydrogens=True))
[H]-[#6](-[H])(-[H])-[#6](=[#8])-[#7](-[H])-[#6](-[H])(-[#6](-[H])(-[H])-[H])-[#6](=[#8])-[#7](-[H])-[#6](-[H])(-[H])-[H]

You can get the smarts string as a :mod:`sire` search term by
passing ``as_search=True``, e.g.

>>> print(mols[0].smarts(as_search=True))
smarts [#6H3]-[#6](=[#8])-[#7H]-[#6H](-[#6H3])-[#6](=[#8])-[#7H]-[#6H3]

>>> print(mols[mols[0].smarts(as_search=True)])
AtomMatch( size=1
0: [10] CH3:2,C:5,O:6,N:7,CA:9...
)

Creating molecules from smarts
------------------------------

You can also create fragment molecules using smarts strings via
the :func:`sire.smarts` function. Note that
these can only be used for substructure searching, as they are lacking
properties such as coordinates.

>>> mol = sr.smarts("[#6H3]-[#6](=[#8])-[#7H]-[#6H](-[#6H3])-[#6](=[#8])-[#7H]-[#6H3]")
>>> print(mol.atoms())
Selector<SireMol::Atom>( size=10
0:  Atom( C1:1 )
1:  Atom( C2:2 )
2:  Atom( O3:3 )
3:  Atom( N4:4 )
4:  Atom( C5:5 )
5:  Atom( C6:6 )
6:  Atom( C7:7 )
7:  Atom( O8:8 )
8:  Atom( N9:9 )
9:  Atom( C10:10 )
)

>>> print(mol.bonds())
SelectorBond( size=9
0: Bond( C1:1 => C2:2 )
1: Bond( C2:2 => O3:3 )
2: Bond( C2:2 => N4:4 )
3: Bond( N4:4 => C5:5 )
4: Bond( C5:5 => C6:6 )
5: Bond( C5:5 => C7:7 )
6: Bond( C7:7 => O8:8 )
7: Bond( C7:7 => N9:9 )
8: Bond( N9:9 => C10:10 )
)

>>> print(mol.properties())
Properties(
    formal_charge => SireMol::AtomCharges( size=10
0: 0 |e|
1: 0 |e|
2: 0 |e|
3: 0 |e|
4: 0 |e|
5: 0 |e|
6: 0 |e|
7: 0 |e|
8: 0 |e|
9: 0 |e|
),
    hybridization => SireMol::AtomHybridizations( size=10
0: unspecified
1: unspecified
2: unspecified
3: unspecified
4: unspecified
5: unspecified
6: unspecified
7: unspecified
8: unspecified
9: unspecified
),
    element => SireMol::AtomElements( size=10
0: Carbon (C, 6)
1: Carbon (C, 6)
2: Oxygen (O, 8)
3: Nitrogen (N, 7)
4: Carbon (C, 6)
5: Carbon (C, 6)
6: Carbon (C, 6)
7: Oxygen (O, 8)
8: Nitrogen (N, 7)
9: Carbon (C, 6)
),
    isotope => SireMol::AtomIntProperty( size=10
0: 0
1: 0
2: 0
3: 0
4: 0
5: 0
6: 0
7: 0
8: 0
9: 0
),
    mass => SireMol::AtomMasses( size=10
0: 12.011 g mol-1
1: 12.011 g mol-1
2: 15.999 g mol-1
3: 14.007 g mol-1
4: 12.011 g mol-1
5: 12.011 g mol-1
6: 12.011 g mol-1
7: 15.999 g mol-1
8: 14.007 g mol-1
9: 12.011 g mol-1
),
    chirality => SireMol::AtomChiralities( size=10
0: unspecified
1: unspecified
2: unspecified
3: unspecified
4: unspecified
5: unspecified
6: unspecified
7: unspecified
8: unspecified
9: unspecified
),
    is_aromatic => SireMol::AtomIntProperty( size=10
0: 0
1: 0
2: 0
3: 0
4: 0
5: 0
6: 0
7: 0
8: 0
9: 0
),
    connectivity => Connectivity: nConnections() == 9.
Connected residues:
Connected atoms:
  * Atom C1:LIG:1 bonded to C2:LIG:1.
  * Atom C2:LIG:1 bonded to O3:LIG:1 N4:LIG:1 C1:LIG:1.
  * Atom O3:LIG:1 bonded to C2:LIG:1.
  * Atom N4:LIG:1 bonded to C2:LIG:1 C5:LIG:1.
  * Atom C5:LIG:1 bonded to N4:LIG:1 C7:LIG:1 C6:LIG:1.
  * Atom C6:LIG:1 bonded to C5:LIG:1.
  * Atom C7:LIG:1 bonded to N9:LIG:1 O8:LIG:1 C5:LIG:1.
  * Atom O8:LIG:1 bonded to C7:LIG:1.
  * Atom N9:LIG:1 bonded to C10:LIG:1 C7:LIG:1.
  * Atom C10:LIG:1 bonded to N9:LIG:1.
)

