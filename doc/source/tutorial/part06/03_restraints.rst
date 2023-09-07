==========
Restraints
==========

It can be useful to add restraints to the system to hold things in place.  For
example, you might want to hold the protein fixed while allowing the ligand to
move around. Or you may want to hold a ligand in place while it is being
decoupled from the simulation. Or you may want to use restraints to control
the distances between atoms or pull a molecule into a particular conformation.

You can specify the restraints you want to add to a system via the functions
in the :mod:`sire.restraints` module. These functions return
:class:`~sire.mm.Restraints` objects that contain all of the information that
needs to be passed to OpenMM to add the restraints to the system.

Positional Restraints
---------------------

The simplest type of restraint is a positional restraint, which holds a set of
atoms at a fixed position.  You can create a positional restraint using the
:func:`~sire.restraints.positional` function.  This function takes the
set of molecules or system you want to simulate, together with a search
string, list of atom indexes, or molecule views holding the atoms that
you want to restrain. For example,

>>> import sire as sr
>>> mols = sr.load(sr.expand(sr.tutorial_url, "ala.top", "ala.crd"))
>>> restraints = sr.restraints.positional(mols, "resname ALA and element C")
>>> print(restraints)
PositionalRestraints( size=3
0: PositionalRestraint( 8 => ( 16.5371, 5.02707, 15.812 ), k=150 kcal mol-1 Å-2 : r0=0 Å )
1: PositionalRestraint( 10 => ( 16.0464, 6.38937, 15.2588 ), k=150 kcal mol-1 Å-2 : r0=0 Å )
2: PositionalRestraint( 14 => ( 15.3698, 4.19397, 16.434 ), k=150 kcal mol-1 Å-2 : r0=0 Å )
)

will create positional restraints for the three carbon atoms in the alanine
residue. These carbon atoms are at indexes 8, 10 and 14 in ``mols.atoms()``,
e.g.

>>> print(mols.atoms([8, 10, 14]))
SireMol::SelectorM<SireMol::Atom>( size=3
0: MolNum(6) Atom( CA:9    [  16.54,    5.03,   15.81] )
1: MolNum(6) Atom( CB:11   [  16.05,    6.39,   15.26] )
2: MolNum(6) Atom( C:15    [  15.37,    4.19,   16.43] )
)

The default force constant is 150 kcal mol-1 Å-2. By default, the
atoms are restrained to their current positions, and are held exactly in
those positions (hence why the ``r0=0 Å`` in the output above).

You can set the force constant and width of the half-harmonic restraint
used to hold the atoms in position using the ``k`` and ``r0`` arguments, e.g.

>>> restraints = sr.restraints.positional(mols, "resname ALA and element C",
...                                       k="100 kcal mol-1 A-2", r0="1.5 A")
>>> print(restraints)
PositionalRestraints( size=3
0: PositionalRestraint( 8 => ( 16.5371, 5.02707, 15.812 ), k=100 kcal mol-1 Å-2 : r0=1.5 Å )
1: PositionalRestraint( 10 => ( 16.0464, 6.38937, 15.2588 ), k=100 kcal mol-1 Å-2 : r0=1.5 Å )
2: PositionalRestraint( 14 => ( 15.3698, 4.19397, 16.434 ), k=100 kcal mol-1 Å-2 : r0=1.5 Å )
)

By default, the atoms are restrained to their current positions. You can
specify a different position using the ``position`` argument, e.g.

>>> restraints = sr.restraints.positional(mols, "resname ALA and element C",
...                                       position=[(0,0,0), (1,1,1), (2,2,2)])
>>> print(restraints)
PositionalRestraints( size=3
0: PositionalRestraint( 8 => ( 0, 0, 0 ), k=150 kcal mol-1 Å-2 : r0=0 Å )
1: PositionalRestraint( 10 => ( 1, 1, 1 ), k=150 kcal mol-1 Å-2 : r0=0 Å )
2: PositionalRestraint( 14 => ( 2, 2, 2 ), k=150 kcal mol-1 Å-2 : r0=0 Å )
)

.. note::

   The number of positions must equal the number of atoms to be restrained,
   or equal to 1. If you only pass in a single position then this will be used
   for all of the restraints.

Sometimes it can be useful to use the same position for all of the restraints,
e.g. in the case of a spherical half-harmonic restraint that holds molecules
within a spherical bubble, e.g.

>>> center = mols[0].coordinates()
>>> restraints = sr.restraints.positional(mols,
...                                       mols[f"atoms within 10 of {center}"],
...                                       position=center,
...                                       r0="10 A")
>>> print(restraints)
PositionalRestraints( size=358
0: PositionalRestraint( 0 => ( 16.5471, 4.50102, 15.6589 ), k=150 kcal mol-1 Å-2 : r0=10 Å )
1: PositionalRestraint( 1 => ( 16.5471, 4.50102, 15.6589 ), k=150 kcal mol-1 Å-2 : r0=10 Å )
2: PositionalRestraint( 2 => ( 16.5471, 4.50102, 15.6589 ), k=150 kcal mol-1 Å-2 : r0=10 Å )
3: PositionalRestraint( 3 => ( 16.5471, 4.50102, 15.6589 ), k=150 kcal mol-1 Å-2 : r0=10 Å )
4: PositionalRestraint( 4 => ( 16.5471, 4.50102, 15.6589 ), k=150 kcal mol-1 Å-2 : r0=10 Å )
...
353: PositionalRestraint( 1892 => ( 16.5471, 4.50102, 15.6589 ), k=150 kcal mol-1 Å-2 : r0=10 Å )
354: PositionalRestraint( 1893 => ( 16.5471, 4.50102, 15.6589 ), k=150 kcal mol-1 Å-2 : r0=10 Å )
355: PositionalRestraint( 1906 => ( 16.5471, 4.50102, 15.6589 ), k=150 kcal mol-1 Å-2 : r0=10 Å )
356: PositionalRestraint( 1907 => ( 16.5471, 4.50102, 15.6589 ), k=150 kcal mol-1 Å-2 : r0=10 Å )
357: PositionalRestraint( 1908 => ( 16.5471, 4.50102, 15.6589 ), k=150 kcal mol-1 Å-2 : r0=10 Å )
)

has created half-harmonic restraints that affect all atoms within 10 Å of the
center of the first molecule, and that restrain those atoms to be within
a sphere of radius 10 Å centered on that molecule.

>>> print(mols.atoms()[[restraint.atom() for restraint in restraints]])
SireMol::SelectorM<SireMol::Atom>( size=358
0: MolNum(6) Atom( HH31:1  [  18.45,    3.49,   12.44] )
1: MolNum(6) Atom( CH3:2   [  18.98,    3.45,   13.39] )
2: MolNum(6) Atom( HH32:3  [  20.05,    3.63,   13.29] )
3: MolNum(6) Atom( HH33:4  [  18.80,    2.43,   13.73] )
4: MolNum(6) Atom( C:5     [  18.48,    4.55,   14.35] )
...
353: MolNum(630) Atom( H1:1893 [  14.92,    1.28,   17.05] )
354: MolNum(630) Atom( H2:1894 [  15.19,   -0.21,   17.07] )
355: MolNum(623) Atom( O:1907  [  21.65,    7.88,    9.79] )
356: MolNum(623) Atom( H1:1908 [  22.33,    8.56,    9.83] )
357: MolNum(623) Atom( H2:1909 [  21.07,    8.08,   10.53] )
)

Distance or Bond Restraints
---------------------------

The :func:`~sire.restraints.distance` or :func:`sire.restraints.bond` functions
are used to create bond or distance restraints. These are identical restraints,
so the functions are just synonyms of each other (they are the same function
with different names).

These functions take the set of molecules or system you want to simulate,
plus two search strings, lists of atom indexes, or molecule views holding
the pairs of atoms that you want to restrain.
