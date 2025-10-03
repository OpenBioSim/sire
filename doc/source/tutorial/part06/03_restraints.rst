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

.. note::

   For all harmonic restraints, the restraint energy is defined as ``E = k * x ** 2``,
   and not as ``E = 0.5 * k * x ** 2``. Hence, the force is ``F = 2 * k * x`` and all
   ``k`` values supplied to the restraints functions are half the value of the force
   constants.

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
>>> restraints = sr.restraints.positional(mols, atoms="resname ALA and element C")
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

You could have specified these indicies directly, e.g.

>>> restraints = sr.restraints.positional(mols, atoms=[8, 10, 14])
>>> print(restraints)
PositionalRestraints( size=3
0: PositionalRestraint( 8 => ( 16.5371, 5.02707, 15.812 ), k=150 kcal mol-1 Å-2 : r0=0 Å )
1: PositionalRestraint( 10 => ( 16.0464, 6.38937, 15.2588 ), k=150 kcal mol-1 Å-2 : r0=0 Å )
2: PositionalRestraint( 14 => ( 15.3698, 4.19397, 16.434 ), k=150 kcal mol-1 Å-2 : r0=0 Å )
)

or just passed in the atoms directly, e.g.

>>> atoms = mols.atoms([8, 10, 14])
>>> restraints = sr.restraints.positional(mols, atoms=atoms)
>>> print(restraints)
PositionalRestraints( size=3
0: PositionalRestraint( 8 => ( 16.5371, 5.02707, 15.812 ), k=150 kcal mol-1 Å-2 : r0=0 Å )
1: PositionalRestraint( 10 => ( 16.0464, 6.38937, 15.2588 ), k=150 kcal mol-1 Å-2 : r0=0 Å )
2: PositionalRestraint( 14 => ( 15.3698, 4.19397, 16.434 ), k=150 kcal mol-1 Å-2 : r0=0 Å )
)

The default half force constant is 150 kcal mol-1 Å-2. By default, the
atoms are restrained to their current positions, and are held exactly in
those positions (hence why the ``r0=0 Å`` in the output above).

You can set the half force constant and width of the half-harmonic restraint
used to hold the atoms in position using the ``k`` and ``r0`` arguments, e.g.

>>> restraints = sr.restraints.positional(mols, atoms="resname ALA and element C",
...                                       k="100 kcal mol-1 A-2", r0="1.5 A")
>>> print(restraints)
PositionalRestraints( size=3
0: PositionalRestraint( 8 => ( 16.5371, 5.02707, 15.812 ), k=100 kcal mol-1 Å-2 : r0=1.5 Å )
1: PositionalRestraint( 10 => ( 16.0464, 6.38937, 15.2588 ), k=100 kcal mol-1 Å-2 : r0=1.5 Å )
2: PositionalRestraint( 14 => ( 15.3698, 4.19397, 16.434 ), k=100 kcal mol-1 Å-2 : r0=1.5 Å )
)

By default, the atoms are restrained to their current positions. You can
specify a different position using the ``position`` argument, e.g.

>>> restraints = sr.restraints.positional(mols, atoms="resname ALA and element C",
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
...                                       atoms=mols[f"atoms within 10 of {center}"],
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

Note that the restraints only contain the indicies of the atoms that are
restrained, not the atoms themselves. You can get the atoms by looking up
those indicies from the ``mols`` molecular container from which the
atoms were selected. For example, here we get all of the atoms that are
subject to that spherical half-harmonic restraint;

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

RMSD Restraints
---------------------

RMSD restraints are similar to positional restraints, but rather than atoms being 
restrained with respect to fixed positions in space, they are restrained relative 
to a reference conformation, which can be defined by the user. 
The function :func:`sire.restraints.rmsd` takes a set of atoms as the input, 
which can be inputted in an identical manner to :func:`~sire.restraints.positional`.
The default reference conformation used to calculate the RMSD is the current state
of the sire :class:`~sire.system.System` object to which the restraints are applied to.
Alternatively, a custom reference system may be provided using the ``ref`` argument.
Note that ``ref`` must have the same number of particles as the input system. For example,

>>> restraints = sr.restraints.rmsd(mols, atoms="atomname CA, C, N")
>>> print(restraints)
RMSDRestraints( name=restraint, size=1
0: RMSDRestraint( [4, 6, 8, 14, 16], k=150 kcal mol-1 Å-2, r0=0 Å )
)

This creates RMSD restraints on the backbone atoms of the system.

>>> print(mols.atoms("atomname CA, C, N"))
SireMol::SelectorM<SireMol::Atom>( size=5
0: MolNum(652) Atom( C:5     [  18.48,    4.55,   14.35] )
1: MolNum(652) Atom( N:7     [  17.22,    4.31,   14.71] )
2: MolNum(652) Atom( CA:9    [  16.54,    5.03,   15.81] )
3: MolNum(652) Atom( C:15    [  15.37,    4.19,   16.43] )
4: MolNum(652) Atom( N:17    [  14.97,    4.59,   17.58] )
)

RMSD restraints can be generated by using the energy-minimised structure 
as the reference conformation, e.g.

>>> restraints = sr.restraints.rmsd(mols, atoms=[4, 6, 8, 14, 16], ref=mols.minimisation().run().commit())

By default, the force constant is set to 150 kcal mol-1 Å-2 and the width of
the flat-bottom potential well is set to 0 Å. 
This can be changed by providing custom ``k`` and ``r0`` values, e.g.

>>> restraints = sr.restraints.rmsd(mols, atoms="atomname CA, C, N",
...                                       k="100 kcal mol-1 A-2", r0="1.5 A", name='rmsd_restraints')
>>> print(restraints)
RMSDRestraints( name=rmsd_restraints, size=1
0: RMSDRestraint( [4, 6, 8, 14, 16], k=100 kcal mol-1 Å-2, r0=1.5 Å )
)

Distance or Bond Restraints
---------------------------

The :func:`sire.restraints.distance` or :func:`sire.restraints.bond` functions
are used to create bond or distance restraints. These are identical restraints,
so the functions are just synonyms of each other (they are the same function
with different names).

These functions take the set of molecules or system you want to simulate,
plus two search strings, lists of atom indexes, or molecule views holding
the pairs of atoms that you want to restrain, e.g.

>>> restraints = sr.restraints.distance(mols, atoms0=0, atoms1=1)
>>> print(restraints)
BondRestraints( size=1
0: BondRestraint( 0 <=> 1, k=150 kcal mol-1 Å-2 : r0=1.09 Å )
)

or

>>> restraints = sr.restraints.bond(mols, atoms0=0, atoms1=1)
>>> print(restraints)
BondRestraints( size=1
0: BondRestraint( 0 <=> 1, k=150 kcal mol-1 Å-2 : r0=1.09 Å )
)


creates a single harmonic distance (or bond) restraint that acts between
atoms 0 and 1. By default, the equilibrium distance (r0) is the current
distance between the atoms (1.09 Å), and the half force constant (k) is
150 kcal mol-1 Å-2.

You can set these via the ``k`` and ``r0`` arguments, e.g.

>>> restraints = sr.restraints.bond(mols, atoms0=0, atoms1=1,
...                                 k="100 kcal mol-1 A-2", r0="1.5 A")
>>> print(restraints)
BondRestraints( size=1
0: BondRestraint( 0 <=> 1, k=100 kcal mol-1 Å-2 : r0=1.5 Å )
)

You can specify as many atoms pairs as you like, e.g.

>>> restraints = sr.restraints.bond(mols, atoms0=[0, 1, 2], atoms1=[3, 4, 5])
>>> print(restraints)
BondRestraints( size=3
0: BondRestraint( 0 <=> 3, k=150 kcal mol-1 Å-2 : r0=1.71245 Å )
1: BondRestraint( 1 <=> 4, k=150 kcal mol-1 Å-2 : r0=1.54643 Å )
2: BondRestraint( 2 <=> 5, k=150 kcal mol-1 Å-2 : r0=2.48642 Å )
)

You can specify the atoms using a search string, passing the atoms themselves,
or using the atom index, just as you could for the positional restraints.

>>> restraints = sr.restraints.bond(mols,
...                                 atoms0=mols[0][0],
...                                 atoms1=mols["water and element O"])
>>> print(restraints)
BondRestraints( size=630
0: BondRestraint( 0 <=> 22, k=150 kcal mol-1 Å-2 : r0=13.2847 Å )
1: BondRestraint( 0 <=> 25, k=150 kcal mol-1 Å-2 : r0=10.8445 Å )
2: BondRestraint( 0 <=> 28, k=150 kcal mol-1 Å-2 : r0=15.9183 Å )
3: BondRestraint( 0 <=> 31, k=150 kcal mol-1 Å-2 : r0=13.5108 Å )
4: BondRestraint( 0 <=> 34, k=150 kcal mol-1 Å-2 : r0=9.18423 Å )
...
625: BondRestraint( 0 <=> 1897, k=150 kcal mol-1 Å-2 : r0=9.52934 Å )
626: BondRestraint( 0 <=> 1900, k=150 kcal mol-1 Å-2 : r0=12.7207 Å )
627: BondRestraint( 0 <=> 1903, k=150 kcal mol-1 Å-2 : r0=10.8053 Å )
628: BondRestraint( 0 <=> 1906, k=150 kcal mol-1 Å-2 : r0=6.04142 Å )
629: BondRestraint( 0 <=> 1909, k=150 kcal mol-1 Å-2 : r0=17.1035 Å )
)

.. note::

   If the number of atoms in ``atoms0`` and ``atoms1`` are not equal, then
   the smaller list will be extended to match by appending multiple copies
   of the last atom. In this case, this duplicates the single atom in
   ``atoms0``, meaning that this has created distance restraints between
   the first atom in the first molecule and the oxygen atoms in all of
   the water molecules.

Angle or Dihedral Restraints
---------------------------

The :func:`sire.restraints.angle` or :func:`sire.restraints.dihedral` functions
are used to create angle or distance restraints.

Just like the other restraint functions, these functions take 
the set of molecules or system you want to simulate,
plus a search string, lists of atom indexes, or molecule views
holding the atoms that you want to restrain., e.g.

>>> restraints = sr.restraints.angle(mols=mols, atoms=[0, 1, 2])
>>> print(restraints)
AngleRestraints( name=restraint, size=1
0: AngleRestraint( [0, 1, 2], theta0=112.8°, ktheta=0.0304617 kcal mol-1 °-2 )
 )

or

>>> restraints = sr.restraints.dihedral(mols=mols, atoms=[0, 1, 2, 3])
>>> print(restraints)
DihedralRestraints( name=restraint, size=1
0: DihedralRestraint( [0, 1, 2, 3], phi0=244.528°, kphi=0.0304617 kcal mol-1 °-2 )
 )

creates a single harmonic angle or dihedral restraint that acts between
the specified atoms. By default, the equilibrium angle (theta0 or phi0)
is the current angle between the atoms (112.8° or 244.528°), 
and the force constant (ktheta or kphi) is 100 kcal mol-1 rad-2.

You can set these via the ``ktheta`` or ``kphi`` and ``theta0`` or ``phi0`` arguments depending
on the restraint used, e.g.

>>> restraints = sr.restraints.angle(mols=mols, atoms=[0, 1, 2],
...                                  theta0="1.5 rad", ktheta="50 kcal mol-1 rad-2")
>>> print(restraints)
AngleRestraints( name=restraint, size=1
0: AngleRestraint( [0, 1, 2], theta0=85.9437°, ktheta=0.0152309 kcal mol-1 °-2 )
 )

or

>>> restraints = sr.restraints.dihedral(mols=mols, atoms=[0, 1, 2, 3],
...                                     phi0="2 rad", kphi="10 kcal mol-1 rad-2")
>>> print(restraints)
DihedralRestraints( name=restraint, size=1
0: DihedralRestraint( [0, 1, 2, 3], phi0=114.592°, kphi=0.00304617 kcal mol-1 °-2 )
 )

You can specify the atoms using a search string, passing the atoms themselves,
or passing the atoms from a molecular container.

>>> ang = mols.angles()[0]
>>> restraints = sr.restraints.angle(mols=mols, atoms=ang.atoms())
>>> print(restraints)
AngleRestraints( name=restraint, size=1
0: AngleRestraint( [0, 1, 2], theta0=112.8°, ktheta=0.0304617 kcal mol-1 °-2 )
 )

or

>>> dih = mols.dihedrals()[0]
>>> restraints = sr.restraints.dihedral(mols=mols, atoms=dih.atoms())
>>> print(restraints)
DihedralRestraints( name=restraint, size=1
0: DihedralRestraint( [0, 1, 4, 5], phi0=243.281°, kphi=0.0304617 kcal mol-1 °-2 )
 )

Morse Potential Restraints
---------------------------

The :func:`sire.restraints.morse_potential` function is used to create Morse potential restraints,
which can be used to carry harmonic bond annihilations alchemical relative binding free energy calculations.

To create a Morse potential restraint, you need to specify the two atoms to be restrained. Like the distance restraints,
the atoms can be specified using a search string, passing lists of atom indexes, or 
molecule views holding the atoms. You have to specify the bond force constants,
equilibrium bond distance value and the dissociation energy for the restraints.
If not supplied, automatic parametrisation feature can be used, which will detect the bond being alchemically
annihilated and set the parameters accordingly (dissociation energy value still needs to be provided). For example,

>>> mols = sr.load_test_files("cyclopentane_cyclohexane.bss")
>>> morse_restraints = sr.restraints.morse_potential(mols,
                                                   atoms0=mols["molecule property is_perturbable and atomidx 0"],
                                                   atoms1=mols["molecule property is_perturbable and atomidx 4"],
                                                   k="100 kcal mol-1 A-2",
                                                   r0="1.5 A",
                                                   de="50 kcal mol-1")
>>> morse_restraint = morse_restraints[0]
>>> print(morse_restraint)
MorsePotentialRestraint( 0 <=> 4, k=100 kcal mol-1 Å-2 : r0=1.5 Å : de=50 kcal mol-1 )

creates a Morse potential restraint between atoms 0 and 4 using the specified parameters.
Alternatively, if the molecule contains a bond that is being alchemically annihilated, e.g.

>>> morse_restraints = sr.restraints.morse_potential(mols, auto_parametrise=True, de="50 kcal mol-1")
>>> morse_restraint = morse_restraints[0]
>>> print(morse_restraint)
MorsePotentialRestraint( 0 <=> 4, k=228.89 kcal mol-1 Å-2 : r0=1.5354 Å : de=50 kcal mol-1 )

creates a Morse potential restraint between atoms 0 and 4 by attempting to match the bond parameters of the bond being alchemically annihilated.

Boresch Restraints
---------------------------

The :func:`sire.restraints.boresch` function is used to create Boresch restraints,
which are commonly used for restraining the relative positions and orientations
of a receptor and ligand during alchemical absolute binding free energy calculations.
They restrain the six relative external degrees of freedom of the receptor and ligand
by restraining one distance, two angles, and three dihedrals angles which are
defined according to three anchor points in the receptor and three anchor points
in the ligand. For more detail, please see J. Phys. Chem. B 2003, 107, 35, 9535–9551. 

To create a Boresch restraint, you need to specify the receptor and ligand anchor
atoms (note that the order of the atoms is important). Like the distance restraints,
the atoms can be specified using a search string, passing lists of atom indexes, or 
molecule views holding the atoms. You can also specify the half force constants and equilibrium
values for the restraints. If not supplied, default half force constants of 5 kcal mol-1 Å-2
and 50 kcal mol-1 rad-2 are used for the distance and angle restraints, respectively,
and the equilibrium values are set to the current values of the distances and angles in
the system supplied. For example,

>>> boresch_restraints = sr.restraints.boresch(mols, receptor=[1574, 1554, 1576], ligand=[4,3,5])
>>> boresch_restraint = boresch_restraints[0]
>>> print(boresch_restraint)
BoreschRestraint( [1574, 1554, 1576] => [4, 3, 5],
                  k=[5 kcal mol-1 Å-2, 0.0152309 kcal mol-1 °-2, 0.0152309 kcal mol-1 °-2, 
                  0.0152309 kcal mol-1 °-2, 0.0152309 kcal mol-1 °-2, 0.0152309 kcal mol-1 °-2]
                  r0=15.1197 Å, theta0=[80.5212°, 59.818°],
                  phi0=[170.562°Ⱐ128.435°Ⱐ192.21°] )

creates a Boresch restraint where the receptor anchor atoms are r1 = 1574, r2 = 1554, and r3 = 1576,
and the ligand anchor atoms are l1 = 4, l2 = 3, and l3 = 5. The default half force constants have been set
and the equilibrium values have been set to the current values of the distances and angles in the
system supplied.

.. note::

   Boresch restraints can be subject to instabilities if any three contiguous anchor points
   approach collinearity (J. Chem. Theory Comput. 2023, 19, 12, 3686–3704). It is important to
   prevent this by ensuring the associated angles are sufficiently far from 0 or 180 degrees,
   and that the `ktheta` half force constants are high enough. Sire will raise a warning if the 
   `theta0` values are too close to 0 or 180 degrees for the given temperature and half force constants.

Alternatively, we could have explicitly set the half force constants and equilibrium values, e.g.

>>> boresch_restraints = sr.restraints.boresch(mols, receptor=[1574, 1554, 1576], ligand=[4,3,5],
                                               kr = "6.2012 kcal mol-1 A-2",
                                               ktheta = ["28.7685 kcal mol-1 rad-2", "24.8204 kcal mol-1 rad-2"],
                                               kphi = ["59.8626 kcal mol-1 rad-2", "0.7923 kcal mol-1 rad-2", 
                                                       "55.1775 kcal mol-1 rad-2"],
                                               r0 = "16 A", 
                                               theta0 = ["1.2 rad", "1.3 rad"],
                                               phi0 = ["2.2 rad", "2.5 rad", "1.5 rad"],
                                               name = "boresch_restraint_1")
>>> boresch_restraint = boresch_restraints[0]
>>> print(boresch_restraint)
BoreschRestraint( [1574, 1554, 1576] => [4, 3, 5],
                  k=[6.2012 kcal mol-1 Å-2, 0.00876339 kcal mol-1 °-2, 0.00756073 kcal mol-1 °-2, 0.0182352 kcal mol-1 °-2, 0.000241348 kcal mol-1 °-2, 0.016808 kcal mol-1 °-2]
                  r0=16 Å, theta0=[68.7549°, 74.4845°],
                  phi0=[126.051°Ⱐ143.239°Ⱐ85.9437°] )

.. note::

   :func:`sire.restraints.boresch` returns a list of Boresch restraints. If you are only
   interested in a single Boresch restraint, you can extract it with the index, e.g.
   ``boresch_restraint = boresch_restraints[0]``.

When performing an alchemical absolute binding free energy calculation, it is necessary to 
calculate the free energy of releasing the decoupled ligand to the standard state volume.
The analytical Boresch correction is almost always accurate if stable restraints have been
selected (see 10.26434/chemrxiv-2023-8s9dz-v2). This can be calculated with
:func:`sire.restraints.standard_state_correction`:

>>> boresch_restraint = boresch_restraints[0]
>>> from sire import u
>>> correction = sr.restraints.standard_state_correction(boresch_restraint, temperature=u("298 K"))
>>> print(correction)
-6.2399 kcal mol-1


Using restraints in minimisation or dynamics
--------------------------------------------

You can use restraints in a minimisation or dynamics simulation by
passing them in via the ``restraints`` argument, e.g.

>>> restraints = sr.restraints.positional(mols, atoms="resname ALA and element C")
>>> d = mols.dynamics(timestep="4fs", temperature="25oC",
...                   restraints=restraints)
>>> d.run("10ps")
>>> mols = d.commit()

runs a dynamics simulation using positional restraints on the carbon atoms
of the ``ALA`` residue, while

>>> restraints = sr.restraints.distance(mols, atoms0=mols[0][0],
...                                     atoms1=mols[1][0], r0="5 A")
>>> mols = mols.minimisation(restraints=restraints).run().commit()

performs a minimisation where a distance restraint is applied between the
first atoms of the first two molecules, pulling them to be 5 Å apart.

You can pass in a list of restraints, meaning that you can use as
many as you wish during a simulation.

>>> pos_rests = sr.restraints.positional(mols, atoms="resname ALA and element C")
>>> dst_rests = sr.restraints.distance(mols, atoms0=mols[0][0], atoms1=mols[1][0])
>>> mols = mols.minimisation(restraints=[pos_rests, dst_rests]).run().commit()
>>> d = mols.dynamics(timestep="4fs", temperature="25oC",
...                   restraints=[pos_rests, dst_rests])
>>> d.run("10ps")
>>> mols = d.commit()

.. note::

   It is a good idea to run minimisation before dynamics whenever you add
   restraints to a system. This is because the restraints could put a lot
   of energy into the system, causing blow-ups or ``NaN`` errors.

Using restraints in the low-level API
-------------------------------------

You can use restraints in the low-level API by passing them in as a
``map`` parameter using the key ``restraints``, e.g.

>>> omm = sr.convert.to(mols, "openmm",
...                     map={"restraints": [pos_rests, dst_rests]})
>>> print(omm)
openmm::Context( num_atoms=1915 integrator=VerletIntegrator timestep=1.0 fs platform=HIP )

More information about the mappable options for the low-level API can be
found in the :doc:`detailed guide <../../cheatsheet/openmm>`.

Fixing atoms in place
---------------------

As well as restraining atoms, you also have the option of fixing atoms in
space during the simulation. Atoms which are fixed are not moved by
minimisation or the dynamics integrator. This can be useful if you, e.g.
want to freeze atoms outside a radius of a ligand binding site, or if you
want to minimise the solvent while keeping the solute fixed.

You can fix atoms by passing in the set of atoms to fix as the
``fixed`` argument to the :meth:`~sire.mol.SelectorMol.minimisation` or
:meth:`~sire.mol.SelectorMol.dynamics` functions, e.g.

>>> mols = sr.load(sr.expand(sr.tutorial_url, "ala.top", "ala.crd"))
>>> d = mols.dynamics(timestep="4fs", temperature="25oC",
...                   fixed=mols[0])
>>> d.run("10ps")
>>> mols = d.commit()

will run dynamics where all of the atoms in ``mols[0]`` (the solute)
are fixed. You can pass in a molecular container containing all of the
atoms that should be fixed, or a search string, or the atom indexes
of the atoms. Here,

>>> d = mols.dynamics(timestep="1fs", temperature="25oC",
...                   fixed="element C")
>>> d.run("10ps")
>>> mols = d.commit()

we have run dynamics where all of the carbon atoms are fixed, while

>>> d = mols.dynamics(timestep="1fs", temperature="25oC",
...                   fixed=[0, 2, 4, 6, 8])
>>> d.run("10ps")
>>> mols = d.commit()

would run dynamics where atoms at indicies 0, 2, 4, 6 and 8 are fixed.

You could even use search strings to fix atoms by distances. If you do this,
it is best to also add half-harmonic positional restraints that hold nearby
molecules within the sphere of mobile atoms, e.g.

>>> center = mols[0].coordinates()
>>> radius = "10 A"
>>> restraints = sr.restraints.positional(mols,
...                                       atoms=mols[f"molecules within 10 of {center}"],
...                                       position=center,
...                                       r0=radius)
>>> d = mols.dynamics(timestep="1fs", temperature="25oC",
...                   fixed=f"not molecules within {radius} of {center}")
>>> d.run("10ps")
>>> mols = d.commit()

This implements the restraints in OpenMM by setting the masses of the fixed
atoms to zero. This signals the OpenMM integrator to skip over these
atoms and not move them.

.. note::

   Be careful using constraints with fixed atoms. If a constraint involves
   a fixed atom and a mobile atom, then there is a high risk that the
   constraint won't be able to be satisfied during dynamics, and
   ``Particle coordinate is NaN`` error or similar will occur.
   It it safest to use a small timestep (e.g. 1 fs) when constraining
   only parts of molecules.

.. note::

   While the fixed atoms are not moved by the integrator, they are still
   included in the energy calculation. This means that the computational
   cost of the simulation will still be high, and also that energies
   will not be conserved. Make sure you use an integrator that provides
   a heat bath (e.g. NVT or NPT integrators).
