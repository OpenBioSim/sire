=================
Residue mutations
=================

So far, perturbations have involved changing an entire molecule.
However, we can also change individual residues in a molecule. This
can be useful to, e.g. study the effect of residue mutation on a protein,
and how this could impact ligand binding or protein folding.

Sub-structure matching
----------------------

To start, we need to define the sub-structure that we want to change.

Let's first load the kigaki protein system.

>>> import sire as sr
>>> mols = sr.load_test_files("kigaki.gro", "kigaki.top")
>>> protein = mols[0]
>>> print(protein)
Molecule( Protein:6 num_atoms=302 num_residues=19 )
>>> print(protein.residues().names())
[ResName('LYS'), ResName('ILE'), ResName('GLY'), ResName('ALA'),
 ResName('LYS'), ResName('ILE'), ResName('LYS'), ResName('ILE'),
 ResName('GLY'), ResName('ALA'), ResName('LYS'), ResName('ILE'),
 ResName('LYS'), ResName('ILE'), ResName('GLY'), ResName('ALA'),
 ResName('LYS'), ResName('ILE'), ResName('NH2')]

To start with, let's mutate the first alanine residue that we find in the
protein into a lysine. First, let's select the first alanine...

>>> ala = protein["resname ALA"][0]
>>> print(ala)
Residue( ALA:4   num_atoms=10 )

To mutate this into a lysine, we need to tell :mod:`sire` what a
lysine looks like. Fortunately, our protein already contains lysine
residues, so let's select on of those...

>>> lys = protein["resname LYS"][1]
>>> print(lys)
Residue( LYS:5   num_atoms=22 )

.. note::

   We have selected the second lysine residue in the protein,
   as the first lysine residue is the N-terminal residue. Since
   our alanine is a mid-chain residue, it makes more sense to
   use a mid-chain lysine.

Next, we need to match the alanine residue to the lysine residue,
just as we did when perturbing entire molecules. Calling
:func:`sire.morph.match` on sub-views of molecules will only match
the atoms in those views.

>>> mapping = sr.morph.match(ala, lys, match_light_atoms=True)
>>> print(mapping)
AtomMapping( size=9, unmapped0=1, unmapped1=13
0: MolNum(6) Atom( HA:54 ) <=> MolNum(6) Atom( HA:64 )
1: MolNum(6) Atom( CA:53 ) <=> MolNum(6) Atom( CA:63 )
2: MolNum(6) Atom( O:60 ) <=> MolNum(6) Atom( O:82 )
3: MolNum(6) Atom( C:59 ) <=> MolNum(6) Atom( C:81 )
4: MolNum(6) Atom( HB2:57 ) <=> MolNum(6) Atom( HB2:67 )
5: MolNum(6) Atom( H:52 ) <=> MolNum(6) Atom( H:62 )
6: MolNum(6) Atom( N:51 ) <=> MolNum(6) Atom( N:61 )
7: MolNum(6) Atom( HB1:56 ) <=> MolNum(6) Atom( HB1:66 )
8: MolNum(6) Atom( CB:55 ) <=> MolNum(6) Atom( CB:65 )
)

.. warning::

   This used the default MCS matching algorithm built into :mod:`sire`.
   This is not supported on Windows. You may want to use a different
   matching algorithm, e.g. using
   `Kartograf <https://kartograf.readthedocs.io/en/latest/>`__

We can see that this mapping has nicely matched the atoms in alanine
to the corresponding atoms in lysine.

Next, we will align the lysine in the mapping onto the alanine.

>>> mapping = mapping.align()

And finally, we will create the merged molecule.

>>> merged = mapping.merge(as_new_molecule=False)
>>> print(merged)
Molecule( Protein:6 num_atoms=315 num_residues=19 )

This has created the merged molecule as before. To see what has changed,
we can check the perturbation...

>>> p_omm = merged.perturbation().to_openmm()
>>> print(p_omm.changed_atoms())
       atom  charge0  charge1    sigma0    sigma1  epsilon0  epsilon1  alpha0  alpha1  kappa0  kappa1
0      N:51  -0.4157  -0.3479  0.325000  0.325000  0.711280  0.711280     0.0     0.0     0.0     0.0
1      H:52   0.2719   0.2747  0.106908  0.106908  0.065689  0.065689     0.0     0.0     0.0     0.0
2     CA:53   0.0337  -0.2400  0.339967  0.339967  0.457730  0.457730     0.0     0.0     0.0     0.0
3     HA:54   0.0823   0.1426  0.247135  0.247135  0.065689  0.065689     0.0     0.0     0.0     0.0
4     CB:55  -0.1825  -0.0094  0.339967  0.339967  0.457730  0.457730     0.0     0.0     0.0     0.0
5    HB1:56   0.0603   0.0362  0.264953  0.264953  0.065689  0.065689     0.0     0.0     0.0     0.0
6    HB2:57   0.0603   0.0362  0.264953  0.264953  0.065689  0.065689     0.0     0.0     0.0     0.0
7    HB3:58   0.0603   0.0000  0.264953  0.264953  0.065689  0.000000     0.0     1.0     1.0     1.0
8      C:59   0.5973   0.7341  0.339967  0.339967  0.359824  0.359824     0.0     0.0     0.0     0.0
9      O:60  -0.5679  -0.5894  0.295992  0.295992  0.878640  0.878640     0.0     0.0     0.0     0.0
10  Xxx:303   0.0000   0.0187  0.339967  0.339967  0.000000  0.457730     1.0     0.0     1.0     1.0
11  Xxx:304   0.0000   0.0103  0.264953  0.264953  0.000000  0.065689     1.0     0.0     1.0     1.0
12  Xxx:305   0.0000   0.0103  0.264953  0.264953  0.000000  0.065689     1.0     0.0     1.0     1.0
13  Xxx:306   0.0000  -0.0479  0.339967  0.339967  0.000000  0.457730     1.0     0.0     1.0     1.0
14  Xxx:307   0.0000   0.0621  0.264953  0.264953  0.000000  0.065689     1.0     0.0     1.0     1.0
15  Xxx:308   0.0000   0.0621  0.264953  0.264953  0.000000  0.065689     1.0     0.0     1.0     1.0
16  Xxx:309   0.0000  -0.0143  0.339967  0.339967  0.000000  0.457730     1.0     0.0     1.0     1.0
17  Xxx:310   0.0000   0.1135  0.195998  0.195998  0.000000  0.065689     1.0     0.0     1.0     1.0
18  Xxx:311   0.0000   0.1135  0.195998  0.195998  0.000000  0.065689     1.0     0.0     1.0     1.0
19  Xxx:312   0.0000  -0.3854  0.325000  0.325000  0.000000  0.711280     1.0     0.0     1.0     1.0
20  Xxx:313   0.0000   0.3400  0.106908  0.106908  0.000000  0.065689     1.0     0.0     1.0     1.0
21  Xxx:314   0.0000   0.3400  0.106908  0.106908  0.000000  0.065689     1.0     0.0     1.0     1.0
22  Xxx:315   0.0000   0.3400  0.106908  0.106908  0.000000  0.065689     1.0     0.0     1.0     1.0

This shows how the atoms in alanine have been perturbed into the equivalent
atoms in lysine, plus how a number of ghost atoms have been added that
correspond to the additional atoms in lysine. In addition, the HB3 atom
of alanine is turned into a ghost.

One-stop function
-----------------

Just as before, the entire set of steps above can be performed in one
step via the :func:`sire.morph.merge` function.

>>> merged = sr.morph.merge(ala, lys)
>>> print(merged)
Molecule( Protein:3627 num_atoms=315 num_residues=19 )
>>> p_omm = merged.perturbation().to_openmm()
>>> print(p_omm.changed_atoms())
       atom  charge0  charge1    sigma0    sigma1  epsilon0  epsilon1  alpha0  alpha1  kappa0  kappa1
0      N:51  -0.4157  -0.3479  0.325000  0.325000  0.711280  0.711280     0.0     0.0     0.0     0.0
1      H:52   0.2719   0.2747  0.106908  0.106908  0.065689  0.065689     0.0     0.0     0.0     0.0
2     CA:53   0.0337  -0.2400  0.339967  0.339967  0.457730  0.457730     0.0     0.0     0.0     0.0
3     HA:54   0.0823   0.1426  0.247135  0.247135  0.065689  0.065689     0.0     0.0     0.0     0.0
4     CB:55  -0.1825  -0.0094  0.339967  0.339967  0.457730  0.457730     0.0     0.0     0.0     0.0
5    HB1:56   0.0603   0.0362  0.264953  0.264953  0.065689  0.065689     0.0     0.0     0.0     0.0
6    HB2:57   0.0603   0.0362  0.264953  0.264953  0.065689  0.065689     0.0     0.0     0.0     0.0
7    HB3:58   0.0603   0.0000  0.264953  0.264953  0.065689  0.000000     0.0     1.0     1.0     1.0
8      C:59   0.5973   0.7341  0.339967  0.339967  0.359824  0.359824     0.0     0.0     0.0     0.0
9      O:60  -0.5679  -0.5894  0.295992  0.295992  0.878640  0.878640     0.0     0.0     0.0     0.0
10  Xxx:303   0.0000   0.0187  0.339967  0.339967  0.000000  0.457730     1.0     0.0     1.0     1.0
11  Xxx:304   0.0000   0.0103  0.264953  0.264953  0.000000  0.065689     1.0     0.0     1.0     1.0
12  Xxx:305   0.0000   0.0103  0.264953  0.264953  0.000000  0.065689     1.0     0.0     1.0     1.0
13  Xxx:306   0.0000  -0.0479  0.339967  0.339967  0.000000  0.457730     1.0     0.0     1.0     1.0
14  Xxx:307   0.0000   0.0621  0.264953  0.264953  0.000000  0.065689     1.0     0.0     1.0     1.0
15  Xxx:308   0.0000   0.0621  0.264953  0.264953  0.000000  0.065689     1.0     0.0     1.0     1.0
16  Xxx:309   0.0000  -0.0143  0.339967  0.339967  0.000000  0.457730     1.0     0.0     1.0     1.0
17  Xxx:310   0.0000   0.1135  0.195998  0.195998  0.000000  0.065689     1.0     0.0     1.0     1.0
18  Xxx:311   0.0000   0.1135  0.195998  0.195998  0.000000  0.065689     1.0     0.0     1.0     1.0
19  Xxx:312   0.0000  -0.3854  0.325000  0.325000  0.000000  0.711280     1.0     0.0     1.0     1.0
20  Xxx:313   0.0000   0.3400  0.106908  0.106908  0.000000  0.065689     1.0     0.0     1.0     1.0
21  Xxx:314   0.0000   0.3400  0.106908  0.106908  0.000000  0.065689     1.0     0.0     1.0     1.0
22  Xxx:315   0.0000   0.3400  0.106908  0.106908  0.000000  0.065689     1.0     0.0     1.0     1.0

This merged molecule can be used in a free energy simulation in the same
way as any other merged molecule.

Protein Mutation
----------------

Just as for other merged molecules, we can extract the end states using
the :func:`sire.morph.extract_reference` and
:func:`sire.morph.extract_perturbed` functions.

>>> ref_prot = sr.morph.extract_reference(merged)
>>> pert_prot = sr.morph.extract_perturbed(merged)
>>> print(ref_prot.residues().names())
[ResName('LYS'), ResName('ILE'), ResName('GLY'), ResName('ALA'),
 ResName('LYS'), ResName('ILE'), ResName('LYS'), ResName('ILE'),
 ResName('GLY'), ResName('ALA'), ResName('LYS'), ResName('ILE'),
 ResName('LYS'), ResName('ILE'), ResName('GLY'), ResName('ALA'),
 ResName('LYS'), ResName('ILE'), ResName('NH2')]
>>> print(pert_prot.residues().names())
[ResName('LYS'), ResName('ILE'), ResName('GLY'), ResName('LYS'),
 ResName('LYS'), ResName('ILE'), ResName('LYS'), ResName('ILE'),
 ResName('GLY'), ResName('ALA'), ResName('LYS'), ResName('ILE'),
 ResName('LYS'), ResName('ILE'), ResName('GLY'), ResName('ALA'),
 ResName('LYS'), ResName('ILE'), ResName('NH2')]

This shows that the first alanine residue has been mutated into a
lysine residue.

Interestingly, if you are only interested in mutation, and are not
interested in the merged molecule, then mutation is just the process
of performing a merge, and then extracting the perturbed end state.

To make this easier, there is a one-stop function for this,
:func:`sire.morph.mutate`.

>>> mutated_protein = sr.morph.mutate(ala, lys)
>>> print(mutated_protein.residues().names())
[ResName('LYS'), ResName('ILE'), ResName('GLY'), ResName('LYS'),
 ResName('LYS'), ResName('ILE'), ResName('LYS'), ResName('ILE'),
 ResName('GLY'), ResName('ALA'), ResName('LYS'), ResName('ILE'),
 ResName('LYS'), ResName('ILE'), ResName('GLY'), ResName('ALA'),
 ResName('LYS'), ResName('ILE'), ResName('NH2')]

.. warning::

   By default this uses the (potentially slow) internal MCS matching
   algorithm (which also is not supported on Windows). You may want to
   use a different matching algorithm, either by passing in the
   dictionary of the mapping you want to use via the ``match``
   argument, or by using `Kartograf <https://kartograf.readthedocs.io/en/latest/>`__.
   The ``match`` argument works identically in this function as it
   does in the ``merge`` and ``match`` functions.

In this case, we copied and pasted the lysine residue from one part of the
protein over the alanine. But, you can copy and paste in this way between
different molecules. For example, you could have a template library of
different residues that you could use for mutation.

Assuming ``template["lys"]`` contained your template for a lysine residue,
then you could have run the following;

>>> mutated_protein = sr.morph.mutate(ala, template["lys"])

This works for any kind of molecules - not just proteins. In future
versions of :mod:`sire` we will add functionality to make it easier to
manage libraries of templates, and to perform molecular editing by
copying and pasting between molecules, and between fragments constructed
via, e.g. smiles strings.
