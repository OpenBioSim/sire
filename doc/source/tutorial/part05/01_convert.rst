====================
Converting Molecules
====================

There are many available excellent molecular modelling packages for Python.
Rather than re-implementing their functionality, a design goal of
:mod:`sire` is to support easy interconversion so that functionality
from those packages can be built on and re-used.

The :mod:`sire.convert` module contains functions that support the
conversion of the :mod:`sire` molecular objects into their equivalents
in a number of other packages.

BioSimSpace
-----------

`BioSimSpace <https://biosimspace.openbiosim.org>`__
is an interoperable Python framework for biomolecular simulation that
makes it easy for you to write simulation workflow components for, e.g.
parameterising molecules, solvating systems, or running molecular
dynamics or free energy simulations using a number of different
backends (e.g. Amber, Gromacs, Namd etc.).

To convert a :mod:`sire` object to a BioSimSpace object, you need to
ensure that you have imported BioSimSpace first in your script, before
you import :mod:`sire`.

>>> import BioSimSpace as BSS
>>> import sire as sr

To convert to BioSimSpace, you just need to pass the :mod:`sire` molecular
object to the function :func:`sire.convert.to`, e.g.

>>> mols = sr.load(sr.expand(sr.tutorial_url, "ala.crd", "ala.top"))
>>> mol = mols[0]
>>> bss_mol = sr.convert.to(mol, "BioSimSpace")
>>> print(bss_mol)
<BioSimSpace.Molecule: number=2, nAtoms=22, nResidues=3>

You can now use the BioSimSpace Molecule exactly as you would any other
BioSimSpace Molecule, e.g. calling BioSimSpace functions to generate
forcefield parameters, solvate it, or run molecular dynamics or
free energy simulations.

.. note::

   Note that the format argument to ``sr.convert.to`` is case-insensitive.
   You could have use ``bss_mol = sr.convert.to(mol, "biosimspace")``.

:func:`sire.convert.to` will convert any sub-view of a
molecule (e.g. a :class:`~sire.mol.Residue` or :class:`~sire.mol.Atom`)
into the BioSimSpace Molecule that contains that view.

>>> atom = mol[0]
>>> bss_mol = sr.convert.to(atom, "BioSimSpace")
>>> print(bss_mol)
<BioSimSpace.Molecule: number=2, nAtoms=22, nResidues=3>

It will convert a list of :class:`~sire.mol.Molecule` or sub-view
objects into a list of equivalent BioSimSpace Molecules, e.g.

>>> bss_mols = sr.convert.to(mols[0:10], "BioSimSpace")
>>> print(bss_mols)
[<BioSimSpace.Molecule: number=2, nAtoms=22, nResidues=3>,
 <BioSimSpace.Molecule: number=3, nAtoms=3, nResidues=1>,
 <BioSimSpace.Molecule: number=4, nAtoms=3, nResidues=1>,
 <BioSimSpace.Molecule: number=5, nAtoms=3, nResidues=1>,
 <BioSimSpace.Molecule: number=6, nAtoms=3, nResidues=1>,
 <BioSimSpace.Molecule: number=7, nAtoms=3, nResidues=1>,
 <BioSimSpace.Molecule: number=8, nAtoms=3, nResidues=1>,
 <BioSimSpace.Molecule: number=9, nAtoms=3, nResidues=1>,
 <BioSimSpace.Molecule: number=10, nAtoms=3, nResidues=1>,
 <BioSimSpace.Molecule: number=11, nAtoms=3, nResidues=1>]

 There is also logic to convert a :class:`~sire.system.System` object,
 which is the collection of molecules plus associated metadata read
 from an input file, into the equivalent BioSimSpace System class.

 >>> bss_sys = sr.convert.to(mols, "BioSimSpace")
 >>> print(bss_sys)
 <BioSimSpace.System: nMolecules=631>

You can convert back from BioSimSpace to :mod:`sire` using the same
function. For example,

>>> mol = sr.convert.to(bss_mol, "sire")
>>> print(mol)
Molecule( ACE:2   num_atoms=22 num_residues=3 )

Note that any sub-view of a BioSimSpace object will be converted to the
:class:`~sire.mol.Molecule` that contains that view, e.g.

>>> mol = sr.convert.to(bss_mol.getAtoms()[0], "sire")
>>> print(mol)
Molecule( ACE:2   num_atoms=22 num_residues=3 )

A list of BioSimSpace molecule (or sub-view) objects will be converted to
a list of :class:`~sire.mol.Molecule` objects.

>>> mols = sr.convert.to(bss_mols, "sire")
>>> print(mols)
SelectorMol( size=10
0: Molecule( ACE:2   num_atoms=22 num_residues=3 )
1: Molecule( WAT:3   num_atoms=3 num_residues=1 )
2: Molecule( WAT:4   num_atoms=3 num_residues=1 )
3: Molecule( WAT:5   num_atoms=3 num_residues=1 )
4: Molecule( WAT:6   num_atoms=3 num_residues=1 )
5: Molecule( WAT:7   num_atoms=3 num_residues=1 )
6: Molecule( WAT:8   num_atoms=3 num_residues=1 )
7: Molecule( WAT:9   num_atoms=3 num_residues=1 )
8: Molecule( WAT:10  num_atoms=3 num_residues=1 )
9: Molecule( WAT:11  num_atoms=3 num_residues=1 )
)

And a BioSimSpace System will be automatically converted to a
:class:`~sire.system.System` object.

>>> mols = sr.convert.to(bss_sys, "sire")
>>> print(mols)
System( name=ACE num_molecules=631 num_residues=633 num_atoms=1912 )

RDKit
-----

`RDKit <https://rdkit.org>`__ is a collection of cheminformatics and
machine-learning software written in C++ and Python. Assuming you have
RDKit installed, you can convert :mod:`sire` molecule and molecule view
objects to and from RDKit Molecule objects.

The :func:`sire.convert.supported_formats` function lists the formats that
:mod:`sire.convert` supports for the current installation. This will
depend on whether or not you have the package installed in the same conda
environment as :mod:`sire`, and whether or not :mod:`sire` was compiled
with support for that package.

>>> print(sr.convert.supported_formats())
['biosimspace', 'gemmi', 'openmm', 'rdkit', 'sire']

.. note::

   If ``rdkit`` isn't listed, then you should quit Python and install
   it, e.g. using the command ``mamba install -c conda-forge rdkit``.
   If it still isn't listed then please raise an issue on the
   `sire GitHub repository <https://github.com/OpenBioSim/sire/issues>`__.

You can convert to RDKit by passing ``rdkit`` as the format argument to
:func:`sire.convert.to`, e.g.

>>> rdkit_mol = sr.convert.to(mol, "rdkit")
>>> print(rdkit_mol)
<rdkit.Chem.rdchem.Mol object at 0x10283da10>

You can now use this RDKit Mol object identically to any other
RDKit Mol object, e.g. generating smiles strings, performing
sub-structure searches, maximum common substructure alignments,
generating 2D views etc.

Just as for BioSimSpace, :func:`sire.convert.to` will return the RDKit Mol
for the entire molecule that contains any sub-views that are passed.
For example,

>>> rdkit_mol = sr.convert.to(mol[0], "rdkit")
>>> print(rdkit_mol.GetNumAtoms())
22

Passing in a list of molecules or molecule views to convert will return
a list of RDKit Mol objects.

>>> rdkit_mols = sr.convert.to(mols[0:10], "rdkit")
>>> print(rdkit_mols)
[<rdkit.Chem.rdchem.Mol object at 0x102c6a180>, <rdkit.Chem.rdchem.Mol object at 0x102c6a340>,
<rdkit.Chem.rdchem.Mol object at 0x102c6a3b0>, <rdkit.Chem.rdchem.Mol object at 0x102c69d90>,
<rdkit.Chem.rdchem.Mol object at 0x102c6a1f0>, <rdkit.Chem.rdchem.Mol object at 0x102c69af0>,
<rdkit.Chem.rdchem.Mol object at 0x102c698c0>, <rdkit.Chem.rdchem.Mol object at 0x102c69bd0>,
<rdkit.Chem.rdchem.Mol object at 0x102c69e00>, <rdkit.Chem.rdchem.Mol object at 0x102c69cb0>]

.. note::

   RDKit does not have an equivalent of a :class:`~sire.system.System` object,
   so these will be converted to a list of RDKit Mol objects.

You can also convert RDKit Mol objects back to :class:`~sire.mol.Molecule`
objects, e.g.

>>> mol = sr.convert.to(rdkit_mol, "sire")
>>> print(mol)
Molecule( ACE:633 num_atoms=22 num_residues=1 )

>>> mols = sr.convert.to(rdkit_mols, "sire")
>>> print(mols)
SelectorMol( size=10
0: Molecule( ACE:634 num_atoms=22 num_residues=1 )
1: Molecule( WAT:635 num_atoms=3 num_residues=1 )
2: Molecule( WAT:636 num_atoms=3 num_residues=1 )
3: Molecule( WAT:637 num_atoms=3 num_residues=1 )
4: Molecule( WAT:638 num_atoms=3 num_residues=1 )
5: Molecule( WAT:639 num_atoms=3 num_residues=1 )
6: Molecule( WAT:640 num_atoms=3 num_residues=1 )
7: Molecule( WAT:641 num_atoms=3 num_residues=1 )
8: Molecule( WAT:642 num_atoms=3 num_residues=1 )
9: Molecule( WAT:643 num_atoms=3 num_residues=1 )
)

This is useful, e.g. if you have created the molecule using RDKit's
smiles functionality, and then want to convert to a :class:`~sire.mol.Molecule`
object for continued manipulation.

OpenMM
------

`OpenMM <https://openmm.org>`__ is a high-performance toolkit for molecular
simulation, which is particularly suited to running GPU-accelerated
molecular dynamics (and related) simulations.

The :func:`sire.convert.supported_formats` function lists the formats that
:mod:`sire.convert` supports for the current installation. This will
depend on whether or not you have the package installed in the same conda
environment as :mod:`sire`, and whether or not :mod:`sire` was compiled
with support for that package.

>>> print(sr.convert.supported_formats())
['biosimspace', 'gemmi', 'openmm', 'rdkit', 'sire']

.. note::

   If ``openmm`` isn't listed, then you should quit Python and install
   it, e.g. using the command ``mamba install -c conda-forge openmm``.
   If it still isn't listed then please raise an issue on the
   `sire GitHub repository <https://github.com/OpenBioSim/sire/issues>`__.

:func:`sire.convert.to` can convert a :class:`~sire.mol.Molecule` or
molecule view into the equivalent for OpenMM.

>>> mols = sr.load(sr.expand(sr.tutorial_url, "ala.crd", "ala.top"))
>>> omm = sr.convert.to(mols[0], "openmm")
>>> print(omm)
<openmm.openmm.Context; proxy of <Swig Object of type 'OpenMM::Context *' at 0x14e95b510> >

The result is an OpenMM Context object which contains just the first molecule
from ``mols``. This can be used just like any
other OpenMM Context object, e.g. for running minimisation or dynamics.

An OpenMM Context object is returned because it contains within it:

* Representations of the potentials and connectivity of the molecule(s)
  in the OpenMM System object (obtained via ``omm.getSystem()``)
* The current coordinates and (optionally) the velocities of the
  molecule(s) in the OpenMM Integrator object (obtained via
  ``omm.getIntegrator()``).

The context has placed these two object onto an OpenMM Platform object
(obtained via ``omm.getPlatform()``), so that the Context is ready for
simulation. You can change the platform or choose a new Integrator by
using the ``omm.getSystem()`` or ``omm.getIntegrator()`` to extract these
objects and then recombine them with the platform or integrator of your choice.

More detail on how you can control what platform and integrator is chosen
for this conversion is :doc:`available here <../../cheatsheet/openmm>`.

You can convert a single molecule, list of molecules or an entire
:class:`~sire.system.System` to an OpenMM context in the same way, e.g.

>>> omm = sr.convert.to(mols, "openmm")
>>> print(omm)
<openmm.openmm.Context; proxy of <Swig Object of type 'OpenMM::Context *' at 0x14e9ee220> >

We do plan to add code to allow conversion back from an OpenMM Context to
the equivalent :mod:`sire` object, but this is not yet ready for release.

Instead, we have lower-level functions that extract coordinates, velocities
and :class:`~sire.vol.Space` objects from OpenMM State objects that are
extracted from the Context. Please do get in touch with us if you would like
to learn about these functions, and would like to contribute to coding
a more complete OpenMM to sire converter.

Gemmi
-----

`Gemmi <https://gemmi.readthedocs.io/en/latest/>`__ is a Python library
developed primarily for use in macromolecular crystallography (MX).
In particular it can be used to parse PDBx/mmCIF files, refinement
restraints, reflection data, 3D grid data and dealing with
crystallographic symmetry. This is useful for structural bioinformatics.

The :func:`sire.convert.supported_formats` function lists the formats that
:mod:`sire.convert` supports for the current installation. This will
depend on whether or not you have the package installed in the same conda
environment as :mod:`sire`, and whether or not :mod:`sire` was compiled
with support for that package.

>>> print(sr.convert.supported_formats())
['biosimspace', 'gemmi', 'openmm', 'rdkit', 'sire']

.. note::

   If ``gemmi`` isn't listed, then you should quit Python and install
   it, e.g. using the command ``conda install -c conda-forge gemmi``.
   If it still isn't listed then please raise an issue on the
   `sire GitHub repository <https://github.com/OpenBioSim/sire/issues>`__.

:func:`sire.convert.to` can convert a :class:`~sire.system.System`, list
of molecules, or single molecule into a
`Gemmi Structure <https://gemmi.readthedocs.io/en/latest/mol.html#structure>`__
object.

>>> mols = sr.load(sr.expand(sr.tutorial_url, "ala.crd", "ala.top"))
>>> gemmi_struct = sr.convert.to(mols, "gemmi")
>>> print(gemmi_struct)
<gemmi.Structure  with 1 model(s)>
>>> print(gemmi_struct[0].get_all_residue_names())
['ACE', 'ALA', 'NME', 'WAT']

Passing in a single molecule or subset of molecules will return a
`Gemmi Structure <https://gemmi.readthedocs.io/en/latest/mol.html#structure>`__
with just those molecules, e.g.

>>> gemmi_struct = sr.convert.to(mols[0], "gemmi")
>>> print(gemmi_struct)
<gemmi.Structure  with 1 model(s)>
>>> print(gemmi_struct[0].get_all_residue_names())
['ACE', 'ALA', 'NME']

You can convert a
`Gemmi Structure <https://gemmi.readthedocs.io/en/latest/mol.html#structure>`__
back to a :class:`~sire.system.System` object, e.g.

>>> mols = sr.convert.to(gemmi_struct, "sire")
>>> print(mols)
System( name= num_molecules=1 num_residues=3 num_atoms=22 )

Anything to Anything
--------------------

Above you have seen how :func:`sire.convert.to` can convert to and from
:mod:`sire` objects and other molecular modelling package objects.

It is actually more powerful than that! It recognises the object being
passed and can convert between any two object types that are supported
by :mod:`sire`, using the :mod:`sire` object format as an
intermediary.

For example, you can convert RDKit objects to BioSimSpace objects...

>>> import BioSimSpace as BSS
>>> import sire as sr
>>> from rdkit import Chem
>>> rdkit_mol = Chem.MolFromSmiles("Cc1ccccc1")
>>> bss_mol = sr.convert.to(rdkit_mol, "BioSimSpace")
>>> bss_mol
<BioSimSpace.Molecule: number=2, nAtoms=7, nResidues=1>

.. note::

   Remember that you may need to exit Python and then restart to ensure
   that BioSimSpace is imported before sire.

Or you could convert BioSimSpace molecules back to RDKit...

>>> rdkit_mol = sr.convert.to(bss_mol, "rdkit")
>>> print(Chem.MolToSmiles(rdkit_mol))
[C-3]c1[c][c][c][c][c]1

or you could setup and parameterise a molecule in BioSimSpace and
convert it to an OpenMM Context ready for minimisation or dynamics...

>>> url = BSS.tutorialUrl()
>>> bss_system = BSS.IO.readMolecules([f"{url}/ala.top", f"{url}/ala.crd"])
>>> omm = sr.convert.to(bss_system, "openmm")
>>> integrator = omm.getIntegrator()
>>> integrator.step(10)
>>> print(omm.getState().getTime())
0.010000000000000002 ps

or you could load a PDBx file from Gemmi and convert a "MAN" moelcule within it
into an RDKit structure...

>>> import gemmi
>>> import sire as sr
>>> import rdkit
>>> from rdkit import Chem
>>> import urllib
>>> urllib.request.urlretrieve("https://files.rcsb.org/download/3NSS.cif.gz",
...                            filename="3NSS.cif.gz")
>>> struct = gemmi.read_structure("3NSS.cif.gz")
>>> mol = gemmi.Selection("(MAN)").copy_structure_selection(struct)
>>> rdkit_mol = sr.convert.to(mol, "rdkit")
>>> print(Chem.MolToSmiles(rdkit_mol))
O=C=C1OC(OC2=C(=O)C(=C=O)O[C-]=C2[O-])=C(=O)=C(=O)C1=O

Supporting other formats
------------------------

We are actively looking for other molecular modelling packages to support.
Please get in touch if you would like to suggest a package we should
look at, or if you want to provide some help with implementation.
