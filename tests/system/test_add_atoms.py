import pytest
import sire as sr


def test_add_atoms(tmpdir, ala_mols):
    mols = ala_mols.clone()

    mol = mols[0]

    editor = mol.edit()

    n_to_add = 10
    natoms = mol.num_atoms()

    for i in range(0, n_to_add):
        # Add the atom and put it into the first
        # residue (update as needed if you want
        # this in a different residue)
        editor = (
            editor.add(sr.mol.AtomName("Re"))
            .renumber(sr.mol.AtomNum(natoms + i + 1))
            .reparent(sr.mol.ResIdx(0))
            .molecule()
        )

    mol = editor.commit()

    # Second, we will set the properties of the
    # atoms that we have added. We will do this using
    # the new cursor syntax
    cursor = mol.cursor()["atomname Re"]

    # This cursor can be used to edit all of the "Re" atoms
    # at once, or can be indexed to edit them individually

    # We will iterate over them and set some default properties
    for atom in cursor.atoms():
        atom["coordinates"] = sr.maths.Vector(0)
        atom["charge"] = 0 * sr.units.mod_electron
        atom["element"] = sr.mol.Element(0)
        atom["mass"] = 0 * sr.units.g_per_mol
        atom["atomtype"] = "DM"

    # Adding atoms has invalidated the "parameters" property, so
    # this should be removed
    cursor = cursor.molecule()

    if "parameters" in cursor:
        del cursor["parameters"]

    mol = cursor.commit()

    # You can also add bonds to these atoms, e.g.
    connectivity = mol.property("connectivity").edit()

    # here we connect the first Re atom to the first atom in the molecule...
    connectivity.connect(mol["atomname Re"][0].index(), mol.atoms()[0].index())

    # also should set a bond potential (even if this is zero)
    bonds = mol.property("bond")
    bonds.set(mol["atomname Re"][0].index(), mol.atoms()[0].index(), 0.0)

    mol = (
        mol.edit()
        .set_property("connectivity", connectivity.commit())
        .set_property("bond", bonds)
        .commit()
    )

    # Finally(!) we update the system with the new version of this molecule...
    mols.update(mol)

    # ...then we write this out to a new PRMTOP/RST file
    d = tmpdir.mkdir("test_add_atoms")
    f = sr.save(mols, d.join("test"), format=["PRMTOP", "RST"])

    # Reload the file to check everything is ok
    mols = sr.load(f)

    mol = mols[0]

    assert mol.num_atoms() == ala_mols[0].num_atoms() + n_to_add
    assert len(mol.bonds("atomname Re")) == 1


def test_add_dummies(tmpdir, ala_mols):
    mols = ala_mols.clone()

    n_existing_atoms = mols.num_atoms()
    n_existing_residues = mols.num_residues()

    newmol = sr.mol.Molecule("dummies")

    n_to_add = 10

    editor = newmol.edit()

    # Create a residue
    editor = (
        editor.add(sr.mol.ResName("Re"))
        .renumber(sr.mol.ResNum(n_existing_residues + 1))
        .molecule()
    )

    for i in range(0, n_to_add):
        editor = (
            editor.add(sr.mol.AtomName("Re"))
            .renumber(sr.mol.AtomNum(n_existing_residues + i + 1))
            .reparent(sr.mol.ResIdx(0))
            .molecule()
        )

    mol = editor.commit()

    cursor = mol.cursor()["atomname Re"]

    # need to set the properties to the correct type...
    cursor[0]["charge"] = 1 * sr.units.mod_electron
    cursor[0]["mass"] = 1 * sr.units.g_per_mol

    for atom in cursor.atoms():
        atom["coordinates"] = sr.maths.Vector(0)
        atom["charge"] = 0 * sr.units.mod_electron
        atom["element"] = sr.mol.Element(0)
        atom["mass"] = 0 * sr.units.g_per_mol
        atom["atomtype"] = "DM"
        atom["LJ"] = sr.mm.LJParameter(
            1 * sr.units.angstrom, 0 * sr.units.kcal_per_mol
        )

    mol = cursor.molecule().commit()

    # This line is a bit janky as there isn't yet a "modern API"
    # way to do this
    mols = mols.molecules()
    mols.append(mol)

    d = tmpdir.mkdir("test_add_dummies")
    f = sr.save(mols, d.join("test"), format=["PRM7", "RST7"])

    # load to check
    mols = sr.load(f)

    mol = mols[-1]

    assert mol.num_atoms() == 10

    for atom in mol.atoms():
        assert atom.name().value() == "Re"
