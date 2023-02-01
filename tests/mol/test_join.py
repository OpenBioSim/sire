import sire as sr


def test_join_single(ala_mols):
    # in the old API, len(mol) would return 1
    # This tripped up sire and caused errors in
    # waterswap 2023.0.0 to 2023.1.1
    orig__len__ = sr.legacy.Mol.Molecule.__len__

    try:
        mols = ala_mols
        mol = mols[0]
        molnum = mol.number()

        sr.legacy.Mol.Molecule.__len__ = lambda obj: 1

        group = sr.legacy.Mol.MoleculeGroup("test")
        group.add(mol)

        m = group[molnum]

        ff1 = sr.legacy.MM.InternalFF("ff1")
        ff1.add(sr.legacy.Mol.Molecule.join(m))

        ff2 = sr.legacy.MM.InternalFF("ff2")
        ff2.add(mol)

        assert ff1.energy() == ff2.energy()

        ff1 = sr.legacy.MM.InternalFF("ff1")
        ff1.add(sr.legacy.Mol.Molecule.join(mol))

        ff2 = sr.legacy.MM.InternalFF("ff2")
        ff2.add(mol)

        assert ff1.energy() == ff2.energy()

        ff1 = sr.legacy.MM.InternalFF("ff1")
        ff1.add(sr.legacy.Mol.Molecule.join(mol))
        assert ff1.energy() == ff2.energy()
    except Exception:
        sr.legacy.Mol.Molecule.__len__ = orig__len__
        raise

    sr.legacy.Mol.Molecule.__len__ = orig__len__


def test_join_all(ala_mols):
    # in the old API, len(mol) would return 1
    # This tripped up sire and caused errors in
    # waterswap 2023.0.0 to 2023.1.1
    orig__len__ = sr.legacy.Mol.Molecule.__len__

    try:
        mols = ala_mols
        mol = mols[0]
        molnum = mol.number()

        sr.legacy.Mol.Molecule.__len__ = lambda obj: 1

        group1 = sr.legacy.Mol.MoleculeGroup("test1")
        group2 = sr.legacy.Mol.MoleculeGroup("test2")

        group1.add(mol)
        group2.add(mol)

        group = sr.legacy.Mol.MoleculeGroups()
        group.add(group1)
        group.add(group2)

        m = group[molnum]

        ff1 = sr.legacy.MM.InternalFF("ff1")
        ff1.add(sr.legacy.Mol.Molecule.join(m))

        ff2 = sr.legacy.MM.InternalFF("ff2")
        ff2.add(mol)

        assert ff1.energy() == ff2.energy()

    except Exception:
        sr.legacy.Mol.Molecule.__len__ = orig__len__
        raise

    sr.legacy.Mol.Molecule.__len__ = orig__len__


def test_join_subset(ala_mols):
    # in the old API, len(mol) would return 1
    # This tripped up sire and caused errors in
    # waterswap 2023.0.0 to 2023.1.1
    orig__len__ = sr.legacy.Mol.Molecule.__len__

    try:
        mols = ala_mols
        mol = mols[0]
        molnum = mol.number()

        sr.legacy.Mol.Molecule.__len__ = lambda obj: 1

        for subset in [mol, mol["residx 0"], mol["atomidx < 5"]]:
            # add the atoms, one by one to the group
            group1 = sr.legacy.Mol.MoleculeGroup("test1")
            group2 = sr.legacy.Mol.MoleculeGroup("test2")

            for i, atom in enumerate(subset.atoms()):
                if i % 2 == 0:
                    group1.add(atom)
                else:
                    group2.add(atom)

            group = sr.legacy.Mol.MoleculeGroups()
            group.add(group1)
            group.add(group2)

            m = group[molnum]

            ff1 = sr.legacy.MM.InternalFF("ff1")
            ff1.add(sr.legacy.Mol.Molecule.join(m))

            ff2 = sr.legacy.MM.InternalFF("ff2")
            ff2.add(subset)

            assert ff1.energy() == ff2.energy()
    except Exception:
        sr.legacy.Mol.Molecule.__len__ = orig__len__
        raise

    sr.legacy.Mol.Molecule.__len__ = orig__len__
