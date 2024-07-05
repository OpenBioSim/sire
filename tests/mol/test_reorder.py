def test_reorder_atoms(ala_mols):
    mols = ala_mols

    mol = mols[0]

    # check that reorder preserves residue cutting
    mol = mol.edit()
    atom = mol.atom(2)
    atom = atom.reindex(0)
    mol = atom.molecule().commit()

    assert mol.num_residues() == 3
    assert mol.num_cutgroups() == mol.num_residues()

    atomidx = 0

    for cutgroup in mol.cutgroups():
        for atom in cutgroup.atoms():
            assert atom.index().value() == atomidx
            atomidx += 1

    # now check with reordering that break residue cutting
    mol = mols[0]

    mol = mol.edit()
    atom = mol.atom("HA")
    atom = atom.reindex(0)
    mol = atom.molecule().commit()

    assert mol.num_cutgroups() == 1

    atomidx = 0

    for atom in mol.cutgroups()[0].atoms():
        assert atom.index().value() == atomidx
        atomidx += 1
