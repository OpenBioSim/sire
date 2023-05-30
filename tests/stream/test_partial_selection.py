import sire as sr

import pytest


def test_partial_selection(tmpdir, ala_mols):
    mols = ala_mols

    mol = mols[0]

    res = mol.residues()[0]

    assert res.selection().num_selected_residues() == 1

    dir = tmpdir.mkdir("test_partial_selection")

    s3file = str(dir.join("output.s3"))

    sr.stream.save(res, s3file)

    res2 = sr.stream.load(s3file)

    assert res2.selection().num_selected_residues() == 1
    assert res2.residue().number() == res.residue().number()

    res = sr.legacy.Mol.PartialMolecule(res)

    assert res.selection().num_selected_residues() == 1

    s3file = str(dir.join("output2.s3"))

    sr.stream.save(res, s3file)

    res2 = sr.stream.load(s3file)

    assert res2.selection().num_selected_residues() == 1
    assert res2.residue().number() == res.residue().number()
