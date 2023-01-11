import pytest


def test_find(ala_mols):

    mols = ala_mols

    assert mols.atoms().find(mols.atoms()[0]) == 0

    assert mols.atoms().find(mols.atoms()[1000]) == 1000

    assert mols.atoms().find(mols.atoms()[10:100]) == list(range(10, 100))

    idxs = mols.atoms().find(mols["element C"])

    assert len(idxs) == len(mols["element C"])

    for idx in idxs:
        assert mols.atoms()[idx].element().num_protons() == 6

    assert mols.residues().find(mols["resname ALA"]) == 1

    assert mols.find(mols[600]) == 600

    assert mols.find(mols[3:50:3]) == list(range(3, 50, 3))

    with pytest.raises(IndexError):
        mols[1:].find(mols[0])

    with pytest.raises(IndexError):
        mols[1:].atoms().find(mols["element C"])
