import pytest

import sire as sr


def test_index_angles(ala_mols):
    mols = ala_mols[1:10]

    angs0 = mols.angles("element H", "element O", "element H")

    assert len(angs0) == len(mols)

    angs1 = mols.angles(mols["element H"], mols["element O"], mols["element H"])

    assert len(angs0) == len(angs1)

    for a0, a1 in zip(angs0, angs1):
        assert a0 == a1

    angs0 = mols.angles(mols.atoms()[0], mols.atoms()[1], mols.atoms()[2])
    angs1 = mols[0].angles("atomidx 0", "atomidx 1", "atomidx 2")

    assert len(angs0) == len(angs1) == 1

    assert angs0[0] == angs1[0]


def test_index_dihedrals(ala_mols):
    mol = ala_mols[0]

    dihs0 = mol.dihedrals("element H", "element C", "element N", "element H")

    assert len(dihs0) == 4

    for dih in dihs0:
        assert len(dih["element H"]) == 2
        assert len(dih.atoms("element C")) == 1
        assert len(dih.atoms("element N")) == 1

    dihs1 = mol.dihedrals(
        mol["element H"], mol["element C"], mol["element N"], mol["element H"]
    )

    assert len(dihs0) == len(dihs1)

    for d0, d1 in zip(dihs0, dihs1):
        assert d0 == d1


def test_index_impropers(ala_mols):
    mol = ala_mols[0]

    imps0 = mol.impropers("element C", "element C", "element N", "element O")

    assert len(imps0) == 2

    for imp in imps0:
        assert len(imp["element C"]) == 2
        assert len(imp.atoms("element N")) == 1
        assert len(imp.atoms("element O")) == 1

    imps1 = mol.impropers(
        mol["element C"], mol["element C"], mol["element N"], mol["element O"]
    )

    assert len(imps0) == len(imps1)

    for i0, i1 in zip(imps0, imps1):
        assert i0 == i1
