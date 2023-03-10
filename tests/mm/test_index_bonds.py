
import pytest


def _assert_same_bonds(b0, b1):
    assert len(b0) == len(b1)

    for bond0 in b0:
        same = False

        for bond1 in b1:
            if bond0 == bond1:
                same = True
                break

        if not same:
            print(f"Not the same!\n{b0}\n{b1}")
            assert False


def test_index_bonds(ala_mols):
    mols = ala_mols

    mol = mols[0]

    # assert we only find the 5 bonds wholly within the residue
    assert len(mol["bonds in resnum 1"]) == 5

    # assert we also find the extra bond that is to resnum 2
    # (i.e. all bonds involving any atom in the residue)
    assert len(mol["bonds with (atoms in resnum 1)"]) == 6

    # check that we don't need brackets
    assert len(mol["bonds with atoms in resnum 1"]) == 6

    _assert_same_bonds(mol["bonds in resnum 1"], mol["resnum 1"].bonds())
    _assert_same_bonds(mol["bonds with (atoms in resnum 1)"], mol.bonds("resnum 1"))
    _assert_same_bonds(mol["bonds with atoms in resnum 1"], mol.bonds("resnum 1"))

    # this is a single bond
    assert mol["bonds to resnum 1"] == mol.bond("resnum 1", "resnum 2")

    # still just a single bond
    assert mol["bonds from resnum 1 to resnum 2"] == mol.bond("resnum 1", "resnum 2")

    cx = mol["bonds to element C"]
    ccx = mol["bonds with element C"]
    cc = mol["bonds from element C to element C"]

    assert len(cc) + len(cx) == len(ccx)

    assert len(cc) == 3
    assert len(cx) == 16
    assert len(ccx) == 19

    _assert_same_bonds(ccx["bonds from element C to element C"], cc)


def test_index_mols_bonds(ala_mols):
    mols = ala_mols

    import sire as sr

    bnds = mols["bonds from element O to element H"]

    assert len(bnds) == 2 * (mols.num_molecules() - 1)
    assert len(bnds) == 1260  # just in case something else failed!

    assert bnds.num_molecules() == mols.num_molecules() - 1
    assert bnds.num_atoms() == bnds.num_molecules() * 3

    assert len(mols[sr.bondid(0, 1)]) == mols.num_molecules()

    for bond in mols[sr.bondid(0, 1)]:
        assert bond.atom0().index().value() == 0
        assert bond.atom1().index().value() == 1

    assert len(mols[sr.bondid("O", "H1")]) == mols.num_molecules() - 1

    for bond in mols[sr.bondid("O", "H1")]:
        assert bond.atom0().name().value() in ["O", "H1"]
        assert bond.atom1().name().value() in ["O", "H1"]

    # mols[1:] are the water molecules
    assert mols["bonds from element O to element H"].mass().value() == \
	pytest.approx((mols[1:]["element O"].mass() + mols[1:]["element H"].mass()).value(), 0.0001)

if __name__ == "__main__":
    test_index_bonds()
    test_index_mols_bonds()


