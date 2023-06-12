import pytest
import sire as sr


@pytest.mark.veryslow
def test_cljpairs(neura_mols):
    mols = neura_mols

    for mol in mols[0:5]:
        scl = mol.property("intrascale")

        # Make sure that the excluded atoms gained
        # by acccessing pairs individually is identical
        # to that gained by accessing them by CutGroup
        for cg in mol.cutgroups():
            p = scl.excluded_atoms(cg.index())

            for atom in p.keys():
                e1 = p[atom]
                e2 = scl.excluded_atoms(atom)

                assert e1 == e2
