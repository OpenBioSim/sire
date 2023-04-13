import sire as sr
import pytest


def test_ambtraj(tmpdir, kigaki_mols):
    mols = kigaki_mols

    d = tmpdir.mkdir("test_ambtraj")

    f = sr.save(mols, d.join("test"), format=["TRAJ", "PRM7"])

    newmols = sr.load(f, show_warnings=False)

    def compare(v1, v2):
        assert v1.x().value() == pytest.approx(v2.x().value(), 0.01)
        assert v1.y().value() == pytest.approx(v2.y().value(), 0.01)
        assert v1.z().value() == pytest.approx(v2.z().value(), 0.01)

    for mol, newmol in zip(mols[0:20], newmols[0:20]):
        for atom, newatom in zip(mol.atoms(), newmol.atoms()):
            compare(
                atom.property("coordinates"), newatom.property("coordinates")
            )
