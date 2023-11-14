import sire as sr

import pytest


def _assert_equal(v0, v1, tol):
    assert v0.x() == pytest.approx(v1.x(), tol)
    assert v0.y() == pytest.approx(v1.y(), tol)
    assert v0.z() == pytest.approx(v1.z(), tol)


@pytest.mark.skipif(
    "gemmi" not in sr.convert.supported_formats(), reason="gemmi not available"
)
def test_gemmi(testfile_cache_dir, pdbx_3nss):
    mols = pdbx_3nss

    import os

    cif_file = os.path.join(testfile_cache_dir, "3NSS.cif")

    import gemmi

    s = gemmi.read_structure(cif_file)

    mols2 = sr.convert.to(s, "sire")

    assert mols.num_molecules() == mols2.num_molecules()
    assert mols.num_atoms() == mols2.num_atoms()

    for atom0, atom1 in zip(mols.atoms(), mols2.atoms()):
        assert atom0.name() == atom1.name()
        assert atom0.number() == atom1.number()
        assert atom0.element() == atom1.element()
        assert atom0.residue().name() == atom1.residue().name()
        assert atom0.residue().number() == atom1.residue().number()
        _assert_equal(atom0.coordinates(), atom1.coordinates(), 1e-6)

    s2 = sr.convert.to(mols, "gemmi")

    assert len(s) == len(s2)
    assert len(list(s[0].all())) == len(list(s2[0].all()))

    r0 = s[0].get_all_residue_names()
    r1 = s2[0].get_all_residue_names()

    assert len(r0) == len(r1)

    for r in r0:
        assert r in r1
