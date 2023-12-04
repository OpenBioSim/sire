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
    assert mols.num_residues() == mols2.num_residues()
    assert mols.num_chains() == mols2.num_chains()

    s2 = sr.convert.to(mols, "gemmi")

    assert len(s) == len(s2)
    assert len(list(s[0].all())) == len(list(s2[0].all()))

    r0 = s[0].get_all_residue_names()
    r1 = s2[0].get_all_residue_names()

    assert len(r0) == len(r1)

    for r in r0:
        assert r in r1


@pytest.mark.skipif(
    "gemmi" not in sr.convert.supported_formats(), reason="gemmi not available"
)
def test_gemmi_roundtrip(tmpdir, pdbx_3nss):
    mols = pdbx_3nss.clone()

    mols.set_metadata("cat", "meow")
    mols.set_metadata("dog", ["growl", "bark", "woof"])
    mols.set_metadata(
        "vehicle",
        {
            "car": ["mercedes", "ferrari", "jaguar"],
            "bike": ["yamaha", "harley", "honda"],
            "numbers": [2, 4, 6],
        },
    )

    g = sr.convert.to(mols, "gemmi")

    mols2 = sr.convert.to(g, "sire")

    assert mols.num_molecules() == mols2.num_molecules()
    assert mols.num_atoms() == mols2.num_atoms()
    assert mols.num_residues() == mols2.num_residues()

    assert mols.metadata() == mols2.metadata()

    d = tmpdir.mkdir("test_gemmi_roundtrip")

    f = sr.save(mols, d.join("test.pdbx"))

    if isinstance(f, list):
        f = f[0]

    mols3 = sr.load(f)

    assert mols.num_molecules() == mols3.num_molecules()
    assert mols.num_atoms() == mols3.num_atoms()
    assert mols.num_residues() == mols3.num_residues()

    m = mols3.metadata()

    assert m["cat"] == "meow"
    assert m["dog"] == ["growl", "bark", "woof"]

    v = m["vehicle"]

    assert v["car"].as_array() == ["mercedes", "ferrari", "jaguar"]
    assert v["bike"].as_array() == ["yamaha", "harley", "honda"]
    assert v["numbers"].as_array() == ["2", "4", "6"]


@pytest.mark.skipif(
    "gemmi" not in sr.convert.supported_formats(), reason="gemmi not available"
)
def test_gemmi_complex_metadata(tmpdir, ala_mols):
    mols = ala_mols.clone()

    mols.set_metadata("name", "alanine dipeptide")
    mols.set_metadata("residues", ["ACE", "ALA", "NME"])
    mols.set_metadata(
        "atoms",
        {
            "element": ["C", "N", "O"],
            "x_coords": [0.0, 1.0, 2.0],
            "y_coords": [3.0, 4.0, 5.0],
            "z_coords": [6.0, 7.0, 8.0],
        },
    )

    d = tmpdir.mkdir("test_gemmi_complex_metadata")

    f = sr.save(mols, d.join("test.pdbx"))

    if isinstance(f, list):
        f = f[0]

    mols2 = sr.load(f)

    m = mols2.metadata()

    assert m["name"] == "alanine dipeptide"
    assert m["residues"] == ["ACE", "ALA", "NME"]

    a = m["atoms"]

    assert a["element"].as_array() == ["C", "N", "O"]
    assert a["x_coords"].as_array() == ["0", "1", "2"]
    assert a["y_coords"].as_array() == ["3", "4", "5"]
    assert a["z_coords"].as_array() == ["6", "7", "8"]
