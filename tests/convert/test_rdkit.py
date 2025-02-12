import pytest
import sire as sr


@pytest.mark.skipif(
    "rdkit" not in sr.convert.supported_formats(),
    reason="rdkit support is not available",
)
@pytest.mark.parametrize(
    "smiles",
    [
        "C1CCCCC1",
        "C",
        "OC[C@H](O)[C@H](O)[C@@H](O)[C@@H](O)CO",
        "C[C@H](N)C(=O)O",  # L-alanine
        "C[C@@H](N)C(=O)O",  # D-alanine
    ],
)
def test_rdkit(smiles):
    mol = sr.smiles(smiles)
    assert mol.smiles() == smiles


@pytest.mark.skipif(
    "rdkit" not in sr.convert.supported_formats(),
    reason="rdkit support is not available",
)
def test_chirality():
    l_ala = sr.smiles("C[C@H](N)C(=O)O")
    d_ala = sr.smiles("C[C@@H](N)C(=O)O")

    assert d_ala["C2"].property("chirality").is_clockwise()
    assert l_ala["C2"].property("chirality").is_counter_clockwise()


@pytest.mark.skipif(
    "rdkit" not in sr.convert.supported_formats(),
    reason="rdkit support is not available",
)
def test_hybridization():
    l_ala = sr.smiles("C[C@H](N)C(=O)O")

    for c in l_ala["element carbon"]:
        bonds = c["bonds to *"]

        if len(bonds) == 3:
            assert c.property("hybridization").is_sp2()
        elif len(bonds) == 4:
            assert c.property("hybridization").is_sp3()
        else:
            print("WEIRD BONDS", c, bonds)
            assert False

    for n in l_ala["element nitrogen"]:
        bonds = n["bonds to *"]

        if len(bonds) == 2:
            assert n.property("hybridization").is_sp2()
        elif len(bonds) == 3:
            assert n.property("hybridization").is_sp3()
        else:
            print("WEIRD BONDS", n, bonds)
            assert False


@pytest.mark.skipif(
    "rdkit" not in sr.convert.supported_formats(),
    reason="rdkit support is not available",
)
def test_hybridization2():
    mol = sr.smiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")

    expected = [
        "SP3",
        "SP2",
        "SP2",
        "SP2",
        "SP2",
        "SP2",
        "SP2",
        "SP2",
        "SP2",
        "SP2",
        "SP2",
        "SP2",
        "SP3",
        "SP3",
    ]

    for atom, e in zip(mol.atoms(), expected):
        assert atom.property("hybridization").to_rdkit() == e


@pytest.mark.skipif(
    "rdkit" not in sr.convert.supported_formats(),
    reason="rdkit support is not available",
)
def test_rdkit_returns_null():
    # This molecule should construct if sanitization is turned off
    mol = sr.smiles("c3cc[c+]2cc(C1CCCC1)[nH]c2c3", must_sanitize=False)

    assert mol.num_atoms() == 36

    # coordinates should be possible, albeit likely wrong!
    assert mol.has_property("coordinates")

    # We should raise an exception if sanitization fails
    # and `must_sanitize` is True
    with pytest.raises(ValueError):
        mol = sr.smiles("c3cc[c+]2cc(C1CCCC1)[nH]c2c3", must_sanitize=True)

    # this should be the default
    with pytest.raises(ValueError):
        mol = sr.smiles("c3cc[c+]2cc(C1CCCC1)[nH]c2c3")


@pytest.mark.skipif(
    "rdkit" not in sr.convert.supported_formats(),
    reason="rdkit support is not available",
)
def test_rdkit_infer_bonds(ejm55_sdf, ejm55_gro):
    sdf = ejm55_sdf[0].molecule()
    gro = ejm55_gro["not (protein or water)"].molecule()

    from rdkit import Chem

    assert sdf.smiles() == gro.smiles()

    match_sdf = sdf["smarts [NX3][CX3](=[OX1])[#6]"]
    match_gro = gro["smarts [NX3][CX3](=[OX1])[#6]"]

    print(match_sdf)
    print(match_gro)

    assert len(match_sdf) == 1
    assert len(match_gro) == 1

    for s, g in zip(match_sdf, match_gro):
        assert s.number() == g.number()


@pytest.mark.skipif(
    "rdkit" not in sr.convert.supported_formats(),
    reason="rdkit support is not available",
)
def test_rdkit_preserve_info(ala_mols, ejm55_gro):
    mol0 = ala_mols[0]
    mol1 = ejm55_gro["not (protein or water)"].molecule()

    r0 = sr.convert.to(mol0, "rdkit")
    r1 = sr.convert.to(mol1, "rdkit")

    m0 = sr.convert.to(r0, "sire")
    m1 = sr.convert.to(r1, "sire")

    for mol, m in [(mol0, m0), (mol1, m1)]:
        for res0, res1 in zip(mol.residues(), m.residues()):
            assert res0.name() == res1.name()
            assert res0.number() == res1.number()

        for atom0, atom1 in zip(mol.atoms(), m.atoms()):
            assert atom0.name() == atom1.name()
            assert atom0.number() == atom1.number()

            res0 = atom0.residue()
            res1 = atom1.residue()

            assert res0.name() == res1.name()
            assert res0.number() == res1.number()
