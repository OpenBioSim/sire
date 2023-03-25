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
        "OCC(O)C(O)C(O)C(O)CO",
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
