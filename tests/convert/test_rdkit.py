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


def test_chirality():
    l_ala = sr.smiles("C[C@H](N)C(=O)O")
    d_ala = sr.smiles("C[C@@H](N)C(=O)O")

    assert d_ala["C2"].property("chirality").is_clockwise()
    assert l_ala["C2"].property("chirality").is_counter_clockwise()
