import pytest
import sire as sr


@pytest.mark.skipif(
    "rdkit" not in sr.convert.supported_formats(),
    reason="rdkit support is not available",
)
@pytest.mark.parametrize("smiles", ["C1CCCCC1", "C", "OCC(O)C(O)C(O)C(O)CO"])
def test_rdkit(smiles):
    mol = sr.smiles(smiles)
    assert mol.smiles() == smiles
