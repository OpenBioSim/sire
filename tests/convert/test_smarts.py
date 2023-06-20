import pytest
import sire as sr


def test_smarts():
    # Thanks to this blog post which really helped me understand
    # how smiles worked (and from where these below examples are copied)
    # https://russodanielp.github.io/blog/a-brief-introduction-to-smarts/

    some_chemicals = [
        "CC(=O)Nc1ccc(O)cc1",
        "CC(C)NCC(O)COc1ccccc1CC=C",
        "CC(N)Cc1ccccc1",
        "CC(CS)C(=O)N1CCCC1C(=O)O",
        "CN(C)CCCN1c2ccccc2Sc3ccc(Cl)cc13",
        "OC(=O)Cc1ccccc1Nc2c(Cl)cccc2Cl",
        "NCC1(CC(=O)O)CCCCC1",
        "COC(=O)c1ccccc1O",
        "Nc1ccc(N=Nc2ccccc2)c(N)n1",
        "IC(=O)c1ccccc1",
        "CCOP(=S)(OCC)Oc1cc(Cl)cc(Cl)c1",
        "c1c(C)c(O)c(N)cc1",
        "Oc1c(C)cc(N)cc1",
        "Oc1c(C)ccc(N)c1",
        "c1c(C)c(N)c(O)cc1",
    ]

    mols = []
    molnums = []

    for mol in some_chemicals:
        mol = sr.smiles(mol)
        mols.append(mol)
        molnums.append(mol.number())

    mols = sr.mol.SelectorMol(mols)

    for match in mols["smarts [#6]"]:
        for group in match.groups():
            assert group.num_atoms() == 1
            assert group[0].element().num_protons() == 6

    counts = {}

    expected = [2, 10, 4, 9, 5, 3, 9, 3, 0, 2, 4, 2, 2, 2, 2]

    for i, molnum in enumerate(molnums):
        counts[molnum] = expected[i]

    for atom in mols["smarts [#6]!:[#6]"]:
        print(atom)
        counts[atom.molecule().number()] -= 1

    print(counts)

    for count in counts.values():
        assert count == 0

    assert False
