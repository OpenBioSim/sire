import pytest
import sire as sr


@pytest.mark.skipif(
    "rdkit" not in sr.convert.supported_formats(),
    reason="rdkit support is not available",
)
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

    m = mols["smarts [#6]"]

    assert m.num_groups() == mols["element C"].num_atoms()

    for group in m.groups():
        assert group.num_atoms() == 1
        assert group[0].element().num_protons() == 6

    counts = {}

    # The expected values have been precomputed using rdkit directly
    # via the script
    #
    # from rdkit import Chem
    #
    # some_chemicals =[
    #                 'CC(=O)Nc1ccc(O)cc1',
    #                 'CC(C)NCC(O)COc1ccccc1CC=C',
    #                 'CC(N)Cc1ccccc1',
    #                 'CC(CS)C(=O)N1CCCC1C(=O)O',
    #                 'CN(C)CCCN1c2ccccc2Sc3ccc(Cl)cc13',
    #                 'OC(=O)Cc1ccccc1Nc2c(Cl)cccc2Cl',
    #                 'NCC1(CC(=O)O)CCCCC1',
    #                 'COC(=O)c1ccccc1O',
    #                 'Nc1ccc(N=Nc2ccccc2)c(N)n1',
    #                 'IC(=O)c1ccccc1',
    #                 'CCOP(=S)(OCC)Oc1cc(Cl)cc(Cl)c1',
    #                 'c1c(C)c(O)c(N)cc1',
    #                 'Oc1c(C)cc(N)cc1',
    #                 'Oc1c(C)ccc(N)c1',
    #                 'c1c(C)c(N)c(O)cc1',
    #                 ]
    #
    # some_chemicals = list(map(Chem.MolFromSmiles, some_chemicals))
    #
    # smart_mol = Chem.MolFromSmarts('[#6]!:[#6]')
    # highlightAtomLists = [len([atom[0] for atom in mol.GetSubstructMatches(smart_mol)]) for mol in some_chemicals]
    #
    # print(highlightAtomLists)

    expected = [1, 7, 3, 7, 2, 2, 9, 1, 0, 1, 2, 1, 1, 1, 1]

    for i, molnum in enumerate(molnums):
        counts[molnum] = expected[i]

    m = mols["smarts [#6]!:[#6]"]

    for group in m.groups():
        counts[group.molecule().number()] -= 1

        assert group.num_atoms() == 2
        assert group[0].element().num_protons() == 6
        assert group[1].element().num_protons() == 6

    for count in counts.values():
        assert count == 0
