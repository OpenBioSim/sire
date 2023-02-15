import sire as sr


def test_excluded(excluded_mols):
    mols = excluded_mols

    mol = mols[0]

    e = sr.legacy.MM.ExcludedPairs(mol)

    assert e.num_excluded_pairs() == 0
