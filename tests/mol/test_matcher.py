def test_atomcoordmatcher(ala_mols):
    """
    Test the AtomCoordMatcher class.
    """
    from sire.legacy.Mol import AtomCoordMatcher

    matcher = AtomCoordMatcher()

    # Match the first two molecules two each other.
    matches = matcher.match(ala_mols[0], ala_mols[0])

    # Assert that all atoms are matched to themselves.
    for k, v in matches.items():
        assert k == v
