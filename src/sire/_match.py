__all__ = ["match_atoms"]


def match_atoms(mol0, mol1, match_light_atoms=False, map0=None, map1=None):
    """
    Perform a simple match that tries to identify the mapping from
    atoms in 'mol0' to the atoms in 'mol1'.
    """
    from .legacy.Mol import AtomMCSMatcher
    from .base import create_map

    map0 = create_map(map0)
    map1 = create_map(map1)

    matcher = AtomMCSMatcher(
        match_light_atoms=match_light_atoms, verbose=False
    )

    m = matcher.match(mol0, map0, mol1, map1)

    match = {}

    for atom0, atom1 in m.items():
        match[atom0.value()] = atom1.value()

    # now sort this dictionary
    keys = list(match.keys())
    keys.sort()

    match2 = {}

    for key in keys:
        match2[key] = match[key]

    return match2
