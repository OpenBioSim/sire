__all__ = ["match_atoms"]


def match_atoms(
    mol0, mol1, match=None, prematch=None, match_light_atoms=False, map0=None, map1=None
):
    """
    Perform a simple match that tries to identify the mapping from
    atoms in 'mol0' to the atoms in 'mol1'.
    """
    from .mol import AtomMapping
    from .legacy.Mol import AtomMCSMatcher, AtomMatcher
    from .base import create_map

    map0 = create_map(map0)
    map1 = create_map(map1)

    if match is not None:
        if not isinstance(match, AtomMatcher):
            from .legacy.Mol import AtomIDMatcher

            matcher = AtomIDMatcher(match)
        else:
            matcher = match

    elif prematch is not None:
        if not isinstance(prematch, AtomMatcher):
            from .legacy.Mol import AtomIDMatcher

            prematch = AtomIDMatcher(prematch)

        matcher = AtomMCSMatcher(
            prematcher=prematch, match_light_atoms=match_light_atoms, verbose=False
        )
    else:
        matcher = AtomMCSMatcher(match_light_atoms=match_light_atoms, verbose=False)

    m = matcher.match(mol0, map0, mol1, map1)

    atoms0 = []
    atoms1 = []

    for atom0, atom1 in m.items():
        atoms0.append(atom0.value())
        atoms1.append(atom1.value())

    return AtomMapping(
        mol0.atoms(),
        mol1.atoms(),
        mol0.molecule().atoms()[atoms0],
        mol1.molecule().atoms()[atoms1],
        map0,
        map1,
    )
