__all__ = ["match_atoms"]


def match_atoms(
    mol0, mol1, match=None, prematch=None, match_light_atoms=False, map0=None, map1=None
):
    """
    Perform a simple match that tries to identify the mapping from
    atoms in 'mol0' to the atoms in 'mol1'. This uses the `AtomMCSMatcher`
    to match the atoms, using the passed `prematch` argument.

    However, if the `match` argument is provided, this will be used
    as the atom mapping directly (it can either be a dictionary mapping
    atom identifiers, or an `AtomMatcher` object).

    Parameters
    ----------
    mol0 : Molecule view
        The reference state molecule (or part of molecule)
    mol1 : Molecule view
        The perturbed state molecule (or part of molecule)
    match : dict, AtomMatcher, optional
        The atom matcher to use to match atoms. If this is a dictionary
        of atom identifiers, then this will be passed to a
        `AtomIDMatcher` object. If this is an `AtomMatcher` object, then
        this will be used directly.
    prematch : dict, AtomMatcher, optional
        The atom matcher to use to prematch atoms. If `match` is not
        supplied, then this will be used as the `prematch` argument
        to the `AtomMCSMatcher` used to find the maximum common subgraph
        match.
    match_light_atoms : bool, optional
        Whether to match light atoms (i.e. hydrogen atoms) if using the
        default `AtomMCSMatcher`. Default is False.
    map0 : dict, optional
        Property map to find properties in `mol0`
    map1 : dict, optional
        Property map to find properties in `mol1`

    Returns
    -------
    AtomMapping
        The atom mapping between the two molecules (or parts of molecules)
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
