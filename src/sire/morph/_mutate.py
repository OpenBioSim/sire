__all__ = ["mutate"]


from ..mol import AtomMapping as _AtomMapping


def _mutate(mapping: _AtomMapping, as_new_molecule: bool = True, map=None):
    """
    Mutate the reference atoms in this mapping to the perturbed atoms,
    returning the mutated (new) molecule.

    This is equivalent to calling `merge` and then extracting the
    perturbed state from the returned merged molecule.

    This function is most useful for mutating parts of molecules,
    e.g. calling this on a mapping of two residues would mutate
    one residue into another within the larger molecule containing
    the reference mapping. This can be used for mutating residues
    in proteins, or for copying and pasting parts of one molecule
    into another.

    Parameters
    ----------
    as_new_molecule : bool, optional
        Whether to return the mutated molecule as a new molecule,
        or to mutate the original molecule in place. Default is True.
    map : dict, optional
        Property map to assign properties in the returned,
        mutated molecule.

    Returns
    -------
    Molecule
        The mutated molecule
    """
    return (
        mapping.merge(as_new_molecule=as_new_molecule, map=map)
        .perturbation()
        .extract_perturbed(remove_ghosts=True)
    )


def mutate(mol0, mol1, match=None, prematch=None, map=None, map0=None, map1=None):
    """
    Mutate `mol0` to `mol1`, returning the mutated (new) molecule.
    This is equivalent to calling `merge` on the two molecules (or
    parts of molecules) and then extracting the perturbed state.

    This function is most useful for mutating parts of molecules,
    e.g. passing in two residues as `mol0` and `mol1` would mutate
    that residue to the other within the larger molecule containing
    `mol0`. This can be used for mutating residues in proteins, or
    for copying and pasting parts of one molecule into another.

    Parameters
    ----------
    mol0 : Molecule view
        The molecule (or part of molecule) that will be mutated.
    mol1 : Molecule view
        The molecule (or part of molecule) that will be mutated to.
        This will replace the atoms in `mol0`.
    match : dict, AtomMatcher, optional
        If provided, this will be passed as the `match` argument
        to `sr.morph.match_atoms`, to aid in the atom mapping.
    prematch : dict, AtomMatcher, optional
        If provided, this will be passed as the `prematch` argument
        to `sr.morph.match_atoms`, to aid in the atom mapping.
    map : dict, optional
        Property map to assign properties in the returned,
        mutated molecule.
    map0 : dict, optional
        Property map to find properties in `mol0`
    map1 : dict, optional
        Property map to find properties in `mol1`

    Returns
    -------
    Molecule
        The mutated molecule
    """
    from ..base import create_map

    map = create_map(map)

    if map0 is None:
        map0 = map
    else:
        map0 = create_map(map, map0)

    if map1 is None:
        map1 = map
    else:
        map1 = create_map(map, map1)

    from . import match

    mapping = match(mol0=mol0, mol1=mol1, match_light_atoms=True, map0=map0, map1=map1)

    return mapping.mutate(as_new_molecule=True, map=map)


if not hasattr(_AtomMapping, "mutate"):
    _AtomMapping.mutate = _mutate
