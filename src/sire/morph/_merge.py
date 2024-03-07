__all__ = ["merge"]


from ..mol import AtomMapping as _AtomMapping


def _merge(mapping: _AtomMapping, as_new_molecule: bool = True, map=None):
    """
    Merge the atoms in this mapping and return as a single merged
    (perturbable) molecule. This function will conduct a merge and
    return a perturbable molecule such that it is equivalent to
    `mol0` at the reference state (lambda=0) and equivalent to `mol1`
    at the perturbed state (lambda=1).

    Parameters
    ----------

    as_new_molecule : bool, optional
        If True, the merged molecule will be assigned a new molecule
        number and treated as a new molecule. If False, the merged
        molecule will use the molecule number of the reference molecule.
    map : dict, optional
        Property map to assign properties in the returned,
        merged molecule.

    Returns
    -------
    Molecule
        The merged molecule
    """
    from ..legacy.System import merge as _merge_mols
    from ..base import create_map

    map = create_map(map)

    # now align the perturbed state onto the reference state,
    # so that any added atoms have roughly the right coordinates
    aligned_mapping = mapping.align()

    mol = _merge_mols(aligned_mapping, as_new_molecule=as_new_molecule, map=map)

    mol = mol.perturbation().link_to_reference()

    return mol


def merge(mol0, mol1, match=None, prematch=None, map=None, map0=None, map1=None):
    """
    Merge together the atoms in 'mol0' and 'mol1' and return as a single
    merged (perturbable) molecule. This function will conduct a merge
    and return a perturbable molecule such that it is equivalent to
    `mol0` at the reference state (lambda=0) and equivalent to `mol1` at
    the perturbed state (lambda=1).

    The `sr.morph.match_atoms` function will be called with the passed
    `match` and `prematch` arguments to determine the atom mapping between
    the two molecules.

    Parameters
    ----------
    mol0 : Molecule view
        The reference state molecule (or part of molecule)
    mol1 : Molecule view
        The perturbed state molecule (or part of molecule)
    match : dict, AtomMapping, optional
        If provided, this will be passed as the `match` argument
        to `sr.morph.match_atoms`, to aid in the atom mapping.
    prematch : dict, AtomMapping, optional
        If provided, this will be passed as the `prematch` argument
        to `sr.morph.match_atoms`, to aid in the atom mapping.
    map : dict, optional
        Property map to assign properties in the returned,
        merged molecule.
    map0 : dict, optional
        Property map to find properties in `mol0`
    map1 : dict, optional
        Property map to find properties in `mol1`

    Returns
    -------
    Molecule
        The merged molecule
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

    from . import match as _match

    mapping = _match(
        mol0=mol0,
        mol1=mol1,
        match=match,
        prematch=prematch,
        match_light_atoms=True,
        map0=map0,
        map1=map1,
    )

    return mapping.merge(as_new_molecule=True, map=map)


if not hasattr(_AtomMapping, "merge"):
    _AtomMapping.merge = _merge
