__all__ = ["merge"]


from ..mol import AtomMapping as _AtomMapping


def _merge(mapping: _AtomMapping, as_new_molecule: bool = True, map=None):
    from ..legacy.System import merge as _merge_mols
    from ..base import create_map

    map = create_map(map)

    # now align the perturbed state onto the reference state,
    # so that any added atoms have roughly the right coordinates
    aligned_mapping = mapping.align()

    return _merge_mols(aligned_mapping, as_new_molecule=as_new_molecule, map=map)


def merge(mol0, mol1, match=None, prematch=None, map=None, map0=None, map1=None):
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
