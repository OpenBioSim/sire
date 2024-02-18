__all__ = ["mutate"]


from ..mol import AtomMapping as _AtomMapping


def _mutate(mapping: _AtomMapping, map=None):
    from ..legacy.System import mutate as _mutate

    return _mutate(mapping, map)


def mutate(mol0, mol1, map=None):
    mapping = _AtomMapping(mol0.atoms(), mol1.atoms())
    return _mutate(mapping, map=map)


if not hasattr(_AtomMapping, "mutate"):
    _AtomMapping.mutate = _mutate
