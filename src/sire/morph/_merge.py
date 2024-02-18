__all__ = ["merge"]


from ..mol import AtomMapping as _AtomMapping


def _merge(mapping: _AtomMapping, map=None):
    from ..legacy.System import merge as _merge

    return _merge(mapping, map)


def merge(mol0, mol1, map=None):
    mapping = _AtomMapping(mol0.atoms(), mol1.atoms())
    return _merge(mapping, map=map)


if not hasattr(_AtomMapping, "merge"):
    _AtomMapping.merge = _merge
