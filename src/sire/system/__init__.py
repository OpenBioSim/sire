__all__ = [
    "System",
    "ForceFieldInfo",
    "calculate_energy",
    "create_forcefield",
]

from ..legacy import System as _System
from .. import use_new_api as _use_new_api

from ._system import System


class ForceFieldInfo(_System.ForceFieldInfo):
    def __init__(self, val, map=None):
        from ..base import create_map

        map = create_map(map)

        if hasattr(val, "_system"):
            val = val._system

        super().__init__(val, map=map)


_use_new_api()


def calculate_energy(*args, **kwargs):
    from ..mol import _to_molecules

    new_args = []
    new_kwargs = {}

    for arg in args:
        try:
            new_args.append(_to_molecules(arg))
        except Exception:
            new_args.append(arg)

    for key, value in kwargs.items():
        if key == "map":
            from ..base import create_map

            new_kwargs[key] = create_map(value)
        else:
            try:
                new_kwargs[key] = _to_molecules(value)
            except Exception:
                new_kwargs[key] = value

    return _System.calculate_energy(*new_args, **new_kwargs)


def create_forcefield(*args, map=None, **kwargs):
    from ..mol import _to_molecules

    new_args = []
    new_kwargs = {}

    for arg in args:
        try:
            new_args.append(_to_molecules(arg))
        except Exception:
            new_args.append(arg)

    for key, value in kwargs.items():
        try:
            new_kwargs[key] = _to_molecules(value)
        except Exception:
            new_kwargs[key] = value

    from ..base import create_map

    map = create_map(map)

    new_kwargs["map"] = map

    return _System.create_forcefield(*new_args, **new_kwargs)
