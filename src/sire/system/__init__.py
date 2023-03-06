__all__ = ["System", "ForceFieldInfo"]

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
