__all__ = []

# Imported so that it is pythonized
from ..legacy import Qt as _Qt  # noqa: F401

from .. import use_new_api as _use_new_api

_use_new_api()
