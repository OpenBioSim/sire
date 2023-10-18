__all__ = []

# Need to import the module to pythonize it
from ..legacy import Analysis as _Analysis  # noqa: F401

from .. import use_new_api as _use_new_api

_use_new_api()
