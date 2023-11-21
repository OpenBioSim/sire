__all__ = ["EMLEEngine"]

from ..legacy import Convert as _Convert

from .. import use_new_api as _use_new_api

_use_new_api()

EMLEEngine = _Convert._SireOpenMM.EMLEEngine
