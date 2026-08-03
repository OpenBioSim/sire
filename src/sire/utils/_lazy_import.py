"""
A minimal, standard-library-only replacement for the third-party
`lazy_import` package (GPLv3), used by sire/__init__.py to defer importing
its submodules until they are actually used.

This is a private module - not part of sire.utils's public API, just an
internal bootstrapping helper.
"""

import importlib.util as _importlib_util
import sys as _sys
import threading as _threading
import types as _types

__all__ = ["lazy_module", "is_lazy_module", "force_load"]


class _LazyModule(_types.ModuleType):
    """
    A stand-in for a module that hasn't been imported yet. The first
    time any real attribute is accessed (or force_load() is called), the
    actual module is imported and this proxy is replaced everywhere it's
    reachable (sys.modules, and the attribute of whatever parent package
    holds it) with the real thing.
    """

    def __init__(self, name: str):
        super().__init__(name)
        # stored directly in __dict__, so accessing these never itself
        # goes through __getattr__ below (which would recurse)
        self._sire_lazy_name = name
        self._sire_lazy_loaded = False
        # guards _sire_load() below - without this, two threads racing to
        # first-touch the same lazy module could each build and exec_module()
        # their *own* separate module object concurrently, corrupting
        # whatever global C++-side registration state that module's import
        # touches (RegisterMetaType and friends aren't designed to run
        # twice at once)
        self._sire_lazy_lock = _threading.Lock()

    def _sire_load(self):
        if self._sire_lazy_loaded:
            return _sys.modules[self._sire_lazy_name]

        with self._sire_lazy_lock:
            # re-check now that we hold the lock - another thread may have
            # already finished loading while we were waiting for it
            if self._sire_lazy_loaded:
                return _sys.modules[self._sire_lazy_name]

            self._sire_lazy_loaded = True

            # This proxy is itself sitting in sys.modules[name] right now.
            # find_spec()/import_module() both consult sys.modules first, so
            # if we left it there they'd either hand this same proxy straight
            # back (recursing forever) or fail (since it has no real
            # __spec__) - pop it out first so the lookup below does a genuine
            # fresh search.
            _sys.modules.pop(self._sire_lazy_name, None)

            # Build the real module and register it in sys.modules *before*
            # executing its body - exactly what Python's own import system
            # does, and important for the same reason: if executing the
            # module's body re-entrantly imports this same name (a circular
            # import), that lookup must see this (partially-initialised)
            # real module, not a second independent load.
            spec = _importlib_util.find_spec(self._sire_lazy_name)
            real_module = _importlib_util.module_from_spec(spec)

            _sys.modules[self._sire_lazy_name] = real_module

            # also fix up the attribute on the parent package, if there is
            # one, so that e.g. `sire.maths` now points at the real module
            parent_name, _, attr_name = self._sire_lazy_name.rpartition(".")

            if parent_name:
                parent = _sys.modules.get(parent_name)

                if parent is not None:
                    setattr(parent, attr_name, real_module)

            spec.loader.exec_module(real_module)

            return real_module

    def __getattr__(self, attr):
        # only real attribute lookups (i.e. ones that fail against this
        # proxy's own __dict__) should trigger the load - this deliberately
        # does *not* override __dir__/__repr__, so introspection (e.g. by
        # a debugger, or pytest's failure reporting) doesn't silently
        # trigger a real import as a side effect
        if attr.startswith("_sire_lazy"):
            raise AttributeError(attr)

        real_module = self._sire_load()
        return getattr(real_module, attr)

    def __repr__(self):
        if self._sire_lazy_loaded:
            return repr(_sys.modules[self._sire_lazy_name])

        return f"<lazily-loaded module {self._sire_lazy_name!r}>"


def lazy_module(name: str):
    """
    Return a proxy for the module called 'name' that defers the actual
    import until one of its attributes is accessed, or force_load() is
    called on it explicitly. This is a drop-in, GPLv3-free replacement for
    lazy_import.lazy_module().
    """
    existing = _sys.modules.get(name)

    if existing is not None and not isinstance(existing, _LazyModule):
        return existing

    module = _LazyModule(name)
    _sys.modules[name] = module

    return module


def is_lazy_module(module) -> bool:
    """Return whether 'module' is a not-yet-loaded lazy module proxy."""
    return isinstance(module, _LazyModule) and not module._sire_lazy_loaded


def force_load(module):
    """
    If 'module' is a lazy module proxy, force it to load now. This is
    a no-op if 'module' is already a real, fully-loaded module.
    """
    if isinstance(module, _LazyModule):
        module._sire_load()
