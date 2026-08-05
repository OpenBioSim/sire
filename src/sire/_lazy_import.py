"""
A minimal, standard-library-only replacement for the third-party
`lazy_import` package (GPLv3), used to defer importing submodules until
they are actually used. Package-agnostic: nothing here assumes anything
about the particular package it's wrapping, so the same module can be
(and is) shared verbatim between sire and BioSimSpace.

This is a private module - not part of the public API, just an internal
bootstrapping helper.
"""

import importlib.machinery as _importlib_machinery
import importlib.util as _importlib_util
import pkgutil as _pkgutil
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

    Known limitations:

    - This replacement only reaches places that look the name up again
      later (sys.modules, or a parent's own attribute) - a name already
      bound to *this* stub before the load happened (e.g. `import
      pkg.sub as x`, which binds `x` to whatever was in sys.modules at
      that moment) keeps pointing at the stub forever. Reads through
      that binding still work, via __getattr__ delegating to the real
      module below - but a write (`x.SOME_ATTR = value`) lands on the
      stub's own __dict__ instead, invisible to sys.modules["pkg.sub"]
      or anyone else holding a reference obtained afterwards.

    - importlib.reload() on a module that hasn't been touched at all
      yet doesn't work: reload() calls exec_module() directly on
      whatever object sys.modules already holds, bypassing this class's
      own _load() (and everything it does - the lock, the parent-
      attribute fixup, the child-stub sweep) entirely. Call
      force_load()/touch an attribute first, then reload() the (by then
      real) module as normal. This isn't specific to this
      implementation - any lazy-loading scheme built on a sys.modules
      stand-in has the same gap, since reload()'s entry point isn't one
      this class - or any other module object - gets a say in.
    """

    def __init__(self, name: str, search_locations=None, spec=None):
        super().__init__(name)

        # Populate the standard module attributes (__file__, __path__,
        # __loader__, __spec__, __package__) straight from the spec, if we
        # have one, exactly as a real (non-lazy) module would have them.
        # This matters because plenty of code that has no interest in this
        # module's contents still checks for these (e.g. inspect.getmodule()
        # scans *every* entry in sys.modules doing hasattr(module,
        # '__file__') while resolving a source location for something else
        # entirely - seen for real via a third-party library's use of
        # inspect deep inside an unrelated import). Without this, that one
        # hasattr() call would fall through to __getattr__ below and force
        # a full, real load - of every single lazily-registered module in
        # the process, all at once, at whatever arbitrary moment such a
        # scan happens to run.
        #
        # These are read directly off the spec, not by building a
        # module_from_spec(spec) and copying its __dict__: ModuleType's
        # own __init__() above has already pre-seeded __spec__/__loader__/
        # __package__/__doc__ as real dict entries set to None, so copying
        # via dict.setdefault() (which a template's __dict__ would need,
        # to avoid clobbering __name__) would silently no-op for exactly
        # those three. And separately, module_from_spec() calls the
        # loader's create_module(), which for a single-phase-init C
        # extension (what Boost.Python wrappers are) actually runs the
        # module's init function - dlopening it and executing its C-level
        # registration code. Building a template purely to copy attributes
        # off it would do that at registration time, for every compiled
        # submodule reachable in the lazy tree, then a second time when
        # something actually imports it properly - exactly the "twice is
        # corruption" scenario the lock elsewhere in this class exists to
        # prevent. Reading spec attributes directly avoids calling
        # module_from_spec() at registration time at all.
        if spec is not None:
            self.__dict__["__spec__"] = spec
            self.__dict__["__loader__"] = spec.loader
            self.__dict__["__package__"] = spec.parent

            # spec.has_location, not "spec.origin is not None" - that's
            # what CPython's own _init_module_attrs uses to decide this,
            # since some loaders put a non-path marker in origin (e.g.
            # 'frozen', 'built-in') rather than leaving it None. Not
            # reachable via PathFinder for anything in sire's own tree,
            # but __file__ is exactly what the inspect.getmodule
            # protection above leans on - a wrong value there is worse
            # than a merely absent one.
            if spec.has_location:
                self.__dict__["__file__"] = spec.origin

            if spec.submodule_search_locations is not None:
                # list(...) here deliberately freezes a namespace
                # package's _NamespacePath at this moment, rather than
                # leaving it as the dynamic, sys.path-tracking object a
                # real import would use - fine for sire (no namespace
                # packages in play), but worth knowing if that ever
                # changes.
                self.__dict__["__path__"] = list(spec.submodule_search_locations)

        # stored directly in __dict__, so accessing these never itself
        # goes through __getattr__ below (which would recurse)
        self._lazy_name = name
        # Keep the spec we were already handed (both call sites below -
        # _register_submodules() and lazy_module() - already had to look
        # it up to get this far) so _load() can reuse it directly instead
        # of doing a second, redundant find_spec() call at load time. See
        # _load() below for the one behavioural consequence of this: the
        # spec is now pinned at registration time rather than re-resolved
        # at first touch.
        self._lazy_spec = spec
        # The real, fully-executed module, once _load() has completed - not
        # just "have we started loading" (see _load() below for why that
        # distinction matters).
        self._lazy_real = None
        # Set to the loading thread's ident for the duration of _load()'s
        # actual work (spec lookup through exec_module()), and back to None
        # once it's done (successfully or not) - lets a *reentrant* call
        # from the same thread (a circular import touching this module from
        # inside its own exec_module()) through without deadlocking, while
        # still blocking any *other* thread until the load truly finishes.
        self._lazy_owner = None
        # If known (i.e. this is a submodule registered by
        # _register_submodules()), the parent's search path - see _load()
        # below for why this needs to be captured up front rather than
        # derived when the module actually loads.
        self._lazy_search_locations = search_locations
        # Immediate children pre-registered under this module's name by
        # _register_submodules() - recorded up front purely so _load()'s
        # post-exec fixup below has an exact, cheap list to walk, rather
        # than needing to scan the whole of sys.modules (which can run to
        # hundreds of entries) on every single module load.
        self._lazy_children = []
        # guards _load() below - without this, two threads racing to
        # first-touch the same lazy module could each build and exec_module()
        # their *own* separate module object concurrently, corrupting
        # whatever global C++-side registration state that module's import
        # touches (RegisterMetaType and friends aren't designed to run
        # twice at once)
        self._lazy_lock = _threading.Lock()

    def _load(self):
        # Fully loaded already - safe to fast-path with no lock, since
        # _lazy_real is only ever set (by the thread that did the loading)
        # *after* exec_module() has completely finished (see below), so
        # every other thread reading it here is guaranteed to see a fully
        # executed module, never a partial one.
        if self._lazy_real is not None and self._lazy_owner is None:
            return self._lazy_real

        # Reentrant call from the thread currently doing the load (a
        # circular import reached back into this module from inside its
        # own exec_module()) - hand back the partially-initialised module,
        # exactly what CPython itself does for ordinary circular imports.
        # Taking the lock here would deadlock against ourselves.
        if self._lazy_owner == _threading.get_ident():
            if self._lazy_real is None:
                # Re-entered before the module object itself exists yet -
                # only reachable between this thread entering the lock
                # below and self._lazy_real being set further down, a
                # window module_from_spec()'s own create_module() can
                # reach into for a single-phase-init C extension (its
                # module init function can run arbitrary C-level code,
                # including imports, before Python-level exec even
                # starts). There's nothing to hand back yet - and nothing
                # useful to build, since a second, independent module
                # object here would defeat the point of this whole class
                # - so raise something a caller stands a chance of
                # diagnosing, rather than the AttributeError on None that
                # returning self._lazy_real as-is would produce a few
                # lines down.
                raise ImportError(
                    f"circular import: {self._lazy_name!r} was re-entered "
                    "before its module object existed"
                )

            return self._lazy_real

        with self._lazy_lock:
            # re-check now that we hold the lock - another thread may have
            # already finished loading while we were waiting for it
            if self._lazy_real is not None and self._lazy_owner is None:
                return self._lazy_real

            self._lazy_owner = _threading.get_ident()

            # Computed up front, before anything that can raise - the
            # except handler below needs both, and must never risk an
            # UnboundLocalError masking the real failure if it happens
            # before this point is otherwise reached.
            parent_name, _, attr_name = self._lazy_name.rpartition(".")
            parent_patched = False

            try:
                # Both places that construct a _LazyModule (_register_
                # submodules() and lazy_module()) already had to look up a
                # spec to get this far, and handed it to us in __init__ -
                # reuse it rather than paying for a second, redundant
                # lookup here, which also makes the fallback branch below
                # (and its PathFinder-vs-find_spec choice) rare rather than
                # eliminating it - lazy_module() still allows a stub to be
                # built with spec=None if its own find_spec() failed, and
                # that's handled there.
                #
                # This is a small behaviour change, not just a speedup:
                # the spec is now resolved once, at *registration* time
                # (effectively, at `import sire` time for everything under
                # it), rather than at first-touch. A sys.path edit or
                # importlib.invalidate_caches() between those two points
                # no longer affects what actually gets loaded - it matches
                # what `import sire` itself already saw, which is more
                # correct, but it is a real difference from before.
                if self._lazy_spec is not None:
                    spec = self._lazy_spec

                    # module_from_spec() never consults sys.modules, and
                    # sys.modules[self._lazy_name] gets overwritten with
                    # real_module a few lines below regardless - so unlike
                    # the fallback branch, there's nothing here that reads
                    # sys.modules and needs this name popped out of it
                    # first. Skipping the pop also narrows the window
                    # where the name is absent from sys.modules entirely
                    # to nothing, rather than widening it.
                else:
                    # This proxy is itself sitting in sys.modules[name]
                    # right now. find_spec()/import_module() both consult
                    # sys.modules first, so if we left it there they'd
                    # either hand this same proxy straight back (recursing
                    # forever) or fail (since it has no real __spec__) -
                    # pop it out first so the lookup below does a genuine
                    # fresh search.
                    _sys.modules.pop(self._lazy_name, None)

                    # Use PathFinder directly with the parent's
                    # already-known search path, rather than
                    # importlib.util.find_spec(name) - for a dotted name,
                    # that variant automatically imports the parent
                    # package to get its __path__, which would re-enter a
                    # still-lazy parent's load from inside this module's
                    # own load.
                    if self._lazy_search_locations is not None:
                        spec = _importlib_machinery.PathFinder.find_spec(
                            self._lazy_name, self._lazy_search_locations
                        )
                    else:
                        spec = _importlib_util.find_spec(self._lazy_name)

                if spec is None:
                    raise ImportError(f"No module named {self._lazy_name!r}")

                real_module = _importlib_util.module_from_spec(spec)

                _sys.modules[self._lazy_name] = real_module

                # also fix up the attribute on the parent package, if there
                # is one, so that e.g. `sire.maths` now points at the real
                # module
                if parent_name:
                    parent = _sys.modules.get(parent_name)

                    if parent is not None:
                        setattr(parent, attr_name, real_module)
                        parent_patched = True

                # Record the (still-executing) real module *before* running
                # its body - exactly what Python's own import system does,
                # and important for the same reason: if executing the
                # module's body re-entrantly imports this same name (a
                # circular import), the reentrant branch above must see
                # this partially-initialised real module, not a second,
                # independent one.
                self._lazy_real = real_module

                spec.loader.exec_module(real_module)

                # Fix up every pre-registered child of *this* module that
                # its own body didn't already set as a real attribute (see
                # _register_submodules() for why that gap exists in the
                # first place). Doing this here, immediately after this
                # module's own exec_module() call, rather than up-front at
                # registration time, matters: for any access that goes
                # through the import system (`import parent.child`, or a
                # bare attribute access on an already-imported parent),
                # this preserves the guarantee that touching a child still
                # forces its parent's __init__.py to run first, in its
                # normal top-to-bottom order, exactly as a real (non-lazy)
                # import would - some packages' own submodules are only
                # safe to import in a specific order (e.g. because of
                # shared state in a third-party library they both touch),
                # and that guarantee would silently break if a child could
                # be reached directly off this module without ever running
                # its body. Note this guarantee comes from CPython's import
                # statement resolving the parent first, not from anything
                # this module does - a lookup that bypasses the import
                # system entirely (e.g. sys.modules["parent.child"]
                # directly) reaches the child's own stub without the
                # parent ever having been touched, same as it always could.
                #
                # Walk self._lazy_children (recorded up front by
                # _register_submodules()) rather than scanning the whole of
                # sys.modules - cheap, and exact.
                prefix_len = len(self._lazy_name) + 1

                for _child_name in self._lazy_children:
                    _child_module = _sys.modules.get(_child_name)

                    if _child_module is None:
                        continue

                    _child_attr = _child_name[prefix_len:]

                    if _child_attr not in real_module.__dict__:
                        setattr(real_module, _child_attr, _child_module)

                # Mirror the real module's own content onto this stub's
                # __dict__ too, not just sys.modules and the parent's
                # attribute. This matters for `from .sibling import *`
                # where 'sibling' has no __all__: CPython's IMPORT_STAR
                # opcode, for that case, reads getattr(mod, '__dict__')
                # directly - bypassing __getattr__ below entirely, since
                # every object always has a __dict__, so it's never
                # "missing" in the way that would trigger it - and
                # enumerates its keys. At that point 'mod' in the
                # caller's bytecode is still this stub, not real_module,
                # even though a preceding hasattr(mod, '__all__') check
                # just triggered this very _load() as a side effect (to
                # answer that hasattr). Without this, `from .sibling
                # import *` for any __all__-less module reached through
                # a lazy stub would silently import nothing at all.
                for _attr, _value in real_module.__dict__.items():
                    if _attr.startswith("_lazy"):
                        continue

                    self.__dict__[_attr] = _value
            except BaseException:
                # Leave this stub fully retryable - CPython itself leaves a
                # failed import retryable (a broken module, or one with a
                # missing optional dependency, doesn't get permanently
                # wedged just because someone touched it once), and without
                # this a failed load here would otherwise leave sys.modules
                # either missing the entry entirely (a confusing KeyError
                # on next use) or holding a half-initialised module that
                # looks loaded but silently lacks whatever wasn't reached
                # before the failure - in both cases masking the *original*
                # error on every subsequent attempt.
                #
                # Put *this stub* back in sys.modules, rather than just
                # leaving the name missing - a plain `import sire.mol`
                # retried after a failure would otherwise bypass it and go
                # through CPython's own machinery instead, reintroducing
                # the duplicate-module-object problem this whole scheme
                # exists to prevent.
                _sys.modules[self._lazy_name] = self

                # Likewise, undo the parent-attribute fixup above if it
                # already ran (it happens *before* exec_module(), so a
                # failure inside exec_module() would otherwise leave the
                # parent pointing at the broken module directly - a plain
                # attribute, not this stub - so retrying via the parent
                # would silently skip _load() altogether next time). Gated
                # on the parent_patched flag set above, rather than
                # comparing the parent's current attribute against
                # self._lazy_real - that comparison is wrong precisely when
                # it matters most: a failure *before* self._lazy_real gets
                # set (e.g. spec is None) leaves it as None, and a parent
                # with no such attribute at all also returns None from
                # .get(), so None is None would wrongly read as "yes, undo
                # it" and attach the stub to a parent it was never set on.
                if parent_patched:
                    parent = _sys.modules.get(parent_name)

                    if parent is not None:
                        setattr(parent, attr_name, self)

                self._lazy_real = None
                raise
            finally:
                self._lazy_owner = None

            return self._lazy_real

    def __dir__(self):
        # Unlike __repr__ below, this one *does* force the load. IPython's
        # tab-completion (`sr.mol.<TAB>`) goes through dir(), so leaving it
        # lazy here would mean nothing lazily-loaded ever tab-completes
        # before something else happens to have touched it - a real
        # regression from the old lazy_import package, which forced a load
        # on dir() too (_pythonize.py used to lean on exactly that as its
        # own force-load idiom). __repr__ stays lazy on purpose - that one
        # matters for debuggers/pytest failure reporting not wanting to
        # trigger a real import just to print a value.
        return dir(self._load())

    def __getattr__(self, attr):
        # only real attribute lookups (i.e. ones that fail against this
        # proxy's own __dict__) should trigger the load - this deliberately
        # does *not* override __repr__, so introspection that just wants to
        # print/log a value (e.g. a debugger, or pytest's failure
        # reporting) doesn't silently trigger a real import as a side
        # effect
        if attr.startswith("_lazy"):
            raise AttributeError(attr)

        # hasattr(module, "__path__") is the standard "is this a package?"
        # idiom, and before load it's answerable from the spec without
        # loading anything: a module whose spec has no
        # submodule_search_locations has no __path__ yet (module_from_spec()
        # itself only ever sets __path__ conditionally, for the same
        # reason __init__ above only sets it conditionally). Without this,
        # any code sweeping over sys.modules asking "which of these are
        # packages" - the same shape of scan the __file__ handling above
        # exists to protect against - would force-load every plain
        # (non-package) module stub in the tree at once.
        #
        # Gated on self._lazy_real is None (reading that field doesn't
        # itself force anything) because the spec only speaks for the
        # module *before* its body has run - a plain .py module can
        # still assign __path__ itself, making it act as a package (an
        # unusual but legitimate pattern). Once loaded, defer to the real
        # module like every other attribute does; only the pre-load case
        # is answerable without it.
        if (
            attr == "__path__"
            and self._lazy_real is None
            and self._lazy_spec is not None
            and self._lazy_spec.submodule_search_locations is None
        ):
            raise AttributeError(attr)

        real_module = self._load()

        try:
            return getattr(real_module, attr)
        except AttributeError:
            # 'attr' may be a submodule that lazy_module() pre-registered
            # under its own dotted name in sys.modules (see below), but
            # that was never actually set as an attribute of this (real,
            # now-loaded) parent module. That happens when the parent's
            # own body does `from . import sub` (rather than `from .sub
            # import Name`): CPython's own fromlist handling can bind
            # 'sub' straight from sys.modules without ever calling
            # getattr() on our stub, so nothing triggered *its* load or
            # did the usual parent-attribute fixup. Recover here instead,
            # on demand.
            full_name = f"{self._lazy_name}.{attr}"
            sub_module = _sys.modules.get(full_name)

            if sub_module is None:
                raise

            if isinstance(sub_module, _LazyModule):
                sub_module = sub_module._load()

            setattr(real_module, attr, sub_module)
            return sub_module

    def __repr__(self):
        if self._lazy_real is not None and self._lazy_owner is None:
            return repr(self._lazy_real)

        return f"<lazily-loaded module {self._lazy_name!r}>"


def _register_submodules(name, search_locations):
    """
    Pre-register a lazy stub for every *immediate* submodule of 'name',
    discovered purely from what's on disk under 'search_locations' - no
    code from any of these modules is executed to do this.

    This matters for correctness, not just convenience: if some external
    code (e.g. pickle, resolving a class's __module__) imports a
    submodule directly by its dotted path before anything has triggered
    the parent package's lazy load, Python's import machinery takes a
    fast path once it sees the parent already present in sys.modules -
    it fetches the parent's __path__ (which correctly triggers our
    load), but then goes on to build *and execute* the submodule itself
    regardless of whether that load already did so as a side effect,
    silently producing a second, distinct copy of that submodule (and
    therefore of any classes it defines). Pre-registering a stub for
    every immediate submodule under its own exact dotted name closes
    that gap, since the direct import then finds (and loads) that same
    stub instead of racing to rebuild it.

    Driving this from pkgutil.iter_modules()/PathFinder.find_spec()
    rather than a hardcoded module list means it stays correct for
    forks or downstream code that add new submodules, with no extra
    maintenance.

    Recursive: submodules of submodules get pre-registered too, all the
    way down, so a class living arbitrarily deep inside a lazily-loaded
    package is protected against the same bug. Note that pre-registering
    a stub here does *not* attach it as an attribute of its parent - it
    only reserves the name in sys.modules. The attribute fixup happens
    lazily instead, in _load() above, right after a package's own body
    finishes running. That ordering is deliberate: attaching children
    eagerly (before the parent's own __init__.py has ever run) would let
    an access that goes through the import system reach a child straight
    off the parent without ever triggering the parent's own
    initialisation - and some packages' submodules are only safe to touch
    in the order their __init__.py imports them in (e.g. because of shared
    state in a third-party library more than one of them happens to
    touch), an ordering guarantee a real, non-lazy import always preserves
    and which this needs to as well, for that same class of access. This
    doesn't (and can't) protect a lookup that bypasses the import system
    entirely, e.g. sys.modules["parent.child"] directly - that reaches the
    child's own stub without the parent ever being touched, exactly as it
    always could; not this module's problem to solve.
    """
    parent = _sys.modules.get(name)

    for _finder, sub_name, _is_pkg in _pkgutil.iter_modules(search_locations):
        full_name = f"{name}.{sub_name}"

        if full_name in _sys.modules:
            # Already registered (e.g. by some other path) - still worth
            # recording against the parent below, so _load()'s fixup loop
            # knows about it too.
            if isinstance(parent, _LazyModule):
                parent._lazy_children.append(full_name)
            continue

        # Use PathFinder directly with the already-known search_locations,
        # rather than importlib.util.find_spec(full_name) - that variant
        # resolves the parent package to read its __path__, which would
        # touch our still-lazy parent stub's __getattr__ and force it to
        # load prematurely, defeating the point of doing this lazily.
        try:
            spec = _importlib_machinery.PathFinder.find_spec(
                full_name, search_locations
            )
        except (ImportError, AttributeError):
            continue

        if spec is None:
            continue

        _sys.modules[full_name] = _LazyModule(full_name, search_locations, spec=spec)

        if isinstance(parent, _LazyModule):
            parent._lazy_children.append(full_name)

        if spec.submodule_search_locations:
            _register_submodules(full_name, spec.submodule_search_locations)


def lazy_module(name: str):
    """
    Return a proxy for the module called 'name' that defers the actual
    import until one of its attributes is accessed, or force_load() is
    called on it explicitly. This is a drop-in, GPLv3-free replacement for
    lazy_import.lazy_module().

    If 'name' is a package, every immediate submodule is also
    pre-registered as its own lazy stub - see _register_submodules() for
    why.

    Idempotent: calling this twice for the same name returns the *same*
    object both times, whatever state it's in (unloaded stub, mid-load, or
    fully real) - building a second, independent _LazyModule for a name
    already present in sys.modules would leave whoever's holding the first
    one with a second, unrelated loader for the same dotted name, which is
    exactly the class of bug this module exists to prevent.
    """
    existing = _sys.modules.get(name)

    if existing is not None:
        return existing

    # Find the spec *before* registering our stub - find_spec() consults
    # sys.modules first, and our stub has no real __spec__ of its own.
    try:
        spec = _importlib_util.find_spec(name)
    except (ImportError, AttributeError):
        spec = None

    module = _LazyModule(name, spec=spec)
    _sys.modules[name] = module

    if spec is not None and spec.submodule_search_locations:
        _register_submodules(name, spec.submodule_search_locations)

    return module


def is_lazy_module(module) -> bool:
    """Return whether 'module' is a not-yet-loaded lazy module proxy."""
    return isinstance(module, _LazyModule) and module._lazy_real is None


def force_load(module):
    """
    If 'module' is a lazy module proxy, force it to load now. This is
    a no-op if 'module' is already a real, fully-loaded module.
    """
    if isinstance(module, _LazyModule):
        module._load()
