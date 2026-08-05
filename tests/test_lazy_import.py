import pickle


def _make_and_dump(path):
    """Runs in a fresh worker process. Independently triggers sire's
    lazy-loading of sire.mol, then pickles a Molecule to disk."""
    import sire as sr

    mol = sr.mol.Molecule()
    with open(path, "wb") as f:
        pickle.dump(mol, f)


def _load_and_check(path):
    """Runs in another, separate fresh worker process. Independently
    triggers sire.mol's lazy load (before ever touching the pickle),
    then unpickles the Molecule and checks that its class is identical
    to (not just equal by name to) the one this process would construct
    itself."""
    import sire as sr

    Molecule = sr.mol.Molecule

    with open(path, "rb") as f:
        obj = pickle.load(f)

    return isinstance(obj, Molecule), type(obj) is Molecule


def test_lazy_import_pickle_across_processes(tmp_path):
    """
    A class from a lazily-loaded module, independently lazy-loaded in two
    separate worker processes, must be the same class object in both -
    isinstance()/type() checks on an object pickled in one process and
    unpickled in another must succeed.
    """
    from multiprocessing import get_context

    ctx = get_context("spawn")

    path = str(tmp_path / "molecule.pickle")

    # Dump in one, independent, fresh worker process.
    with ctx.Pool(1) as pool:
        pool.apply(_make_and_dump, (path,))

    # Load and check in a completely separate, fresh worker process.
    with ctx.Pool(1) as pool:
        isinstance_ok, type_is_ok = pool.apply(_load_and_check, (path,))

    assert isinstance_ok
    assert type_is_ok


def _check_direct_submodule_import():
    """
    Runs in a fresh process. Imports a submodule of a lazily-loaded
    package directly by its dotted path, *before* anything has triggered
    the parent package's own lazy load, and checks that the class
    obtained this way is identical to the one obtained via the parent's
    (now-loaded) attribute.
    """
    import sys

    import sire  # noqa: F401

    assert type(sys.modules["sire.mol"]).__name__ == "_LazyModule"

    # Note: sire.mol._element / Element, not sire.mol._trajectory, since
    # mol/__init__.py happens to also define its own unrelated function
    # called `_trajectory`, which shadows the submodule attribute
    # regardless of lazy loading (import a.b.c as x walks attributes,
    # not sys.modules, so a later same-named function always wins) -
    # that's an unrelated naming collision, not what this test targets.
    import sire.mol._element as leaf

    parent = sys.modules["sire.mol"]

    return leaf.Element is parent.Element


def test_lazy_import_direct_submodule_import():
    """
    Importing a submodule of a lazily-loaded package directly by its
    dotted path (e.g. `import sire.mol._x`, or pickle resolving a
    class's __module__), before the parent has ever been touched, must
    produce the same class object as accessing it via the parent's own
    (now-loaded) attribute - not a second, independent copy.
    """
    from multiprocessing import get_context

    ctx = get_context("spawn")

    with ctx.Pool(1) as pool:
        ok = pool.apply(_check_direct_submodule_import)

    assert ok


def _make_synthetic_package(root, pkg_name, mod_name, mod_body):
    """Write a minimal, importable package to disk: root/pkg_name/__init__.py
    (empty) and root/pkg_name/mod_name.py (mod_body). Returns str(root), for
    inserting onto sys.path in a worker process."""
    import os

    pkg_dir = os.path.join(root, pkg_name)
    os.makedirs(pkg_dir, exist_ok=True)

    with open(os.path.join(pkg_dir, "__init__.py"), "w") as f:
        f.write("")

    with open(os.path.join(pkg_dir, f"{mod_name}.py"), "w") as f:
        f.write(mod_body)

    return str(root)


def _thread_race_worker(root):
    """Runs in a fresh process. Registers a lazy stub for a submodule
    whose body sleeps before setting an attribute (widening the race
    window deliberately), then hammers that stub's first-touch load from
    many threads at once, synchronised to start together via a Barrier.
    Returns (errors, results, exec_count)."""
    import builtins
    import sys
    import threading

    sys.path.insert(0, root)

    from sire._lazy_import import lazy_module

    lazy_module("synth_race_pkg")
    stub = sys.modules["synth_race_pkg.slow_mod"]

    n_threads = 16
    barrier = threading.Barrier(n_threads)
    results = []
    errors = []
    results_lock = threading.Lock()

    def touch():
        barrier.wait()
        try:
            value = stub.VALUE
        except BaseException as exc:  # noqa: BLE001 - want any failure, not just AttributeError
            with results_lock:
                errors.append(repr(exc))
        else:
            with results_lock:
                results.append(value)

    threads = [threading.Thread(target=touch) for _ in range(n_threads)]

    for t in threads:
        t.start()
    for t in threads:
        t.join()

    exec_count = getattr(builtins, "_LAZY_TEST_RACE_EXEC_COUNT", 0)

    return errors, results, exec_count


def test_lazy_import_thread_race(tmp_path):
    """
    Many threads touching the same not-yet-loaded module for the first
    time, simultaneously, must all block until the load genuinely
    finishes: none may see a partially-initialised module or raise, and
    the module's own body must execute exactly once, not once per
    thread. The submodule's body sleeps before defining its one
    attribute, to widen the race window; builtins is used as a
    cross-thread (same-process) counter to confirm single execution.
    """
    body = (
        "import builtins\n"
        "import time\n"
        "\n"
        "builtins._LAZY_TEST_RACE_EXEC_COUNT = (\n"
        "    getattr(builtins, '_LAZY_TEST_RACE_EXEC_COUNT', 0) + 1\n"
        ")\n"
        "\n"
        "time.sleep(0.3)\n"
        "\n"
        "VALUE = 42\n"
    )
    root = _make_synthetic_package(str(tmp_path), "synth_race_pkg", "slow_mod", body)

    from multiprocessing import get_context

    ctx = get_context("spawn")

    with ctx.Pool(1) as pool:
        errors, results, exec_count = pool.apply(_thread_race_worker, (root,))

    assert errors == []
    assert results == [42] * 16
    assert exec_count == 1


def _failed_load_retry_worker(root):
    """Runs in a fresh process. Registers a lazy stub for a submodule
    whose body raises on its first execution and succeeds on the second,
    and checks that the failure is retryable rather than permanently
    masked."""
    import sys

    sys.path.insert(0, root)

    from sire._lazy_import import lazy_module

    lazy_module("synth_retry_pkg")
    stub = sys.modules["synth_retry_pkg.flaky_mod"]

    first_error = None

    try:
        stub.VALUE
    except RuntimeError as exc:
        first_error = str(exc)

    # A second attempt, on the *same* stub, should retry cleanly.
    second_value = stub.VALUE

    return first_error, second_value


def test_lazy_import_failed_load_is_retryable(tmp_path):
    """
    A module whose body raises on import must be retryable afterwards -
    the original exception should propagate (not a KeyError, and not a
    silently half-initialised module), and a later attempt must be able
    to succeed cleanly if the module's own code would now succeed.
    """
    body = (
        "import builtins\n"
        "\n"
        "builtins._LAZY_TEST_RETRY_ATTEMPTS = (\n"
        "    getattr(builtins, '_LAZY_TEST_RETRY_ATTEMPTS', 0) + 1\n"
        ")\n"
        "\n"
        "if builtins._LAZY_TEST_RETRY_ATTEMPTS == 1:\n"
        "    raise RuntimeError('simulated first-attempt failure')\n"
        "\n"
        "VALUE = 99\n"
    )
    root = _make_synthetic_package(str(tmp_path), "synth_retry_pkg", "flaky_mod", body)

    from multiprocessing import get_context

    ctx = get_context("spawn")

    with ctx.Pool(1) as pool:
        first_error, second_value = pool.apply(_failed_load_retry_worker, (root,))

    assert first_error == "simulated first-attempt failure"
    assert second_value == 99


def _lazy_module_idempotent_worker(root):
    """Runs in a fresh process. Calls lazy_module() twice for the same
    still-unloaded name and checks both calls return the identical
    object, then force-loads through one of the two references and
    confirms the other sees the same real module."""
    import sys

    sys.path.insert(0, root)

    from sire._lazy_import import force_load, lazy_module

    first = lazy_module("synth_idempotent_pkg")
    second = lazy_module("synth_idempotent_pkg")

    same_before_load = first is second

    force_load(first)

    # first/second (the stub) stay as they were - the real module replaces
    # the stub *in sys.modules*, not in variables already holding the stub
    # - so the meaningful check here is that sys.modules now holds exactly
    # the module first._load() produced, not some other, independent one.
    return (
        same_before_load,
        first is second,
        sys.modules["synth_idempotent_pkg"] is first._lazy_real,
    )


def test_lazy_module_is_idempotent(tmp_path):
    """
    lazy_module() called twice for the same still-unloaded name must
    return the identical object both times - not a second, independent
    stub - and that identity must still hold once the module is loaded.
    """
    root = _make_synthetic_package(
        str(tmp_path), "synth_idempotent_pkg", "leaf", "VALUE = 7\n"
    )

    from multiprocessing import get_context

    ctx = get_context("spawn")

    with ctx.Pool(1) as pool:
        same_before_load, same_after_load, matches_sys_modules = pool.apply(
            _lazy_module_idempotent_worker, (root,)
        )

    assert same_before_load
    assert same_after_load
    assert matches_sys_modules


def _find_spec_on_unloaded_stub_worker(root):
    """Runs in a fresh process. Registers a lazy stub, then calls
    importlib.util.find_spec() on it *before* touching it at all, and
    returns whatever that produced (a real ModuleSpec, or the repr of
    whatever exception it raised)."""
    import sys
    import importlib.util

    sys.path.insert(0, root)

    from sire._lazy_import import lazy_module

    lazy_module("synth_find_spec_pkg")

    assert type(sys.modules["synth_find_spec_pkg.leaf"]).__name__ == "_LazyModule"

    try:
        spec = importlib.util.find_spec("synth_find_spec_pkg.leaf")
    except BaseException as exc:  # noqa: BLE001 - want to see exactly what, if anything, was raised
        return None, repr(exc)

    return spec is not None, None


def test_find_spec_on_unloaded_stub(tmp_path):
    """
    importlib.util.find_spec() on a not-yet-loaded lazy module must
    return a real ModuleSpec, not raise - a stub's __spec__ (and
    __loader__, __package__) must be genuinely populated, not None.
    """
    root = _make_synthetic_package(
        str(tmp_path), "synth_find_spec_pkg", "leaf", "VALUE = 1\n"
    )

    from multiprocessing import get_context

    ctx = get_context("spawn")

    with ctx.Pool(1) as pool:
        got_spec, error = pool.apply(_find_spec_on_unloaded_stub_worker, (root,))

    assert error is None, f"find_spec() raised: {error}"
    assert got_spec


def _make_extension_package(root, pkg_name, so_module_name):
    """Write a minimal package (root/pkg_name/__init__.py, empty) whose
    one submodule is a real, copied stdlib C extension (single-phase
    init, e.g. _ctypes) rather than a .py file - named after the
    original so its PyInit_<name> symbol still matches. Returns
    str(root)."""
    import importlib
    import os
    import shutil

    ext_module = importlib.import_module(so_module_name)
    so_path = ext_module.__file__

    pkg_dir = os.path.join(root, pkg_name)
    os.makedirs(pkg_dir, exist_ok=True)

    with open(os.path.join(pkg_dir, "__init__.py"), "w") as f:
        f.write("")

    shutil.copy(so_path, os.path.join(pkg_dir, os.path.basename(so_path)))

    return str(root)


def _extension_module_not_eagerly_initialised_worker(root):
    """Runs in a fresh process. Registers a lazy stub for a package
    containing a real, single-phase-init C extension submodule, and
    checks that registering it doesn't dlopen/initialise it (the stub
    should have none of the extension's real attributes), while touching
    it afterwards does load the real thing into sys.modules."""
    import sys

    sys.path.insert(0, root)

    from sire._lazy_import import lazy_module

    lazy_module("synth_ext_pkg")

    stub = sys.modules["synth_ext_pkg._ctypes"]
    stub_is_lazy_module = type(stub).__name__ == "_LazyModule"

    # Anything not one of our own private '_lazy*' bookkeeping attributes
    # or a dunder here would mean the real extension's attributes were
    # already present at registration time, before anyone touched it.
    leaked_attrs = [
        a
        for a in stub.__dict__
        if not a.startswith("_lazy") and not (a.startswith("__") and a.endswith("__"))
    ]

    # Now actually touch it, and confirm the real module properly took
    # over in sys.modules (not just delegation through the stub).
    real_attr_count = len(dir(stub))
    real_module_in_sys_modules = type(sys.modules["synth_ext_pkg._ctypes"]).__name__

    return (
        stub_is_lazy_module,
        leaked_attrs,
        real_attr_count,
        real_module_in_sys_modules,
    )


def test_extension_module_not_eagerly_initialised(tmp_path):
    """
    A single-phase-init C extension (e.g. a compiled Boost.Python
    module) reachable inside a lazily-registered package must not be
    dlopened/initialised at registration time - the stub must carry none
    of the real extension's attributes until it's actually touched, and
    touching it must then produce the real module in sys.modules, with
    the real module's full attribute set. Uses a real, copied stdlib
    extension (_ctypes) rather than a synthetic pure-Python module,
    since this only shows up against something with real C-level
    initialisation side effects.
    """
    root = _make_extension_package(str(tmp_path), "synth_ext_pkg", "_ctypes")

    from multiprocessing import get_context

    ctx = get_context("spawn")

    with ctx.Pool(1) as pool:
        (
            stub_is_lazy_module,
            leaked_attrs,
            real_attr_count,
            real_module_type,
        ) = pool.apply(_extension_module_not_eagerly_initialised_worker, (root,))

    assert stub_is_lazy_module
    assert leaked_attrs == []
    assert real_attr_count > len(leaked_attrs)
    assert real_module_type == "module"


def _hasattr_path_reflects_load_state_worker(root):
    """Runs in a fresh process. Checks hasattr(stub, '__path__') on a
    not-yet-loaded module stub whose body assigns __path__ itself (an
    unusual but legitimate way for a plain module to act like a
    package): before load, that must be answerable without loading the
    module at all; after load, it must agree with the real module."""
    import sys

    sys.path.insert(0, root)

    from sire._lazy_import import lazy_module

    lazy_module("synth_path_pkg")
    stub = sys.modules["synth_path_pkg.leaf"]

    has_path_before = hasattr(stub, "__path__")
    still_lazy_before = (
        type(sys.modules["synth_path_pkg.leaf"]).__name__ == "_LazyModule"
    )

    stub.VALUE  # force the load

    has_path_after = hasattr(stub, "__path__")
    real_module_after = type(sys.modules["synth_path_pkg.leaf"]).__name__ == "module"

    return has_path_before, still_lazy_before, has_path_after, real_module_after


def test_hasattr_path_reflects_load_state(tmp_path):
    """
    hasattr(module, "__path__") - the standard "is this a package?" idiom
    - on a not-yet-loaded module stub whose spec has no
    submodule_search_locations must return False without loading the
    module. Once the module has loaded, the answer must switch to
    reflect what the real module actually has - including a plain
    module that assigns __path__ in its own body, which the spec alone
    can't predict.
    """
    root = _make_synthetic_package(
        str(tmp_path),
        "synth_path_pkg",
        "leaf",
        "__path__ = ['/nonexistent']\nVALUE = 1\n",
    )

    from multiprocessing import get_context

    ctx = get_context("spawn")

    with ctx.Pool(1) as pool:
        has_path_before, still_lazy_before, has_path_after, real_module_after = (
            pool.apply(_hasattr_path_reflects_load_state_worker, (root,))
        )

    assert has_path_before is False
    assert still_lazy_before
    assert has_path_after is True
    assert real_module_after


def _reentrant_before_real_worker(root):
    """Runs in a fresh process. Simulates the one window a same-thread
    reentrant import can hit before there's a module object to hand
    back (only reachable in practice for a C extension whose init
    routine itself imports something, which module_from_spec() can run
    before _load() ever gets to set self._lazy_real) by monkeypatching
    module_from_spec() to re-enter the stub's own _load() first."""
    import sys

    sys.path.insert(0, root)

    from sire import _lazy_import
    from sire._lazy_import import lazy_module

    lazy_module("synth_reentrant_pkg")
    stub = sys.modules["synth_reentrant_pkg.leaf"]

    original_module_from_spec = _lazy_import._importlib_util.module_from_spec
    caught = []

    def patched_module_from_spec(spec):
        if spec.name == "synth_reentrant_pkg.leaf" and not caught:
            try:
                stub._load()
            except ImportError as exc:
                caught.append(str(exc))
        return original_module_from_spec(spec)

    _lazy_import._importlib_util.module_from_spec = patched_module_from_spec
    try:
        value = stub.VALUE
    finally:
        _lazy_import._importlib_util.module_from_spec = original_module_from_spec

    return caught, value


def test_reentrant_import_before_module_exists_raises_clearly(tmp_path):
    """
    A same-thread reentrant _load() call that lands in the narrow window
    before the module object itself has been built must raise a clear,
    diagnosable error - not silently return None and let the caller hit
    an opaque AttributeError on it further down the line.
    """
    root = _make_synthetic_package(
        str(tmp_path), "synth_reentrant_pkg", "leaf", "VALUE = 1\n"
    )

    from multiprocessing import get_context

    ctx = get_context("spawn")

    with ctx.Pool(1) as pool:
        caught, value = pool.apply(_reentrant_before_real_worker, (root,))

    assert len(caught) == 1
    assert "circular import" in caught[0]
    assert "synth_reentrant_pkg.leaf" in caught[0]
    assert value == 1


def _make_star_import_package(root, pkg_name):
    """Write a package root/pkg_name whose __init__.py does
    `from .leaf import *`, where leaf.py defines a public name and has
    no __all__. Returns str(root)."""
    import os

    pkg_dir = os.path.join(root, pkg_name)
    os.makedirs(pkg_dir, exist_ok=True)

    with open(os.path.join(pkg_dir, "__init__.py"), "w") as f:
        f.write("from .leaf import *\n")

    with open(os.path.join(pkg_dir, "leaf.py"), "w") as f:
        f.write("VALUE = 123\n")

    return str(root)


def _star_import_without_all_worker(root):
    """Runs in a fresh process. Registers a lazy stub for a package
    whose __init__.py star-imports a submodule with no __all__, and
    checks that the star-imported name is reachable on the parent."""
    import sys

    sys.path.insert(0, root)

    from sire._lazy_import import lazy_module

    lazy_module("synth_star_pkg")
    stub = sys.modules["synth_star_pkg"]

    return stub.VALUE


def test_star_import_of_module_without_all(tmp_path):
    """
    `from .sibling import *`, where sibling has no __all__, makes every
    one of sibling's public names reachable on the package that did the
    star import, the same as it would for a module with __all__ or for
    a plain, non-lazy import.
    """
    root = _make_star_import_package(str(tmp_path), "synth_star_pkg")

    from multiprocessing import get_context

    ctx = get_context("spawn")

    with ctx.Pool(1) as pool:
        value = pool.apply(_star_import_without_all_worker, (root,))

    assert value == 123


def _force_load_then_reload_worker(root):
    """Runs in a fresh process. force_load()s a stub, then calls
    importlib.reload() on the resulting real module, and checks the
    module still works afterwards."""
    import importlib
    import sys

    sys.path.insert(0, root)

    from sire._lazy_import import force_load, lazy_module

    lazy_module("synth_reload_pkg")
    stub = sys.modules["synth_reload_pkg.leaf"]

    force_load(stub)
    real_before = sys.modules["synth_reload_pkg.leaf"]

    reloaded = importlib.reload(real_before)

    return (
        type(real_before).__name__,
        reloaded is sys.modules["synth_reload_pkg.leaf"],
        reloaded.VALUE,
    )


def test_reload_after_force_load(tmp_path):
    """
    force_load() followed by importlib.reload() on the resulting module
    works exactly as reload() does for any ordinary, non-lazy module.
    """
    root = _make_synthetic_package(
        str(tmp_path), "synth_reload_pkg", "leaf", "VALUE = 5\n"
    )

    from multiprocessing import get_context

    ctx = get_context("spawn")

    with ctx.Pool(1) as pool:
        real_type, same_module, value = pool.apply(
            _force_load_then_reload_worker, (root,)
        )

    assert real_type == "module"
    assert same_module
    assert value == 5
