__all__ = ["try_import", "try_import_from"]

_module_to_package = {}


def _find_conda():
    import os
    import sys

    conda_base = os.path.abspath(os.path.dirname(sys.executable))

    if os.path.basename(conda_base) == "bin":
        conda_base = os.path.dirname(conda_base)

    conda = None

    if "CONDA_EXE" in os.environ:
        conda = os.environ["CONDA_EXE"]
    else:
        conda = None

    if "CONDA_DEFAULT_ENV" in os.environ:
        conda_env = os.environ["CONDA_DEFAULT_ENV"]
    else:
        conda_env = None

    if os.path.exists(os.path.join(conda_base, "python.exe")):
        # Windows
        conda_bin = os.path.join(conda_base, "Library", "bin")

        if conda is None:
            conda = os.path.join(conda_base, "Scripts", "conda.exe")
    elif os.path.exists(os.path.join(conda_base, "bin", "python")):
        # MacOS and Linux
        conda_bin = os.path.join(conda_base, "bin")

        if conda is None:
            conda = os.path.join(conda_bin, "conda")
    else:
        print(
            "Cannot find a 'python' binary in directory '%s'. "
            "Are you running this script using the python executable "
            "from a valid miniconda or anaconda installation?" % conda_base
        )
        return None

    if conda.endswith(".exe"):
        m = os.path.join(os.path.dirname(conda), "mamba.exe")
    else:
        m = os.path.join(os.path.dirname(conda), "mamba")

    if os.path.exists(m):
        return m
    else:
        return conda


def _install_package(name, package_registry):
    """
    Internal function used to install the module
    called 'name', using the passed 'package_registry'
    to find the package that contains the package
    that contains this module
    """

    conda = _find_conda()

    # ensure that we have the root package name
    try:
        package = name.split(".")[0]
    except Exception:
        package = name

    if package in package_registry:
        package = package_registry[name]

    import os

    try:
        print(
            "\nTrying to install %s from package %s using %s...\n"
            % (name, package, conda)
        )
        ok = os.system("%s install %s -y" % (conda, package))

        if ok == 0:
            # installed ok
            return
    except Exception:
        pass

    print(
        "\nWARNING: Unable to install '%s' from package '%s'\n"
        % (name, package)
    )


def try_import(name, package_registry=_module_to_package):
    """Try to import the module called 'name', returning
    the loaded module as an argument. If the module
    is not available, then it looks up the name of
    the package to install using "package_registry"
    (or if this is not available, using just the name
    of the module). This will then be installed using
    "conda", then "pip" then "easy_install" (first one
    that works will return).

    For example, use this via

    sys = try_import("sys")
    mdtraj = try_import("mdtraj")

    Note that you can also rename modules, e.g. by using

    md = try_import("mdtraj")

    Note that you should use try_import_from if you
    want only specific symbols, e.g.

    (argv, stdout) = try_import_from("sys", ["argv","stdout"])
    """

    try:
        mod = __import__(name)
        return mod
    except Exception:
        pass

    if not (package_registry is None):
        _install_package(name, package_registry)
        return try_import(name, package_registry=None)

    raise ImportError("Failed to install module %s" % name)


def try_import_from(name, fromlist, package_registry=_module_to_package):
    """Try to import from the module called 'name' the passed symbol
    (or list of symbols) contained in 'fromlist', returning
    the symbol (or list of symbols).

    If the module cannot be loaded, then the package containing
    the module is looked up in 'module_to_package' (or just guessed
    from the name if it does not exist in 'module_to_package'.
    An attempt is made to load the package, using first conda,
    then pip, then easy_install.

    Example usage:

    mol = try_import_from("sire", "mol")
    (argv,stdout = try_import_from("sys", ["argv", "stdout"])
    ut = try_import_from("mdtraj", "utils")
    """

    if isinstance(fromlist, str):
        # we are importing only a single module - put
        # this string into a list for the user
        fromlist = [fromlist]

    try:
        nsyms = len(fromlist)
    except Exception:
        return try_import(name, package_registry)

    if nsyms == 0:
        # just import the entire module
        return try_import(name, package_registry)

    is_loaded = False

    try:
        mod = __import__(name, globals(), locals(), fromlist)
        is_loaded = True
    except Exception:
        is_loaded = False

    if not is_loaded:
        if not (package_registry is None):
            _install_package(name, package_registry)
            return try_import_from(name, fromlist, package_registry=None)
        else:
            raise ImportError("Failed to install module '%s'" % name)

    if nsyms == 1:
        try:
            return getattr(mod, fromlist[0])
        except Exception:
            raise ImportError(
                "Cannot find the symbol '%s' in module '%s'"
                % (fromlist[0], name)
            )
    else:
        ret = []
        missing_symbols = []

        for sym in fromlist:
            try:
                ret.append(getattr(mod, sym))
            except Exception:
                missing_symbols.append(sym)

        if len(missing_symbols) > 0:
            raise ImportError(
                "Cannot find the following symbols in module '%s' : [ %s ]"
                % (name, ", ".join(missing_symbols))
            )

        return ret
