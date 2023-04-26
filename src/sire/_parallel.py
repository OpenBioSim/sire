__all__ = [
    "get_max_num_threads",
    "set_default_num_threads",
    "set_max_num_threads",
]


def get_max_num_threads():
    """
    Return the maximum number of C++ threads that will be allowed
    for computation. Note that this is the maximum number of compute
    threads managed via Intel's Threading Building blocks. Sire
    additionally uses a small number of lightweight background
    threads for management, which are unaffected by this setting.
    """
    from .legacy import Base as _Base

    return _Base.get_max_num_threads()


def set_max_num_threads(n: int):
    """
    Set the maximum number of C++ threads that will be allowed
    for computation. Note that this is the maximum number of compute
    threads managed via Intel's Threading Building blocks. Sire
    additionally uses a small number of lightweight background
    threads for management, which are unaffected by this setting.
    """
    from .legacy import Base as _Base

    _Base.set_max_num_threads(n)


def set_default_num_threads():
    """
    Set a reasonable default number of C++ threads that will be
    allowed for computation. Note that this is the maximum number of compute
    threads managed via Intel's Threading Building blocks. Sire
    additionally uses a small number of lightweight background
    threads for management, which are unaffected by this setting.
    """
    import os as _os
    from .legacy import Base as _Base

    try:
        n = int(_os.environ["TBB_NUM_THREADS"])
    except Exception:
        try:
            n = int(_os.environ["OMP_NUM_THREADS"])
        except Exception:
            n = None

    if n is None:
        _Base.set_default_num_threads()
    else:
        if n < 1:
            n = 1

        _Base.set_max_num_threads(n)


# start off setting the default number - this will take into
# account any environment variables that have been set
set_default_num_threads()
