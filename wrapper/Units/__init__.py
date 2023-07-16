from ._Units import *

from .. import Base as _Base


_have_set_internal_units = False


def _set_internal_units():
    global _have_set_internal_units

    if _have_set_internal_units:
        return

    _have_set_internal_units = True

    try:
        GeneralUnit.clearDefaults()
    except AttributeError:
        # we have activated the new api, which will
        # handle setting internal units for us
        return

    for unit in [
        "g",
        "kcal",
        "Å",
        "ps",
        "|e|",
        "°",
        "K",
        "kcal mol-1",
        "kcal mol-1 Å-1",
        "kcal mol-1 Å-2",
        "kcal mol-1 Å-3",
        "kcal mol-1 °-1",
        "kcal mol-1 °-2",
        "kcal mol-1 K-1",
        "kcal mol-1 K-2",
        "kcal s",
        "kcal s-1",
        "kcal mol-1 s",
        "kcal mol-1 s-1",
        "kcal Å-1",
    ]:
        u = GeneralUnit(unit)
        u.setAsDefault(unit)


_set_internal_units()
