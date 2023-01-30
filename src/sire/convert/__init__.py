__all__ = [
    "biosimspace_to_sire",
    "openmm_to_sire",
    "rdkit_to_sire",
    "sire_to_biosimspace",
    "sire_to_rdkit",
    "sire_to_openmm",
    "to",
    "to_biosimspace",
    "to_rdkit",
    "to_openmm",
    "to_sire",
]


def to(obj):
    """
    Convert the passed object from its current object format to a
    sire object format. Typically this will be converting
    from, e.g. a BioSimSpace, OpenMM or rdkit molecule to a sire molecule
    (or from a list of molecules to a SelectorMol).
    """
    return to_sire(obj)


def _convert(obj, converter):
    """
    Internal function that converts 'obj' using the passed converter
    function.
    """
    if type(obj) is list:
        converted = []

        for o in obj:
            converted.append(converter(o))

        try:
            # Try to return this as a molecules collection
            from ..mol import SelectorMol

            return SelectorMol(converted)
        except Exception:
            # this failed, so we can only assume that the
            # objects are not all molecules. Just return the list
            pass

        return converted

    is_selector = False

    try:
        is_selector = obj.is_selector()
    except Exception:
        pass

    if is_selector:
        converted = []

        for o in obj:
            converted.append(converter(o))

        try:
            # Try to return this as a molecules collection
            from ..mol import SelectorMol

            return SelectorMol(converted)
        except Exception:
            # this failed, so we can only assume that the
            # objects are not all molecules. Just return the list
            pass

        return converted

    # assume this is a single molecule
    return converter(obj)


def to_sire(obj):
    """
    Convert the passed object from its current object format to a
    sire object format. Typically this will be converting
    from, e.g. an openmm or rdkit molecule to a sire molecule
    (or from a list of molecules to a SelectorMol).
    """

    if "sire" in str(type(obj)):
        # already a sire object?
        return obj

    def _to_sire(o):
        t = str(type(o))

        if "sire" in t:
            # this is already a sire molecule
            return o.molecule()

        elif "BioSimSpace" in t:
            # this is a BioSimSpace molecule?
            return biosimspace_to_sire(o)

        elif "rdkit" in t:
            return rdkit_to_sire(o)

        elif "openmm" in t:
            return openmm_to_sire(o)

        else:
            raise TypeError(
                f"Cannot convert '{o}' as it is of unrecognised type {type(o)}"
            )

    return _convert(obj, converter=_to_sire)


def to_biosimspace(obj):
    """
    Convert the passed object from its current object format to a
    BioSimSpace object format.
    """
    return sire_to_biosimspace(to_sire(obj))


def to_rdkit(obj):
    """
    Convert the passed object from its current object format to a
    rdkit object format.
    """
    return sire_to_rdkit(to_sire(obj))


def to_openmm(obj):
    """
    Convert the passed object from its current object format to an
    openmm object format.
    """
    return sire_to_openmm(to_sire(obj))


def biosimspace_to_sire(obj):
    """
    Convert the passed BioSimSpace object (either a Molecule or list
    of Molecules) to a sire equivalent
    """
    # will eventually support System too...
    if type(obj) is list:
        converted = []
        for o in obj:
            converted.append(biosimspace_to_sire(o))

        return converted

    if not hasattr(obj, "_sire_object"):
        raise TypeError(
            f"The object {obj} of type {type(obj)} does not look like a "
            "supported BioSimSpace object that can be converted to a "
            "sire object. Supported objects are Molecule."
        )

    obj = obj._sire_object.molecules()

    if obj.num_molecules() == 1:
        return obj[0]
    elif obj.num_molecules() == 0:
        return None
    else:
        return obj


def sire_to_biosimspace(obj):
    """
    Convert the passed sire object (either a molecule or list
    of molecules) to a BioSimSpace equivalent
    """
    # will eventually support System too...
    obj = obj.molecules()

    if obj.num_molecules() == 1:
        obj = obj[0]
    elif obj.num_molecules() == 0:
        return None

    try:
        import BioSimSpace as BSS
    except (ImportError, ModuleNotFoundError):
        raise ModuleNotFoundError(
            "BioSimSpace is not available. Please install via "
            "'mamba install -c openbiosim biosimspace'"
        )
    except Exception as e:
        raise ImportError(
            "There was an error importing BioSimSpace. This can occur "
            "if you import BioSimSpace after sire. Please import "
            "BioSimSpace first and try again. For info, the error "
            f"was {e}"
        )

    # Now convert this to a BioSimSpace object
    def _to_biosimspace(o):
        """
        Convert the passed sire Molecule to a BioSimSpace Molecule
        """
        return BSS._SireWrappers.Molecule(o)

    return _convert(obj, converter=_to_biosimspace)


def openmm_to_sire(obj):
    """
    Convert the passed OpenMM object (either a molecule or
    list of molecules) to the sire equivalent
    """
    # will eventually support System too...
    if type(obj) is list:
        converted = []
        for o in obj:
            converted.append(openmm_to_sire(o))

        return converted

    try:
        from ..legacy.Convert._sire_openmm import from_openmm as _from_openmm
    except Exception:
        raise ModuleNotFoundError(
            "openmm is not available. Please install via "
            "'mamba install -c conda-forge openmm'"
        )

    obj = _from_openmm(obj)

    obj = obj.molecules()

    if obj.num_molecules() == 1:
        return obj[0]
    elif obj.num_molecules() == 0:
        return None
    else:
        return obj


def sire_to_openmm(obj):
    """
    Convert the passed sire object (either a molecule or list
    of molecules) to an OpenMM equivalent
    """
    # will eventually support System too...
    obj = obj.molecules()

    if obj.num_molecules() == 1:
        obj = obj[0]
    elif obj.num_molecules() == 0:
        return None

    try:
        from ..legacy.Convert._sire_openmm import to_openmm as _to_openmm
    except Exception:
        raise ModuleNotFoundError(
            "openmm is not available. Please install via "
            "'mamba install -c conda-forge openmm'"
        )

    return _convert(obj, converter=_to_openmm)


def rdkit_to_sire(obj):
    """
    Convert the passed rdkit object (either a molecule or
    list of molecules) to the sire equivalent
    """
    if type(obj) is list:
        converted = []
        for o in obj:
            converted.append(rdkit_to_sire(o))

        return converted

    try:
        from ..legacy.Convert._sire_rdkit import from_rdkit as _from_rdkit
    except Exception:
        raise ModuleNotFoundError(
            "rdkit is not available. Please install via "
            "'mamba install -c conda-forge rdkit'"
        )

    obj = _from_rdkit(obj)

    obj = obj.molecules()

    if obj.num_molecules() == 1:
        return obj[0]
    elif obj.num_molecules() == 0:
        return None
    else:
        return obj


def sire_to_rdkit(obj):
    """
    Convert the passed sire object (either a molecule or list
    of molecules) to a rdkit equivalent
    """
    obj = obj.molecules()

    if obj.num_molecules() == 1:
        obj = obj[0]
    elif obj.num_molecules() == 0:
        return None

    # Now convert to an rdkit object
    try:
        from ..legacy.Convert._sire_rdkit import to_rdkit as _to_rdkit
    except Exception:
        raise ModuleNotFoundError(
            "rdkit is not available. Please install via "
            "'mamba install -c conda-forge rdkit'"
        )

    return _convert(obj, converter=_to_rdkit)
