from . import config

import os as _os
import warnings as _warnings

from ._pythonize import use_mixed_api, use_new_api, use_old_api
from ._load import (
    load,
    save,
    save_to_string,
    expand,
    tutorial_url,
    load_test_files,
    supported_formats,
    smiles,
    smarts,
)

from ._measure import measure, minimum_distance
from ._colname import colname, colnames
from ._match import match_atoms
from ._parallel import (
    get_max_num_threads,
    set_max_num_threads,
    set_default_num_threads,
)

__all__ = [
    "atomid",
    "chainid",
    "colname",
    "colnames",
    "expand",
    "get_max_num_threads",
    "load",
    "load_test_files",
    "match_atoms",
    "minimum_distance",
    "measure",
    "molid",
    "save",
    "save_to_string",
    "segid",
    "set_default_num_threads",
    "set_max_num_threads",
    "smiles",
    "smarts",
    "sqrt",
    "supported_formats",
    "tutorial_url",
    "u",
    "use_mixed_api",
    "use_new_api",
    "use_old_api",
    "v",
]


# filter out annoying double-wrapped warnings
_warnings.filterwarnings("ignore", "to-Python converter for")


def _fix_openmm_path():
    """We need to fix the OpenMM path on Windows, because the DLL
    is not where we would expect it to be
    """
    import sys

    if sys.platform != "win32":
        return

    import os
    import glob

    condadir = os.path.dirname(sys.executable)

    # The DLL is put in libdir, but other DLLs are put in bindir
    libdir = os.path.join(condadir, "Library", "lib")
    bindir = os.path.join(condadir, "Library", "bin")
    openmm_files = glob.glob(os.path.join(libdir, "OpenMM*.dll"))

    need_path = False

    for file in openmm_files:
        binfile = os.path.join(bindir, os.path.basename(file))

        if not os.path.exists(binfile):
            # We do need to add bindir to the PATH
            need_path = True
            break

    if need_path:
        # importing openmm should add this path
        try:
            import openmm  # noqa: F401 (imported but unused)
        except Exception:
            print("OpenMM import failed!")
            # Copy the files
            try:
                for file in openmm_files:
                    binfile = os.path.join(bindir, os.path.basename(file))

                    if not os.path.exists(binfile):
                        # copy the file into the right place
                        import shutil

                        shutil.copy(file, binfile)
            except Exception:
                print(
                    "Could not resolve OpenMM library location. "
                    "This may cause issues later."
                )

        # also add it manually here
        try:
            os.add_dll_directory(libdir)
        except Exception:
            # os.add_dll_directory is not available on old Python (3.7)
            pass


_fix_openmm_path()


def sqrt(x):
    """Return the square root of the passed value"""
    if hasattr(x, "sqrt"):
        return x.sqrt()
    else:
        import math

        return math.sqrt(x)


def u(unit):
    """
    Return a sire unit created from the passed expression. If this is a
    sire.units.GeneralUnit then it will be returned. If this is a string,
    then it will be parseed and returned as a sire.units.GeneralUnit.
    """
    from .units import GeneralUnit

    error = None

    try:
        return GeneralUnit(unit)
    except Exception as e:
        error = e

    # is this a different unit model?

    # Try BioSimSpace
    if str(type(unit)).find("BioSimSpace") != -1:
        try:
            return GeneralUnit(unit._sire_unit)
        except Exception:
            pass

        return GeneralUnit(f"{unit}")

    # Try Pint
    if str(type(unit)).find("pint") != -1 and hasattr(type(unit), "magnitude"):
        # this is a pint unit - convert to a string (using the long-default
        # format, as sire should be able to read and understand this)
        # (we can't use short default as we use 'A' to mean angstrom, not amp)
        return GeneralUnit(unit.magnitude, f"{unit.units}")

    # Just try representing this as a string and see what happens
    try:
        return GeneralUnit(f"{unit}")
    except Exception as e:
        raise TypeError(
            f"Could not convert {unit} to a sire unit. The original error "
            f"that was raised was: {error}. The error raised when parsing "
            f"the string version was {e}"
        )


def v(x, y=None, z=None, units=None):
    """
    Return a sire vector from the passed expression. If this is a set of
    numbers or lengths (or a combination) then a sire.maths.Vector will
    be returned. If this is a value with velocity or force units then
    a Velocity3D or Force3D will be returned. If there is no
    vector type for data of this value then a simple python vector object
    will be returned.

    Args:
        x: The x-value, or something containing 3 values
        y: The y-value (cannot be specified if x has more than 1 value)
        z: The z-value (cannot be specified if x has more than 1 value)
        units: The units of the passed values (optional - will be guessed
               if not specified). You should not pass this if x, y or z
               already have values.
    """
    if type(x) is dict:
        if y is not None or z is not None:
            raise ValueError("You cannot specify y or z values when passing a dict.")

        y = x["x"]
        z = x["z"]
        x = x["x"]

    elif (not isinstance(x, str)) and hasattr(x, "__len__"):
        if len(x) != 3:
            raise ValueError(
                "The passed list or tuple must have three elements to be "
                f"converted to a Vector - the value '{x}' is not valid."
            )

        if y is not None or z is not None:
            raise ValueError(
                "You cannot specify y or z values when passing a list or tuple."
            )

        x, y, z = (x[0], x[1], x[2])

    else:
        if y is None:
            y = 0

        if z is None:
            z = 0

    if units is not None:
        u_units = u(units)

        if u_units.temperature() == 0:
            x *= u_units
            y *= u_units
            z *= u_units
        else:
            from .units import kelvin

            if u_units.has_same_units(kelvin):
                x = u(f"{x} {units}")
                y = u(f"{y} {units}")
                z = u(f"{z} {units}")
            else:
                raise ValueError(
                    "You can't specify units that include temperature, "
                    "as this can't be mulitplied easily."
                )

    x = u(x)
    y = u(y)
    z = u(z)

    from .maths import Vector

    # find the units of the passed values
    if x.is_zero() and y.is_zero() and z.is_zero():
        return Vector(0)

    units = None

    if not x.is_dimensionless():
        units = x.units()

    if not y.is_dimensionless():
        if units is None:
            units = y.units()
        elif not units.has_same_units(y):
            raise ValueError(
                "The passed y value has units that are aren't compatible with x. "
                f"{x} versus {y}"
            )

    if not z.is_dimensionless():
        if units is None:
            units = z.units()
        elif not units.has_same_units(z):
            raise ValueError(
                "The passed z value has units that are aren't compatible with x or y. "
                f"{x} versus {y} versus {z}"
            )

    if units is None:
        # all dimensionless - will be a simple vector
        return Vector(x.value(), y.value(), z.value())
    else:
        if x.is_dimensionless():
            x = x * units

        if y.is_dimensionless():
            y = y * units

        if z.is_dimensionless():
            z = z * units

    # we have units - need to create a vector with the right type
    from .units import angstrom

    if units.has_same_units(angstrom):
        return Vector(x, y, z)

    from .units import picosecond

    if units.has_same_units(angstrom / picosecond):
        from .legacy.Mol import Velocity3D

        return Velocity3D(x, y, z)

    from .units import newton

    if units.has_same_units(newton):
        from .legacy.Mol import Force3D

        return Force3D(x, y, z)

    # no vector type for this - just return a simple vector
    return (x, y, z)


def molid(
    num: int = None,
    name: str = None,
    idx: int = None,
    case_sensitive: bool = True,
):
    """Construct an identifer for a Molecule from the passed
       name, number and index.

    Args:
        name (str, optional): The molecule name. Defaults to None.
        num (int, optional): The molecule number. Defaults to None.
        idx (int, optional): The molecule index. Defaults to None.
        case_sensitive (bool): Whether or not the name is case sensitive

    Returns:
        MolID : The returned molecule identifier
    """
    ID = None

    if type(num) is str:
        # used in unnamed argument mode
        if name is None:
            name = num
            num = None
        elif type(name) is int:
            num, name = (name, num)
        else:
            raise TypeError("The number cannot be a string.")

    from .mol import MolName, MolNum, MolIdx

    if name is not None:
        if case_sensitive:
            from .id import CaseSensitive

            cs = CaseSensitive
        else:
            from .id import CaseInsensitive

            cs = CaseInsensitive

        ID = MolName(name, cs)

    if num is not None:
        if ID is None:
            ID = MolNum(num)
        else:
            ID = ID + MolNum(num)

    if idx is not None:
        if ID is None:
            ID = MolIdx(idx)
        else:
            ID = ID + MolIdx(idx)

    if ID is None:
        return MolIdx()
    else:
        return ID


def atomid(num: int = None, name: str = None, idx: int = None, case_sensitive=True):
    """Construct an identifer for an Atom from the passed
       name, number and index.

    Args:
        name (str, optional): The atom name. Defaults to None.
        num (int, optional): The atom number. Defaults to None.
        idx (int, optional): The atom index. Defaults to None.
        case_sensitive (bool): Whether the name is case sensitive or not

    Returns:
        AtomID : The returned atom identifier
    """
    ID = None

    if type(num) is str:
        # used in unnamed argument mode
        if name is None:
            name = num
            num = None
        elif type(name) is int:
            num, name = (name, num)
        else:
            raise TypeError("The number cannot be a string.")

    from .mol import AtomName, AtomNum, AtomIdx

    if name is not None:
        if case_sensitive:
            from .id import CaseSensitive

            cs = CaseSensitive
        else:
            from .id import CaseInsensitive

            cs = CaseInsensitive

        ID = AtomName(name, cs)

    if num is not None:
        if ID is None:
            ID = AtomNum(num)
        else:
            ID = ID + AtomNum(num)

    if idx is not None:
        if ID is None:
            ID = AtomIdx(idx)
        else:
            ID = ID + AtomIdx(idx)

    if ID is None:
        return AtomIdx()
    else:
        return ID


def resid(num: int = None, name: str = None, idx: int = None, case_sensitive=True):
    """Construct an identifer for a Residue from the passed
       name, number and index.

    Args:
        name (str, optional): The residue name. Defaults to None.
        number (int, optional): The residue number. Defaults to None.
        index (int, optional): The residue index. Defaults to None.
        case_sensitive (bool): Whether or not the name is case sensitive

    Returns:
        ResID : The returned atom identifier
    """
    ID = None

    if type(num) is str:
        # used in unnamed argument mode
        if name is None:
            name = num
            num = None
        elif type(name) is int:
            num, name = (name, num)
        else:
            raise TypeError("The number cannot be a string.")

    from .mol import ResName, ResNum, ResIdx

    if name is not None:
        if case_sensitive:
            from .id import CaseSensitive

            cs = CaseSensitive
        else:
            from .id import CaseInsensitive

            cs = CaseInsensitive

        ID = ResName(name, cs)

    if num is not None:
        if ID is None:
            ID = ResNum(num)
        else:
            ID = ID + ResNum(num)

    if idx is not None:
        if ID is None:
            ID = ResIdx(idx)
        else:
            ID = ID + ResIdx(idx)

    if ID is None:
        return ResIdx()
    else:
        return ID


def chainid(idx: int = None, name: str = None, case_sensitive: bool = True):
    """Construct an identifer for a Chain from the passed
       name and index.

    Args:
        name (str, optional): The chain name. Defaults to None.
        index (int, optional): The chain index. Defaults to None.
        case_sensitive (bool): Whether or not the name is case sensitive

    Returns:
        ChainID : The returned chain identifier
    """
    ID = None

    if type(idx) is str:
        # used in unnamed argument mode
        if name is None:
            name = idx
            idx = None
        elif type(name) is int:
            idx, name = (name, idx)
        else:
            raise TypeError("The index cannot be a string.")

    from .mol import ChainName, ChainIdx

    if name is not None:
        if case_sensitive:
            from .id import CaseSensitive

            cs = CaseSensitive
        else:
            from .id import CaseInsensitive

            cs = CaseInsensitive

        ID = ChainName(name, cs)

    if idx is not None:
        if ID is None:
            ID = ChainIdx(idx)
        else:
            ID = ID + ChainIdx(idx)

    if ID is None:
        return ChainIdx()
    else:
        return ID


def segid(idx: int = None, name: str = None, case_sensitive: bool = True):
    """Construct an identifer for a Segment from the passed
       name and index.

    Args:
        name (str, optional): The segment name. Defaults to None.
        index (int, optional): The segment index. Defaults to None.
        case_sensitive (bool): Whether or not the name is case sensitive

    Returns:
        SegID : The returned chain identifier
    """
    ID = None

    if type(idx) is str:
        # used in unnamed argument mode
        if name is None:
            name = idx
            idx = None
        elif type(name) is int:
            idx, name = (name, idx)
        else:
            raise TypeError("The index cannot be a string.")

    from .mol import SegName, SegIdx

    if name is not None:
        if case_sensitive:
            from .id import CaseSensitive

            cs = CaseSensitive
        else:
            from .id import CaseInsensitive

            cs = CaseInsensitive

        ID = SegName(name, cs)

    if idx is not None:
        if ID is None:
            ID = SegIdx(idx)
        else:
            ID = ID + SegIdx(idx)

    if ID is None:
        return SegIdx()
    else:
        return ID


def bondid(atom0, atom1, case_sensitive: bool = True):
    """Construct an identifier for a Bond from the passed
       identifiers for the two atoms.

       The atom identifiers can be:

       * integers - in this case they are treated as Atom indexes
       * strings - in this case they are treated as Atom names
       * AtomIDs - these are AtomIDs created via, e.g. the atomid function.

    Args:
        atom0 (int, str, AtomID): The identifier for the first atom.
        atom1 (int, str, AtomID): The identifier for the second atom.
        case_sensitive: Whether or not the name is case sensitive

    Returns:
        BondID: The returned bond identifier
    """

    def _convert(id):
        if type(id) is int:
            from .mol import AtomIdx

            return AtomIdx(id)
        elif type(id) is str:
            if case_sensitive:
                from .id import CaseSensitive

                cs = CaseSensitive
            else:
                from .id import CaseInsensitive

                cs = CaseInsensitive

            from .mol import AtomName

            return AtomName(id, cs)
        else:
            from .mol import AtomID

            if AtomID in type(id).mro():
                return id
            else:
                return atomid(id)

    from .mol import BondID

    return BondID(_convert(atom0), _convert(atom1))


def angleid(atom0, atom1, atom2, case_sensitive: bool = True):
    """Construct an identifier for a Angle from the passed
       identifiers for the three atoms.

       The atom identifiers can be:

       * integers - in this case they are treated as Atom indexes
       * strings - in this case they are treated as Atom names
       * AtomIDs - these are AtomIDs created via, e.g. the atomid function.

    Args:
        atom0 (int, str, AtomID): The identifier for the first atom.
        atom1 (int, str, AtomID): The identifier for the second atom.
        atom2 (int, str, AtomID): The identifier for the third atom.
        case_sensitive: Whether or not the name is case sensitive

    Returns:
        AngleID: The returned angle identifier
    """

    def _convert(id):
        if type(id) is int:
            from .mol import AtomIdx

            return AtomIdx(id)
        elif type(id) is str:
            if case_sensitive:
                from .id import CaseSensitive

                cs = CaseSensitive
            else:
                from .id import CaseInsensitive

                cs = CaseInsensitive

            from .mol import AtomName

            return AtomName(id, cs)
        else:
            from .mol import AtomID

            if AtomID in type(id).mro():
                return id
            else:
                return atomid(id)

    from .mol import AngleID

    return AngleID(_convert(atom0), _convert(atom1), _convert(atom2))


def dihedralid(atom0, atom1, atom2, atom3, case_sensitive: bool = True):
    """Construct an identifier for a Dihedral from the passed
       identifiers for the four atoms.

       The atom identifiers can be:

       * integers - in this case they are treated as Atom indexes
       * strings - in this case they are treated as Atom names
       * AtomIDs - these are AtomIDs created via, e.g. the atomid function.

    Args:
        atom0 (int, str, AtomID): The identifier for the first atom.
        atom1 (int, str, AtomID): The identifier for the second atom.
        atom2 (int, str, AtomID): The identifier for the third atom.
        atom3 (int, str, AtomID): The identifier for the fourth atom.
        case_sensitive: Whether or not the name is case sensitive

    Returns:
        DihedralID: The returned dihedral identifier
    """

    def _convert(id):
        if type(id) is int:
            from .mol import AtomIdx

            return AtomIdx(id)
        elif type(id) is str:
            if case_sensitive:
                from .id import CaseSensitive

                cs = CaseSensitive
            else:
                from .id import CaseInsensitive

                cs = CaseInsensitive

            from .mol import AtomName

            return AtomName(id, cs)
        else:
            from .mol import AtomID

            if AtomID in type(id).mro():
                return id
            else:
                return atomid(id)

    from .mol import DihedralID

    return DihedralID(
        _convert(atom0), _convert(atom1), _convert(atom2), _convert(atom3)
    )


def improperid(atom0, atom1, atom2, atom3, case_sensitive: bool = True):
    """Construct an identifier for an Improper from the passed
       identifiers for the four atoms.

       The atom identifiers can be:

       * integers - in this case they are treated as Atom indexes
       * strings - in this case they are treated as Atom names
       * AtomIDs - these are AtomIDs created via, e.g. the atomid function.

    Args:
        atom0 (int, str, AtomID): The identifier for the first atom.
        atom1 (int, str, AtomID): The identifier for the second atom.
        atom2 (int, str, AtomID): The identifier for the third atom.
        atom3 (int, str, AtomID): The identifier for the fourth atom.
        case_sensitive: Whether or not the name is case sensitive

    Returns:
        ImproperID: The returned improper identifier
    """

    def _convert(id):
        if type(id) is int:
            from .mol import AtomIdx

            return AtomIdx(id)
        elif type(id) is str:
            if case_sensitive:
                from .id import CaseSensitive

                cs = CaseSensitive
            else:
                from .id import CaseInsensitive

                cs = CaseInsensitive

            from .mol import AtomName

            return AtomName(id, cs)
        else:
            from .mol import AtomID

            if AtomID in type(id).mro():
                return id
            else:
                return atomid(id)

    from .mol import ImproperID

    return ImproperID(
        _convert(atom0), _convert(atom1), _convert(atom2), _convert(atom3)
    )


__version__ = config.__version__

__branch__ = config.sire_repository_branch
__repository__ = config.sire_repository_url
__revisionid__ = config.sire_repository_version[0:7]

_can_lazy_import = False

if "SIRE_NO_LAZY_IMPORT" not in _os.environ:
    try:
        import lazy_import as _lazy_import
        import logging as _logging

        _logger = _logging.getLogger("lazy_import")
        _logger.setLevel(_logging.ERROR)

        # Previously needed to filter to remove excessive warnings
        # from 'frozen importlib' when lazy loading.
        # import warnings
        # warnings.filterwarnings("ignore")

        _can_lazy_import = True

    except Exception as e:
        print("Lazy import disabled")
        print(e)
        _can_lazy_import = False


# Lazy import the modules for speed, and also to prevent pythonizing them
# if the users wants to run in legacy mode
if _can_lazy_import:
    analysis = _lazy_import.lazy_module("sire.analysis")
    base = _lazy_import.lazy_module("sire.base")
    cas = _lazy_import.lazy_module("sire.cas")
    convert = _lazy_import.lazy_module("sire.convert")
    cluster = _lazy_import.lazy_module("sire.cluster")
    error = _lazy_import.lazy_module("sire.error")
    ff = _lazy_import.lazy_module("sire.ff")
    id = _lazy_import.lazy_module("sire.id")
    io = _lazy_import.lazy_module("sire.io")
    maths = _lazy_import.lazy_module("sire.maths")
    mm = _lazy_import.lazy_module("sire.mm")
    mol = _lazy_import.lazy_module("sire.mol")
    morph = _lazy_import.lazy_module("sire.morph")
    move = _lazy_import.lazy_module("sire.move")
    options = _lazy_import.lazy_module("sire.options")
    qm = _lazy_import.lazy_module("sire.qm")
    qt = _lazy_import.lazy_module("sire.qt")
    restraints = _lazy_import.lazy_module("sire.restraints")
    search = _lazy_import.lazy_module("sire.search")
    squire = _lazy_import.lazy_module("sire.squire")
    stream = _lazy_import.lazy_module("sire.stream")
    units = _lazy_import.lazy_module("sire.units")
    utils = _lazy_import.lazy_module("sire.utils")
    vol = _lazy_import.lazy_module("sire.vol")


def _version_string():
    """
    Return a nicely formatted string that describes
    the current Sire version
    """
    from .base import (
        get_release_version,
        get_repository_branch,
        get_repository_version_is_clean,
    )

    from .config import sire_repository_version

    return """Sire %s [%s|%s, %s]""" % (
        get_release_version(),
        get_repository_branch(),
        sire_repository_version[0:7],
        ["unclean", "clean"][get_repository_version_is_clean()],
    )


config.version_string = _version_string
