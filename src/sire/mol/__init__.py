__all__ = [
    "get_alignment",
    "is_water",
    "selection_to_atoms",
    "Atom",
    "AtomIdx",
    "AtomMapping",
    "AtomMatch",
    "AtomMatchM",
    "AtomName",
    "AtomNum",
    "BondOrder",
    "Chain",
    "ChainIdx",
    "ChainName",
    "Cursor",
    "Cursors",
    "CursorsM",
    "Dynamics",
    "Element",
    "Minimisation",
    "ResIdx",
    "ResName",
    "ResNum",
    "Residue",
    "MolIdx",
    "Molecule",
    "MolName",
    "MolNum",
    "SegIdx",
    "Segment",
    "SegName",
    "Selector_Atom_",
    "Selector_Chain_",
    "Selector_Residue_",
    "Selector_Segment_",
    "SelectorM_Atom_",
    "SelectorM_Chain_",
    "SelectorM_Residue_",
    "SelectorM_Segment_",
    "SelectorMol",
    "Stereochemistry",
    "TrajectoryIterator",
]


from ..legacy import Mol as _Mol
from .. import use_new_api as _use_new_api

# Imported as we need to ensure that Base is pythonized before
# loading the rest of this module
from ..legacy import Base as _Base  # noqa: F401


from ..legacy.Mol import (
    AtomName,
    AtomNum,
    AtomIdx,
    AtomID,
    ResName,
    ResNum,
    ResIdx,
    ResID,
    ChainName,
    ChainIdx,
    ChainID,
    SegName,
    SegIdx,
    SegID,
    MolName,
    MolNum,
    MolIdx,
    BondID,
    AngleID,
    DihedralID,
    ImproperID,
    Atom,
    Selector_Atom_,
    SelectorM_Atom_,
    CutGroup,
    Selector_CutGroup_,
    SelectorM_CutGroup_,
    Residue,
    Selector_Residue_,
    SelectorM_Residue_,
    Chain,
    Selector_Chain_,
    SelectorM_Chain_,
    Segment,
    Selector_Segment_,
    SelectorM_Segment_,
    Molecule,
    SelectorMol,
    MoleculeView,
    Select,
    BondOrder,
    Stereochemistry,
    AtomCoords,
    AtomMapping,
    AtomMatch,
    AtomMatchM,
)

from ._cursor import Cursor, Cursors, CursorsM
from ._minimisation import Minimisation
from ._dynamics import Dynamics
from ._trajectory import TrajectoryIterator
from ._element import Element
from ._view import view as _viewfunc
from ._smiles import (
    _to_smiles,
    _to_smarts,
    _view2d,
    _selector_to_smiles,
    _selector_to_smarts,
    _selector_view2d,
)

# Imported to make sure that Vector has units attached
from ..maths import Vector as _Vector  # noqa: F401

_use_new_api()


try:
    get_alignment = _Mol.getAlignment
except AttributeError:
    get_alignment = _Mol.get_alignment


def selection_to_atoms(mols, atoms):
    """
    Convert the passed selection to a list of atoms.

    Parameters
    ----------
    mols : A molecule view or collection
        The molecule container from which to select the atoms.
    atoms : str, int, list, molecule view/collection etc.
        Any valid search string, atom index, list of atom indicies,
        or molecule view/container that can be used to select
        atoms from ``mols``

    Returns
    -------
    atoms : A SelectorM_Atoms_ or Selector_Atom_ containing the list
        of atoms.

    Examples
    --------

    >>> import sire as sr
    >>> mols = sr.load(sr.expand(sr.tutorial_url, "ala.top", "ala.crd"))
    >>> sr.mol.selection_to_atoms(mols, "atomname CA")
    Selector<SireMol::Atom>( size=1
    0:  Atom( CA:9    [  16.54,    5.03,   15.81] )
    )

    >>> sr.mol.selection_to_atoms(mols, "resname ALA")
    Selector<SireMol::Atom>( size=10
    0:  Atom( N:7     [  17.22,    4.31,   14.71] )
    1:  Atom( H:8     [  16.68,    3.62,   14.22] )
    2:  Atom( CA:9    [  16.54,    5.03,   15.81] )
    3:  Atom( HA:10   [  17.29,    5.15,   16.59] )
    4:  Atom( CB:11   [  16.05,    6.39,   15.26] )
    5:  Atom( HB1:12  [  15.63,    6.98,   16.07] )
    6:  Atom( HB2:13  [  16.90,    6.89,   14.80] )
    7:  Atom( HB3:14  [  15.24,    6.18,   14.55] )
    8:  Atom( C:15    [  15.37,    4.19,   16.43] )
    9:  Atom( O:16    [  14.94,    3.17,   15.88] )
    )

    >>> sr.mol.selection_to_atoms(mols, [0, 1, 2, 3])
    SireMol::SelectorM<SireMol::Atom>( size=4
    0: MolNum(641) Atom( HH31:1  [  18.45,    3.49,   12.44] )
    1: MolNum(641) Atom( CH3:2   [  18.98,    3.45,   13.39] )
    2: MolNum(641) Atom( HH32:3  [  20.05,    3.63,   13.29] )
    3: MolNum(641) Atom( HH33:4  [  18.80,    2.43,   13.73] )
    )

    >>> sr.mol.selection_to_atoms(mols, mols[0])
    Selector<SireMol::Atom>( size=22
    0:  Atom( HH31:1  [  18.45,    3.49,   12.44] )
    1:  Atom( CH3:2   [  18.98,    3.45,   13.39] )
    2:  Atom( HH32:3  [  20.05,    3.63,   13.29] )
    3:  Atom( HH33:4  [  18.80,    2.43,   13.73] )
    4:  Atom( C:5     [  18.48,    4.55,   14.35] )
    ...
    17:  Atom( H:18    [  15.34,    5.45,   17.96] )
    18:  Atom( CH3:19  [  13.83,    3.94,   18.35] )
    19:  Atom( HH31:20 [  14.35,    3.41,   19.15] )
    20:  Atom( HH32:21 [  13.19,    4.59,   18.94] )
    21:  Atom( HH33:22 [  13.21,    3.33,   17.69] )
    )
    """
    if hasattr(atoms, "atoms"):
        return atoms.atoms()

    from ..legacy.Base import NumberProperty, IntegerArrayProperty

    if type(atoms) is int or type(atoms) is list:
        return mols.atoms()[atoms]
    elif type(atoms) is NumberProperty or type(atoms) is IntegerArrayProperty:
        atoms = [x.as_integer() for x in atoms.as_array()]
        return mols.atoms()[atoms]
    else:
        return mols[atoms].atoms()


def is_water(mol, map=None):
    """
    Return whether or not the passed molecule (or collection of molecules)
    are water. This returns a single True or False if a single molecule
    is passed, or a list of True/False values if a molecule collection
    is passed.
    """
    from ..base import create_map

    map = create_map(map)

    try:
        return _Mol.is_water(mol, map)
    except Exception:
        pass

    if hasattr(type(mol), "molecules"):
        return _Mol.is_water(mol.molecules(), map)

    raise TypeError(
        f"Cannot convert {mol} to a molecule view or molecule collection."
    )


# Here I will define some functions that make accessing
# things from moleculeviews more convenient
def __is_molecule_class(obj):
    mro = type(obj).mro()

    return Molecule in mro or SelectorMol in mro


def __is_bond_class(obj):
    mro = type(obj).mro()

    from sire.mm import Bond, SelectorBond

    return Bond in mro or SelectorBond in mro


def __is_angle_class(obj):
    mro = type(obj).mro()

    from sire.mm import Angle, SelectorAngle

    return Angle in mro or SelectorAngle in mro


def __is_dihedral_class(obj):
    mro = type(obj).mro()

    from sire.mm import Dihedral, SelectorDihedral

    return Dihedral in mro or SelectorDihedral in mro


def __is_improper_class(obj):
    mro = type(obj).mro()

    from sire.mm import Improper, SelectorImproper

    return Improper in mro or SelectorImproper in mro


def __is_atom_class(obj):
    mro = type(obj).mro()

    return Atom in mro or Selector_Atom_ in mro or SelectorM_Atom_ in mro


def __is_residue_class(obj):
    mro = type(obj).mro()

    return (
        Residue in mro or Selector_Residue_ in mro or SelectorM_Residue_ in mro
    )


def __is_chain_class(obj):
    mro = type(obj).mro()

    return Chain in mro or Selector_Chain_ in mro or SelectorM_Chain_ in mro


def __is_segment_class(obj):
    mro = type(obj).mro()

    return (
        Segment in mro or Selector_Segment_ in mro or SelectorM_Segment_ in mro
    )


def __is_cutgroup_class(obj):
    mro = type(obj).mro()

    return (
        CutGroup in mro
        or Selector_CutGroup_ in mro
        or SelectorM_CutGroup_ in mro
    )


def __is_selector_class(obj):
    try:
        t = obj.what()
        return (
            t.find("SireMol::Selector") != -1
            or t.find("SireMM::Selector") != -1
        )
    except Exception:
        return False


def __is_internal_class(obj):
    from ..mm import Bond, Angle, Dihedral, Improper

    return type(obj) in [Bond, Angle, Dihedral, Improper]


def __is_list_class(obj):
    if type(obj) is list:
        return True
    else:
        try:
            return obj.what().find("::Selector") != -1
        except Exception:
            return False


def __fix_obj(obj):
    """This is needed because MolViewPtr objects that hold Selector_T_ types
    do not convert properly.
    """
    w = obj.what()

    if w == Molecule.typename():
        return obj
    elif w == Selector_Atom_.typename():
        return obj.atoms()
    elif w == Selector_Residue_.typename():
        return obj.residues()
    elif w == Selector_Chain_.typename():
        return obj.chains()
    elif w == Selector_Segment_.typename():
        return obj.segments()
    elif w == Selector_CutGroup_.typename():
        return obj.cutgroups()
    else:
        return obj


def __from_select_result(obj):
    """Convert the passed SelectResult from a search into the
    most appropriate MoleculeView-derived class
    """
    if hasattr(obj, "listCount") and not hasattr(obj, "list_count"):
        # Sometimes the SelectResult hasn't been converted, i.e. because
        # it has come from an old api or mixed version of Sire, or
        # BioSimSpace has done something weird...
        raise SystemError(
            "Something has gone wrong with sire. Despite being loaded "
            "with the new or mixed API, it is being passed the object "
            f"'{obj}' of type {type(obj)} which only has the old API active. "
            "Has Sire been loaded with support for old module names?"
        )

    if obj.list_count() == 0:
        raise KeyError("Nothing matched the search.")

    typ = obj.get_common_type()

    if obj.list_count() == 1:
        obj = __fix_obj(obj.list_at(0))

        if obj.what() in [
            "SireMM::SelectorBond",
            "SireMM::SelectorAngle",
            "SireMM::SelectorDihedral",
            "SireMM::SelectorImproper",
        ]:
            if obj.count() == 1:
                obj = obj[0]
        elif obj.what() != typ:
            if typ == Molecule.typename():
                return obj.molecule()
            elif typ == Segment.typename():
                return obj.segments(auto_reduce=True)
            elif typ == Chain.typename():
                return obj.chains(auto_reduce=True)
            elif typ == Residue.typename():
                return obj.residues(auto_reduce=True)
            elif typ == Atom.typename():
                return obj.atoms(auto_reduce=True)

        return obj

    if Atom.typename() != "SireMol::Atom":
        raise AssertionError(
            "The typename of an atom should be 'SireMol::Atom', but it "
            f"is instead {Atom.typename()}. The mro is {Atom.mro()}. "
            "This suggests something is broken with the boost wrappers, "
            "and that further strange bugs will be present!"
        )

    if typ == Molecule.typename():
        return SelectorMol(obj)
    elif typ == Atom.typename():
        return SelectorM_Atom_(obj)
    elif typ == Residue.typename():
        return SelectorM_Residue_(obj)
    elif typ == Chain.typename():
        return SelectorM_Chain_(obj)
    elif typ == Segment.typename():
        return SelectorM_Segment_(obj)
    elif typ == CutGroup.typename():
        return SelectorM_CutGroup_(obj)
    elif typ == AtomMatch.typename():
        return AtomMatchM(obj)
    else:
        from ..mm import SelectorBond, SelectorMBond

        if SelectorBond in type(obj.list_at(0)).mro():
            return SelectorMBond(obj)

        from ..mm import SelectorAngle, SelectorMAngle

        if SelectorAngle in type(obj.list_at(0)).mro():
            return SelectorMAngle(obj)

        from ..mm import SelectorDihedral, SelectorMDihedral

        if SelectorDihedral in type(obj.list_at(0)).mro():
            return SelectorMDihedral(obj)

        from ..mm import SelectorImproper, SelectorMImproper

        if SelectorImproper in type(obj.list_at(0)).mro():
            return SelectorMImproper(obj)

        # We really shouldn't get here, and if we do, then difficult
        # to diagnose errors will propogate
        from ..utils import Console

        Console.warning(
            f"Unrecognised type {typ}. "
            "Expecting an Atom, Residue or other MoleculeView "
            "type. Something may be wrong with the wrappers, and "
            "further bugs from this point onwards are likely. "
            f"Here is the container: {type(obj)}"
        )

        # return this as a raw list
        return obj.to_list()


def __select_call__(obj, molecules, map=None):
    """Search for the desired objects in the passed molecules,
    optionally passing in a property map to identify the properties
    """
    from ..system import System

    if type(molecules) is System:
        molecules = molecules._system

    from ..base import create_map

    map = create_map(map)

    return __from_select_result(obj.__orig_call__(molecules, map))


if not hasattr(Select, "__orig_call__"):
    Select.__orig_call__ = Select.__call__
    Select.__call__ = __select_call__


def __fixed__getitem__(obj, key):
    map = None

    from ..base import create_map

    try:
        if len(key) == 2:
            map = create_map(key[1])
            key = key[0]
    except Exception:
        # the second value is not a property map
        pass

    if type(key) is int:
        if __is_selector_class(obj):
            return obj.__orig__getitem__(key)
        elif __is_chain_class(obj):
            return obj.residue(key)
        elif __is_internal_class(obj):
            return obj.atom(obj.id()[key])
        else:
            return obj.atom(key)
    elif type(key) is str:
        # is this a search object - if so, then return whatever is
        # most relevant from the search
        if map is None:
            map = create_map(map)

        try:
            # try to search for the object - this will raise
            # a SyntaxError if this is not a search term
            # (and is instead a name)
            return __from_select_result(obj.search(key, map=map))
        except SyntaxError:
            # ignore SyntaxErrors as this is a name
            pass
    elif AtomID in type(key).mro():
        return obj.atoms(key, auto_reduce=True)
    elif ResID in type(key).mro():
        return obj.residues(key, auto_reduce=True)
    elif ChainID in type(key).mro():
        return obj.chains(key, auto_reduce=True)
    elif SegID in type(key).mro():
        return obj.segments(key, auto_reduce=True)
    elif BondID in type(key).mro():
        return obj.bonds(key, auto_reduce=True)
    elif AngleID in type(key).mro():
        return obj.angles(key, auto_reduce=True)
    elif DihedralID in type(key).mro():
        return obj.dihedrals(key, auto_reduce=True)
    elif ImproperID in type(key).mro():
        return obj.impropers(key, auto_reduce=True)

    if __is_selector_class(obj):
        return obj.__orig__getitem__(key)
    elif __is_chain_class(obj):
        return obj.residues(key, auto_reduce=True)
    else:
        return obj.atoms(key, auto_reduce=True)


def __fixed__atoms__(obj, idx=None, auto_reduce=False, map=None):
    from ..base import create_map

    if idx is None:
        result = obj.__orig__atoms()
    elif type(idx) is range:
        result = obj.__orig__atoms(list(idx), map=create_map(map))
    else:
        result = obj.__orig__atoms(idx, map=create_map(map))

    if auto_reduce and len(result) == 1:
        return result[0]
    else:
        return result


def __fixed__bonds__(obj, idx=None, idx1=None, auto_reduce=False, map=None):
    if idx is None and idx1 is not None:
        idx = idx1
        idx1 = None

    from . import MoleculeView
    from ..base import create_map

    if issubclass(type(obj), MoleculeView):
        # this is a single-molecule view
        from ..mm import SelectorBond

        C = SelectorBond

        def _fromBondID(obj, bondid):
            return SelectorBond(obj, bondid, map=create_map(map))

    else:
        # this is a multi-molecule container
        from ..mm import SelectorMBond

        C = SelectorMBond

        def _fromBondID(obj, bondid):
            return SelectorMBond(
                obj.to_select_result(), bondid, map=create_map(map)
            )

    if idx is None:
        try:
            result = C(obj, map=create_map(map))
        except Exception:
            result = C(obj.to_select_result(), map=create_map(map))
    elif idx1 is None:
        if BondID in type(idx).mro():
            result = _fromBondID(obj, idx)
        else:
            result = C(obj.atoms(idx, map=create_map(map)))
    else:
        result = C(obj.atoms(idx), obj.atoms(idx1), map=create_map(map))

    if auto_reduce and len(result) == 1:
        return result[0]
    else:
        return result


def __fixed__angles__(
    obj, idx=None, idx1=None, idx2=None, auto_reduce=False, map=None
):
    if idx1 is None and idx2 is not None:
        idx1 = idx2
        idx2 = None

    if idx is None and idx1 is not None:
        idx = idx1
        idx1 = None

    from . import MoleculeView
    from ..base import create_map

    if issubclass(type(obj), MoleculeView):
        # this is a single-molecule view
        from ..mm import SelectorAngle

        C = SelectorAngle

        def _fromAngleID(obj, angid):
            return SelectorAngle(obj, angid, map=create_map(map))

    else:
        # this is a multi-molecule container
        from ..mm import SelectorMAngle

        C = SelectorMAngle

        def _fromAngleID(obj, angid):
            return SelectorMAngle(
                obj.to_select_result(), angid, map=create_map(map)
            )

    if idx is None:
        try:
            result = C(obj, map=create_map(map))
        except Exception:
            result = C(obj.to_select_result(), map=create_map(map))
    elif idx1 is None:
        if AngleID in type(idx).mro():
            result = _fromAngleID(obj, idx)
        else:
            result = C(obj.atoms(idx, map=create_map(map)))
    elif idx2 is None:
        result = C(obj.atoms(idx), obj.atoms(idx1), map=create_map(map))
    else:
        result = C(
            obj.atoms(idx),
            obj.atoms(idx1),
            obj.atoms(idx2),
            map=create_map(map),
        )

    if auto_reduce and len(result) == 1:
        return result[0]
    else:
        return result


def __fixed__dihedrals__(
    obj, idx=None, idx1=None, idx2=None, idx3=None, auto_reduce=False, map=None
):
    if idx2 is None and idx3 is not None:
        idx2 = idx3
        idx3 = None

    if idx1 is None and idx2 is not None:
        idx1 = idx2
        idx2 = None

    if idx is None and idx1 is not None:
        idx = idx1
        idx1 = None

    from . import MoleculeView
    from ..base import create_map

    if issubclass(type(obj), MoleculeView):
        # this is a single-molecule view
        from ..mm import SelectorDihedral

        C = SelectorDihedral

        def _fromDihedralID(obj, dihid):
            return SelectorDihedral(obj, dihid, map=create_map(map))

    else:
        # this is a multi-molecule container
        from ..mm import SelectorMDihedral

        C = SelectorMDihedral

        def _fromDihedralID(obj, dihid):
            return SelectorMDihedral(
                obj.to_select_result(), dihid, map=create_map(map)
            )

    if idx is None:
        try:
            result = C(obj, map=create_map(map))
        except Exception:
            result = C(obj.to_select_result(), map=create_map(map))
    elif idx1 is None:
        if DihedralID in type(idx).mro():
            result = _fromDihedralID(obj, idx)
        else:
            result = C(obj.atoms(idx), map=create_map(map))
    elif idx2 is None:
        result = C(obj.atoms(idx), obj.atoms(idx1), map=create_map(map))
    elif idx3 is None:
        result = C(
            obj.atoms(idx),
            obj.atoms(idx1),
            obj.atoms(idx2),
            map=create_map(map),
        )
    else:
        result = C(
            obj.atoms(idx),
            obj.atoms(idx1),
            obj.atoms(idx2),
            obj.atoms(idx3),
            map=create_map(map),
        )

    if auto_reduce and len(result) == 1:
        return result[0]
    else:
        return result


def __fixed__impropers__(
    obj, idx=None, idx1=None, idx2=None, idx3=None, auto_reduce=False, map=None
):
    if idx2 is None and idx3 is not None:
        idx2 = idx3
        idx3 = None

    if idx1 is None and idx2 is not None:
        idx1 = idx2
        idx2 = None

    if idx is None and idx1 is not None:
        idx = idx1
        idx1 = None

    from . import MoleculeView
    from ..base import create_map

    if issubclass(type(obj), MoleculeView):
        # this is a single-molecule view
        from ..mm import SelectorImproper

        C = SelectorImproper

        def _fromImproperID(obj, impid):
            return SelectorImproper(obj, impid, map=create_map(map))

    else:
        # this is a multi-molecule container
        from ..mm import SelectorMImproper

        C = SelectorMImproper

        def _fromImproperID(obj, impid):
            return SelectorMImproper(
                obj.to_select_result(), impid, map=create_map(map)
            )

    if idx is None:
        try:
            result = C(obj, map=create_map(map))
        except Exception:
            result = C(obj.to_select_result(), map=create_map(map))
    elif idx1 is None:
        if ImproperID in type(idx).mro():
            result = _fromImproperID(obj, idx)
        else:
            result = C(obj.atoms(idx), map=create_map(map))
    elif idx2 is None:
        result = C(obj.atoms(idx), obj.atoms(idx1), map=create_map(map))
    elif idx3 is None:
        result = C(
            obj.atoms(idx),
            obj.atoms(idx1),
            obj.atoms(idx2),
            map=create_map(map),
        )
    else:
        result = C(
            obj.atoms(idx),
            obj.atoms(idx1),
            obj.atoms(idx2),
            obj.atoms(idx3),
            map=create_map(map),
        )

    if auto_reduce and len(result) == 1:
        return result[0]
    else:
        return result


def __fixed__bond__(obj, idx=None, idx1=None, map=None):
    bonds = __fixed__bonds__(obj, idx, idx1, auto_reduce=False, map=map)

    if len(bonds) == 0:
        raise KeyError("There is no matching bond in this view.")
    elif len(bonds) > 1:
        raise KeyError(
            f"More than one bond matches. Number of matches is {len(bonds)}."
        )

    return bonds[0]


def __fixed__angle__(obj, idx=None, idx1=None, idx2=None, map=None):
    angles = __fixed__angles__(
        obj, idx, idx1, idx2, auto_reduce=False, map=map
    )

    if len(angles) == 0:
        raise KeyError("There is no matching angle in this view.")
    elif len(angles) > 1:
        raise KeyError(
            f"More than one angle matches. Number of matches is {len(angles)}."
        )

    return angles[0]


def __fixed__dihedral__(
    obj, idx=None, idx1=None, idx2=None, idx3=None, map=None
):
    dihedrals = __fixed__dihedrals__(
        obj, idx, idx1, idx2, idx3, auto_reduce=False, map=map
    )

    if len(dihedrals) == 0:
        raise KeyError("There is no matching dihedral in this view.")
    elif len(dihedrals) > 1:
        raise KeyError(
            "More than one dihedral matches. Number of "
            f"matches is {len(dihedrals)}."
        )

    return dihedrals[0]


def __fixed__improper__(
    obj, idx=None, idx1=None, idx2=None, idx3=None, map=None
):
    impropers = __fixed__impropers__(
        obj, idx, idx1, idx2, idx3, auto_reduce=False, map=map
    )

    if len(impropers) == 0:
        raise KeyError("There is no matching improper in this view.")
    elif len(impropers) > 1:
        raise KeyError(
            "More than one improper matches. Number of "
            f"matches is {len(impropers)}."
        )

    return impropers[0]


def __fixed__residues__(obj, idx=None, auto_reduce=False, map=None):
    from ..base import create_map

    if idx is None:
        result = obj.__orig__residues()
    elif type(idx) is range:
        result = obj.__orig__residues(list(idx), map=create_map(map))
    else:
        result = obj.__orig__residues(idx, map=create_map(map))

    if auto_reduce and len(result) == 1:
        return result[0]
    else:
        return result


def __fixed__chains__(obj, idx=None, auto_reduce=False, map=None):
    from ..base import create_map

    if idx is None:
        result = obj.__orig__chains()
    elif type(idx) is range:
        result = obj.__orig__chains(list(idx), map=create_map(map))
    else:
        result = obj.__orig__chains(idx, map=create_map(map))

    if auto_reduce and len(result) == 1:
        return result[0]
    else:
        return result


def __fixed__segments__(obj, idx=None, auto_reduce=False, map=None):
    from ..base import create_map

    if idx is None:
        result = obj.__orig__segments()
    elif type(idx) is range:
        result = obj.__orig__segments(list(idx), map=create_map(map))
    else:
        from ..base import create_map

        result = obj.__orig__segments(idx, map=create_map(map))

    if auto_reduce and len(result) == 1:
        return result[0]
    else:
        return result


def __fixed__molecules__(obj, idx=None, auto_reduce=False, map=None):
    from ..base import create_map

    if idx is None:
        result = obj.__orig__molecules()
    elif type(idx) is range:
        result = obj.__orig__molecules(list(idx), map=create_map(map))
    else:
        result = obj.__orig__molecules(idx, map=create_map(map))

    if auto_reduce and len(result) == 1:
        return result[0]
    else:
        return result


def __fix_getitem(C):
    if not hasattr(C, "__orig__getitem__"):
        C.__orig__getitem__ = C.__getitem__

    if not hasattr(C, "__orig__atoms"):
        C.__orig__atoms = C.atoms

    if not hasattr(C, "__orig__residues"):
        C.__orig__residues = C.residues

    if not hasattr(C, "__orig__chains"):
        C.__orig__chains = C.chains

    if not hasattr(C, "__orig__segments"):
        C.__orig__segments = C.segments

    C.__getitem__ = __fixed__getitem__
    C.atoms = __fixed__atoms__
    C.residues = __fixed__residues__
    C.chains = __fixed__chains__
    C.segments = __fixed__segments__

    C.count = C.__len__

    if hasattr(C, "measure"):
        # make sure we use the right `size` function
        C.size = C.measure
    else:
        C.size = C.__len__

    if hasattr(C, "molecules"):
        if not hasattr(C, "__orig__molecules"):
            C.__orig__molecules = C.molecules

        C.molecules = __fixed__molecules__

    C.bonds = __fixed__bonds__
    C.bond = __fixed__bond__
    C.angles = __fixed__angles__
    C.angle = __fixed__angle__
    C.dihedrals = __fixed__dihedrals__
    C.dihedral = __fixed__dihedral__
    C.impropers = __fixed__impropers__
    C.improper = __fixed__improper__


try:
    Residue.__len__ = Residue.nAtoms
    Chain.__len__ = Chain.nResidues
    Segment.__len__ = Segment.nAtoms
    CutGroup.__len__ = CutGroup.nAtoms
    Molecule.__len__ = Molecule.nAtoms
except AttributeError:
    Residue.__len__ = Residue.num_atoms
    Chain.__len__ = Chain.num_residues
    Segment.__len__ = Segment.num_atoms
    CutGroup.__len__ = CutGroup.num_atoms
    Molecule.__len__ = Molecule.num_atoms

for C in [
    Atom,
    CutGroup,
    Residue,
    Chain,
    Segment,
    Molecule,
    Selector_Atom_,
    Selector_Residue_,
    Selector_Chain_,
    Selector_Segment_,
    Selector_CutGroup_,
    SelectorMol,
    SelectorM_Atom_,
    SelectorM_Residue_,
    SelectorM_Chain_,
    SelectorM_Segment_,
    SelectorM_CutGroup_,
]:
    __fix_getitem(C)


def _atom_coordinates(atom, map=None):
    if map is None:
        return atom.property("coordinates")

    from ..base import create_map

    map = create_map(map)

    return atom.property(map["coordinates"])


Atom.element = lambda x: x.property("element")
Atom.lj = lambda x: x.property("LJ")
Atom.coordinates = _atom_coordinates
Atom.coords = Atom.coordinates
Atom.x = lambda atom, map=None: atom.coordinates(map=map).x()
Atom.y = lambda atom, map=None: atom.coordinates(map=map).y()
Atom.z = lambda atom, map=None: atom.coordinates(map=map).z()


def __atomcoords__str__(obj):
    n = len(obj)

    if n == 0:
        return "AtomCoords::empty"

    parts = []

    if n <= 10:
        for i in range(0, n):
            parts.append(f"{i}: {obj[i]}")
    else:
        for i in range(0, 5):
            parts.append(f"{i}: {obj[i]}")

        parts.append("...")

        for i in range(n - 5, n):
            parts.append(f"{i}: {obj[i]}")

    joined = "\n".join(parts)

    return f"AtomCoords( size={n}\n{joined}\n)"


AtomCoords.__str__ = __atomcoords__str__
AtomCoords.__repr__ = __atomcoords__str__


def _add_evals(obj):
    from ..base import create_map

    obj.mass = lambda x, map=None: x.evaluate().mass(map=create_map(map))
    obj.charge = lambda x, map=None: x.evaluate().charge(map=create_map(map))
    obj.coordinates = lambda x, map=None: x.evaluate().center_of_mass(
        map=create_map(map)
    )
    obj.coords = obj.coordinates
    obj.x = lambda x, map=None: x.coordinates(map=map).x()
    obj.y = lambda x, map=None: x.coordinates(map=map).y()
    obj.z = lambda x, map=None: x.coordinates(map=map).z()


def _get_container_property(x, key):
    vals = []

    from ..base import ProgressBar

    with ProgressBar(len(x), "Extract property") as bar:
        for v in x:
            prop = _get_property(v, key)

            if type(prop) is list:
                vals += prop
            else:
                vals.append(prop)

            bar.tick()

    return vals


def _get_property(x, key):
    if hasattr(x, "__orig__property"):
        # get the property directly
        try:
            return x.__orig__property(key)
        except Exception as e:
            saved_exception = e
    else:
        saved_exception = None

    # we couldn't get the property directly, so get
    # the property at the molecule level...
    try:
        mol = x.molecule()
    except Exception:
        # this must be a SelectorMol or other container
        return _get_container_property(x, key)

    prop = mol.property(key)

    # Now extract the bits we want
    import sire

    if issubclass(prop.__class__, sire.legacy.Mol.AtomProp):
        vals = []
        for atom in x.atoms():
            vals.append(atom.property(key))

        return vals
    elif issubclass(prop.__class__, sire.legacy.Mol.ResProp):
        vals = []
        for res in x.residues():
            vals.append(res.property(key))

        return vals
    elif issubclass(prop.__class__, sire.legacy.Mol.ChainProp):
        vals = []
        for chain in x.chains():
            vals.append(chain.property(key))

        return vals
    elif issubclass(prop.__class__, sire.legacy.Mol.SegProp):
        vals = []
        for seg in x.segments():
            vals.append(seg.property(key))

        return vals
    elif issubclass(prop.__class__, sire.legacy.Mol.MolViewProperty):
        return prop
    elif saved_exception is not None:
        raise saved_exception
    else:
        raise KeyError(f"Could not find property {key} in {x}")


def _apply(objs, func, *args, **kwargs):
    """
    Call the passed function on all views in the container,
    appending the result to a list of results, which
    is returned.

    The function can be either;

    1. a string containing the name of the function to call, or
    2. an actual function (either a normal function or a lambda expression)

    You can optionally pass in positional and keyword arguments
    here that will be passed to the function.

    Args:
        objs (self): The container itself (this is self)
        func (str or function): The function to be called, or the name
                                of the function to be called.

    Returns:
        list: A list of the results, with one result per view in the container.
    """
    result = []

    from ..base import ProgressBar

    if str(func) == func:
        # we calling a named function
        with ProgressBar(
            total=len(objs), text="Looping through views"
        ) as progress:
            for i, obj in enumerate(objs):
                result.append(getattr(obj, func)(*args, **kwargs))
                progress.set_progress(i + 1)

    else:
        # we have been passed the function to call
        with ProgressBar(
            total=len(objs), text="Looping through views"
        ) as progress:
            for i, obj in enumerate(objs):
                result.append(func(obj, *args, **kwargs))
                progress.set_progress(i + 1)

    return result


def _apply_reduce(objs, func, reduction_func=None, *args, **kwargs):
    """
    Call the passed function on all views in the container,
    reducing the result into a single value via the 'reduce' function.

    This is equivalent to calling

    ```
    reduce(reduction_func, objs.apply(func, *args, **kwargs))

    The function can be either;

    1. a string containing the name of the function to call, or
    2. an actual function (either a normal function or a lambda expression)

    The reduction function should be a function that can be passed
    to `reduce`. If this isn't passed, then it is assumed to
    be operator.add.

    You can optionally pass in positional and keyword arguments
    here that will be passed to the applied function.

    Args:
        objs (self): The container itself (this is self)
        func (str or function): The function to be called, or the name
                                of the function to be called.
        reduction_func: The function used to reduce the result. This
                        is operator.add by default.

    Returns:
        result: The reduced result
    """
    if reduction_func is None:
        from operator import add

        reduction_func = add

    from functools import reduce

    return reduce(reduction_func, objs.apply(func, *args, **kwargs))


def _add_apply_func(obj):
    if hasattr(obj, "apply"):
        return

    obj.apply = _apply
    obj.apply_reduce = _apply_reduce


def _add_property_func(obj):
    if hasattr(obj, "__orig__property"):
        return

    if hasattr(obj, "property"):
        obj.__orig__property = obj.property

    obj.property = _get_property


for C in [
    MoleculeView,
    SelectorMol,
    SelectorM_Atom_,
    SelectorM_Residue_,
    SelectorM_Chain_,
    SelectorM_CutGroup_,
    SelectorM_Segment_,
]:
    _add_evals(C)
    _add_property_func(C)
    _add_apply_func(C)

for C in [Residue, Chain, Segment]:
    _add_property_func(C)


def _molecule(obj, id=None):
    """
    Return the molecule that contains this view. If 'id' is specified
    then this will check that the ID matches this molecule before
    returning it. This will raise an exception if the ID doesn't match.

    Returns
    -------

    sire.mol.Molecule: The molecule containing this view
    """
    if id is None:
        return obj.__orig__molecule()
    else:
        from . import SelectorMol

        return SelectorMol(obj).molecule(id)


def _molecules(obj, id=None):
    """
    Return the molecule that contains this view as a list of molecules,
    containing just a single molecule.

    If 'id' is specified then this checks that this molecule matches
    the ID before returning it. This will raise an exception if the
    ID doesn't match.

    Returns
    -------

    sire.mol.SelectorMol: A collection containing only a single molecule
                          that contains this view.
    """
    from . import SelectorMol

    if id is None:
        return SelectorMol(obj).molecules()
    else:
        return SelectorMol(obj).molecules(id)


if not hasattr(MoleculeView, "__orig__molecule"):
    MoleculeView.__orig__molecule = MoleculeView.molecule
    MoleculeView.molecule = _molecule
    MoleculeView.molecules = _molecules


def _get_atom_mass(x):
    if x.has_property("mass"):
        return x.property("mass")
    elif x.has_property("element"):
        return x.property("element").mass()
    else:
        return 0


Atom.mass = _get_atom_mass


def _get_atom_charge(x):
    if x.has_property("charge"):
        return x.property("charge")
    elif x.has_property("formal_charge"):
        return x.property("formal_charge")
    else:
        return 0


Atom.charge = _get_atom_charge

Molecule.connectivity = lambda x: x.property("connectivity")


def __molecule__add__(mol0, mol1):
    """
    Combine two molecules together into a system.
    """
    from ..system import System

    mols = System()
    mols.add(mol0)
    mols.add(mol1)

    from ..vol import Cartesian

    mols.set_space(Cartesian())
    mols.set_time(0)

    return mols


Molecule.__add__ = __molecule__add__


# Here are some extra classes / functions defined as part of the
# public API


def _cursor(view, map=None):
    """
    Return a Cursor that can be used to edit the properties
    of this view
    """
    from ..base import create_map

    return Cursor(view, map=create_map(map))


Atom.cursor = _cursor
Residue.cursor = _cursor
Chain.cursor = _cursor
Segment.cursor = _cursor
Molecule.cursor = _cursor


def _dynamics(
    view,
    cutoff=None,
    cutoff_type=None,
    timestep=None,
    save_frequency=None,
    energy_frequency=None,
    frame_frequency=None,
    save_velocities=None,
    constraint=None,
    perturbable_constraint=None,
    include_constrained_energies: bool = True,
    schedule=None,
    lambda_value=None,
    swap_end_states=None,
    ignore_perturbations=None,
    temperature=None,
    pressure=None,
    vacuum=None,
    shift_delta=None,
    coulomb_power=None,
    restraints=None,
    fixed=None,
    platform=None,
    device=None,
    precision=None,
    map=None,
):
    """
    Return a Dynamics object that can be used to perform
    dynamics of the molecule(s) in this view

    cutoff: Length
        The size of the non-bonded cutoff

    cutoff_type: str
        The type of cutoff to use, e.g. "PME", "RF" etc.
        See https://sire.openbiosim.org/cheatsheet/openmm.html#choosing-options
        for the full list of options

    timestep: time
        The size of the dynamics timestep

    save_frequency: time
        The amount of simulation time between saving energies and frames.
        This can be overridden using `energy_frequency` or `frame_frequency`,
        or by these options in individual dynamics runs. Set this
        to zero if you don't want any saves.

    energy_frequency: time
        The amount of time between saving energies. This overrides the
        value in `save_frequency`. Set this to zero if you don't want
        to save energies during the trajectory. This can be overridden
        by setting energy_frequency during an individual run.

    frame_frequency: time
        The amount of time between saving frames. This overrides the
        value in `save_frequency`. Set this to zero if you don't want
        to save frames during the trajectory. This can be overridden
        by setting frame_frequency during an individual run.

    save_velocities: bool
        Whether or not to save velocities when saving trajectory frames
        during the simulation. This defaults to False, as velocity
        trajectories aren't often needed, and they double the amount
        of space that is required for a trajectory.

    constraint: str
        The type of constraint to use for bonds and/or angles, e.g.
        `h-bonds`, `bonds` etc.
        See https://sire.openbiosim.org/cheatsheet/openmm.html#choosing-options
        for the full list of options. This will be automatically
        guessed from the timestep if it isn't set.

    perturbable_constraint: str
        The type of constraint to use for perturbable bonds and/or angles,
        e.g. `h-bonds`, `bonds` etc.
        See https://sire.openbiosim.org/cheatsheet/openmm.html#choosing-options
        for the full list of options. This equal the value of `constraint`
        if it isn't set.

    include_constrained_energies: bool
        Whether or not to include the energies of the constrained bonds
        and angles. If this is False, then the internal bond or angle
        energy of the constrained degrees of freedom are not included
        in the total energy, and their forces are not evaluated.

    schedule: sire.cas.LambdaSchedule
        The schedule used to control how perturbable forcefield parameters
        should be morphed as a function of lambda. If this is not set
        then a sire.cas.LambdaSchedule.standard_morph() is used.

    lambda_value: float
        The value of lambda at which to run dynamics. This only impacts
        perturbable molecules, whose forcefield parameters will be
        scaled according to the lambda schedule for the specified
        value of lambda.

    swap_end_states: bool
        Whether or not to swap the end states. If this is True, then
        the perturbation will run from the perturbed back to the
        reference molecule (the perturbed molecule will be at lambda=0,
        while the reference molecule will be at lambda=1). This will
        use the coordinates of the perturbed molecule as the
        starting point.

    ignore_perturbations: bool
        Whether or not to ignore perturbations. If this is True, then
        the perturbation will be ignored, and the simulation will
        be run using the properties of the reference molecule
        (or the perturbed molecule if swap_end_states is True). This
        is useful if you just want to run standard molecular dynamics
        of the reference or perturbed states.

    temperature: temperature
        The temperature at which to run the simulation. A
        microcanonical (NVE) simulation will be run if you don't
        specify the temperature.

    pressure: pressure
        The pressure at which to run the simulation. A
        microcanonical (NVE) or canonical (NVT) simulation will be
        run if the pressure is not set.

    vacuum: bool (optional)
        Whether or not to run the simulation in vacuum. If this is
        set to `True`, then the simulation space automatically be
        replaced by a `sire.vol.Cartesian` space, and the
        simulation run in vacuum.

    shift_delta: length
        The shift_delta parameter that controls the electrostatic
        and van der Waals softening potential that smooths the
        creation and deletion of ghost atoms during a potential.
        This defaults to 2.0 A.

    coulomb_power: int
        The coulomb power parmeter that controls the electrostatic
        softening potential that smooths the creation and deletion
        of ghost atoms during a potential. This defaults to 0.

    restraints: sire.mm.Restraints or list[sire.mm.Restraints]
        A single set of restraints, or a list of sets of
        restraints that will be applied to the atoms during
        the simulation.

    fixed: molecule(s) view, search string, int, list[int] etc
        Anything that can be used to identify the atom or atoms
        that should be fixed in place during the simulation. These
        atoms will not be moved by dynamics.

    platform: str
        The name of the OpenMM platform on which to run the dynamics,
        e.g. "CUDA", "OpenCL", "Metal" etc.

    device: str or int
        The ID of the GPU (or accelerator) used to accelerate
        the simulation. This would be CUDA_DEVICE_ID or similar
        if CUDA was used. This can be any valid OpenMM device string

    precision: str
        The desired precision for the simulation (e.g. `single`,
        `mixed` or `double`)

    map: dict
        A dictionary of additional options. Note that any options
        set in this dictionary that are also specified via one of
        the arguments above will be overridden by the argument
        value
    """
    from ..base import create_map
    from ..system import System
    from .. import u

    map = create_map(map)

    if vacuum is True:
        from ..vol import Cartesian

        map.set("space", Cartesian())

        view = view.clone()

        try:
            view.set_property("space", Cartesian())
        except Exception:
            pass

    if not map.specified("space"):
        map = create_map(map, {"space": "space"})

    if (
        System.is_system(view)
        and map.specified("space")
        and not map["space"].has_value()
        and not view.shared_properties().has_property(map["space"])
    ):
        # space is not a shared property, so may be lost when we
        # convert to molecules. Make sure this doesn't happen by
        # adding the space directly to the property map
        try:
            map.set("space", view.property(map["space"]))
        except Exception:
            from ..vol import Cartesian

            map.set("space", Cartesian())

    # Set default values if these have not been set
    if cutoff is None and not map.specified("cutoff"):
        from ..units import angstrom

        cutoff = 7.5 * angstrom

    if cutoff_type is None and not map.specified("cutoff_type"):
        try:
            if view.property(map["space"]).is_periodic():
                cutoff_type = "PME"
            else:
                cutoff_type = "RF"
        except Exception:
            # no space, use RF
            cutoff_type = "RF"

    if timestep is None and not map.specified("timestep"):
        from ..units import femtosecond

        timestep = 1 * femtosecond
    else:
        timestep = u(timestep)

    if save_frequency is None and not map.specified("save_frequency"):
        from ..units import picosecond

        map.set("save_frequency", 25 * picosecond)
    elif save_frequency is not None:
        map.set("save_frequency", u(save_frequency))

    if energy_frequency is not None:
        map.set("energy_frequency", u(energy_frequency))

    if frame_frequency is not None:
        map.set("frame_frequency", u(frame_frequency))

    if save_velocities is not None:
        map.set("save_velocities", save_velocities)

    if constraint is None and not map.specified("constraint"):
        map.set("constraint", "auto")

    if perturbable_constraint is not None:
        perturbable_constraint = str(perturbable_constraint).lower()

    if temperature is not None:
        temperature = u(temperature)
        map.set("temperature", temperature)

    if pressure is not None:
        pressure = u(pressure)
        map.set("pressure", pressure)

    if device is not None:
        map.set("device", str(device))

    if precision is not None:
        map.set("precision", str(precision).lower())

    if include_constrained_energies is not None:
        map.set("include_constrained_energies", include_constrained_energies)

    if platform is not None:
        map.set("platform", str(platform))

    return Dynamics(
        view,
        cutoff=cutoff,
        cutoff_type=cutoff_type,
        timestep=timestep,
        constraint=str(constraint),
        perturbable_constraint=perturbable_constraint,
        schedule=schedule,
        lambda_value=lambda_value,
        shift_delta=shift_delta,
        coulomb_power=coulomb_power,
        swap_end_states=swap_end_states,
        ignore_perturbations=ignore_perturbations,
        restraints=restraints,
        fixed=fixed,
        map=map,
    )


def _minimisation(
    view,
    cutoff=None,
    cutoff_type=None,
    constraint=None,
    perturbable_constraint=None,
    include_constrained_energies: bool = True,
    schedule=None,
    lambda_value=None,
    swap_end_states=None,
    ignore_perturbations=None,
    vacuum=None,
    shift_delta=None,
    coulomb_power=None,
    platform=None,
    device=None,
    precision=None,
    restraints=None,
    fixed=None,
    map=None,
):
    """
    Return a Minimisation object that can be used to perform
    minimisation of the molecule(s) in this view

    cutoff: Length
        The size of the non-bonded cutoff

    cutoff_type: str
        The type of cutoff to use, e.g. "PME", "RF" etc.
        See https://sire.openbiosim.org/cheatsheet/openmm.html#choosing-options
        for the full list of options

    constraint: str
        The type of constraint to use for bonds and/or angles, e.g.
        `h-bonds`, `bonds` etc.
        See https://sire.openbiosim.org/cheatsheet/openmm.html#choosing-options
        for the full list of options. This is `none` if it is not set.

    perturbable_constraint: str
        The type of constraint to use for perturbable bonds and/or angles,
        e.g. `h-bonds`, `bonds` etc.
        See https://sire.openbiosim.org/cheatsheet/openmm.html#choosing-options
        for the full list of options. This equal the value of `constraint`
        if it isn't set.

    include_constrained_energies: bool
        Whether or not to include the energies of the perturbable bonds
        and angles. If this is False, then the internal bond or angle
        energy of the perturbable degrees of freedom are not included
        in the total energy, and their forces are not evaluated.

    schedule: sire.cas.LambdaSchedule
        The schedule used to control how perturbable forcefield parameters
        should be morphed as a function of lambda. If this is not set
        then a sire.cas.LambdaSchedule.standard_morph() is used.

    lambda_value: float
        The value of lambda at which to run minimisation. This only impacts
        perturbable molecules, whose forcefield parameters will be
        scaled according to the lambda schedule for the specified
        value of lambda.

    swap_end_states: bool
        Whether or not to swap the end states. If this is True, then
        the perturbation will run from the perturbed back to the
        reference molecule (the perturbed molecule will be at lambda=0,
        while the reference molecule will be at lambda=1). This will
        use the coordinates of the perturbed molecule as the
        starting point.

    ignore_perturbations: bool
        Whether or not to ignore perturbations. If this is True, then
        the perturbation will be ignored, and the simulation will
        be run using the properties of the reference molecule
        (or the perturbed molecule if swap_end_states is True). This
        is useful if you just want to run standard molecular dynamics
        of the reference or perturbed states.

    vacuum: bool (optional)
        Whether or not to run the simulation in vacuum. If this is
        set to `True`, then the simulation space automatically be
        replaced by a `sire.vol.Cartesian` space, and the
        simulation run in vacuum.

    shift_delta: length
        The shift_delta parameter that controls the electrostatic
        and van der Waals softening potential that smooths the
        creation and deletion of ghost atoms during a potential.
        This defaults to 2.0 A.

    coulomb_power: int
        The coulomb power parmeter that controls the electrostatic
        softening potential that smooths the creation and deletion
        of ghost atoms during a potential. This defaults to 0.

    restraints: sire.mm.Restraints or list[sire.mm.Restraints]
        A single set of restraints, or a list of sets of
        restraints that will be applied to the atoms during
        the simulation.

    fixed: molecule(s) view, search string, int, list[int] etc
        Anything that can be used to identify the atom or atoms
        that should be fixed in place during the simulation. These
        atoms will not be moved by minimisation.

    platform: str
        The name of the OpenMM platform on which to run the dynamics,
        e.g. "CUDA", "OpenCL", "Metal" etc.

    device: str or int
        The ID of the GPU (or accelerator) used to accelerate
        minimisation. This would be CUDA_DEVICE_ID or similar
        if CUDA was used. This can be any valid OpenMM device string

    precision: str
        The desired precision for the simulation (e.g. `single`,
        `mixed` or `double`)

    map: dict
        A dictionary of additional options. Note that any options
        set in this dictionary that are also specified via one of
        the arguments above will be overridden by the argument
        value
    """
    from ..base import create_map
    from ..system import System

    map = create_map(map)

    if vacuum is True:
        from ..vol import Cartesian

        map.set("space", Cartesian())

        view = view.clone()

        try:
            view.set_property("space", Cartesian())
        except Exception:
            pass

    if not map.specified("space"):
        map = create_map(map, {"space": "space"})

    if (
        System.is_system(view)
        and map.specified("space")
        and not map["space"].has_value()
        and not view.shared_properties().has_property(map["space"])
    ):
        # space is not a shared property, so may be lost when we
        # convert to molecules. Make sure this doens't happen by
        # adding the space directly to the property map
        try:
            map.set("space", view.property(map["space"]))
        except Exception:
            from ..vol import Cartesian

            map.set("space", Cartesian())

    # Set default values if these have not been set
    if cutoff is None and not map.specified("cutoff"):
        from ..units import angstrom

        cutoff = 7.5 * angstrom

    if cutoff_type is None and not map.specified("cutoff_type"):
        try:
            if view.property(map["space"]).is_periodic():
                cutoff_type = "PME"
            else:
                cutoff_type = "RF"
        except Exception:
            # no space, use RF
            cutoff_type = "RF"

    if device is not None:
        map.set("device", str(device))

    if precision is not None:
        map.set("precision", str(precision).lower())

    if constraint is not None:
        map.set("constraint", str(constraint).lower())

    if perturbable_constraint is not None:
        map.set("perturbable_constraint", str(perturbable_constraint).lower())

    if include_constrained_energies is not None:
        map.set("include_constrained_energies", include_constrained_energies)

    if platform is not None:
        map.set("platform", str(platform))

    return Minimisation(
        view,
        cutoff=cutoff,
        cutoff_type=cutoff_type,
        schedule=schedule,
        lambda_value=lambda_value,
        swap_end_states=swap_end_states,
        ignore_perturbations=ignore_perturbations,
        shift_delta=shift_delta,
        coulomb_power=coulomb_power,
        restraints=restraints,
        fixed=fixed,
        map=map,
    )


MoleculeView.dynamics = _dynamics
SelectorM_Atom_.dynamics = _dynamics
SelectorM_Residue_.dynamics = _dynamics
SelectorM_Chain_.dynamics = _dynamics
SelectorM_Segment_.dynamics = _dynamics
SelectorM_CutGroup_.dynamics = _dynamics
SelectorMol.dynamics = _dynamics

MoleculeView.minimisation = _minimisation
SelectorM_Atom_.minimisation = _minimisation
SelectorM_Residue_.minimisation = _minimisation
SelectorM_Chain_.minimisation = _minimisation
SelectorM_Segment_.minimisation = _minimisation
SelectorM_CutGroup_.minimisation = _minimisation
SelectorMol.minimisation = _minimisation


def _cursors(views, map=None):
    """Return the Cursors object that contains cursors for all
    of the views in this collection. Note that the `parent`
    cursor of this list will be the molecule
    """
    cursor = views.molecule().cursor(map=map)
    return cursor._from_views(views)


Selector_Atom_.cursor = _cursors
Selector_Residue_.cursor = _cursors
Selector_Chain_.cursor = _cursors
Selector_Segment_.cursor = _cursors


def _trajectory(
    obj, align=None, smooth=None, wrap=None, mapping=None, frame=None, map=None
):
    """
    Return an iterator over the trajectory of frames of this view.

    align:
        Pass in a selection string to select atoms against which
        every frame will be aligned. These atoms will be moved
        to the center of the periodic box (if a periodic box
        is used). If "True" is passed, then this will attempt
        to align *ALL* of the coordinates in the view.

        You can also choose to pass in a molecular container,
        and it will align against the atoms in that container,
        assuming they are contained in this view. If not, then
        you need to supply a mapping that maps from the
        atoms in the align container, to the atoms in this view.

    frame:
        The frame of the trajectory against which the alignment
        should be based. For example, `frame=3` would align based
        on the coordinates of the aligned atoms in frame 3 of
        the trajectory. If this is `None` (the default) then the
        first frame will be used.

    mapping: AtomMapping
        An AtomMapping object that maps from atoms in the alignment
        container to atoms in this view. You only need to supply
        this if all of the alignment atoms are not contained
        in this view.

    smooth:
      Pass in the number of frames to smooth (average) the view
      over. If 'True' is passed, then the recommended number
      of frames will be averaged over

    wrap: bool
      Whether or not to wrap the coordinates into the periodic box

    """
    from ._trajectory import TrajectoryIterator

    return TrajectoryIterator(
        obj,
        align=align,
        frame=frame,
        smooth=smooth,
        wrap=wrap,
        mapping=mapping,
        map=map,
    )


MoleculeView.trajectory = _trajectory
SelectorM_Atom_.trajectory = _trajectory
SelectorM_Residue_.trajectory = _trajectory
SelectorM_Chain_.trajectory = _trajectory
SelectorM_Segment_.trajectory = _trajectory
SelectorM_CutGroup_.trajectory = _trajectory
SelectorMol.trajectory = _trajectory


def _cursorsm(obj, map=None):
    from ._cursor import CursorsM

    return CursorsM(parent=obj, map=map)


SelectorM_Atom_.cursor = _cursorsm
SelectorM_Residue_.cursor = _cursorsm
SelectorM_Chain_.cursor = _cursorsm
SelectorM_Segment_.cursor = _cursorsm
SelectorM_CutGroup_.cursor = _cursorsm
SelectorMol.cursor = _cursorsm


def _to_molecules(obj):
    if obj is None:
        return None

    if hasattr(obj, "to_molecule_group"):
        return obj.to_molecule_group().molecules()
    else:
        from ..legacy.Mol import Molecules

        m = Molecules(obj)
        return m


def _energy(obj, other=None, map=None):
    """
    Calculate the total energy of the molecule view(s) in this
    collection.

    other: (view or views, optional)
        An optional second view (or collection of views).
        If this is passed, then the energy between the
        views in this collections and others will be calculated.

    map: (dictionary, optional)
        An optional property map that will be used to find
        or map the properties used for this energy calculation.

    Returns
    -------

    sire.units.GeneralUnit (energy per quantity):
        Returns an energy, with attached components for the
        sub-components (if any) for this energy.
    """
    from ..system import calculate_energy

    if map is None:
        if other is None:
            return calculate_energy(obj)
        else:
            return calculate_energy(obj, _to_molecules(other))
    elif other is None:
        return calculate_energy(obj, map=map)
    else:
        return calculate_energy(obj, _to_molecules(other), map=map)


def _energies(obj, other=None, map=None):
    """
    Calculate the total energy of the molecule view(s) in this
    collection.

    other: (view or views, optional)
        An optional second view (or collection of views).
        If this is passed, then the energy between the
        views in this collections and others will be calculated.

    map: (dictionary, optional)
        An optional property map that will be used to find
        or map the properties used for this energy calculation.

    Returns
    -------

    sire.units.GeneralUnit (energy per quantity):
        Returns an energy, with attached components for the
        sub-components (if any) for this energy.
    """
    return obj.apply("energy", other=other, map=map)


def _atom_energy(obj, other=None, map=None):
    """
    Calculate the total energy of the molecule view(s) in this
    collection.

    other: (view or views, optional)
        An optional second view (or collection of views).
        If this is passed, then the energy between the
        views in this collections and others will be calculated.

    map: (dictionary, optional)
        An optional property map that will be used to find
        or map the properties used for this energy calculation.

    Returns
    -------

    sire.units.GeneralUnit (energy per quantity):
        Returns an energy, with attached components for the
        sub-components (if any) for this energy.
    """
    if other is None:
        from ..units import GeneralUnit

        return GeneralUnit(0)
    elif map is None:
        from ..system import calculate_energy

        return calculate_energy(obj, _to_molecules(other))
    else:
        from ..system import calculate_energy

        return calculate_energy(obj, _to_molecules(other), map=map)


def _total_energy(obj, other=None, map=None):
    """
    Calculate the total energy of the molecule view(s) in this
    collection.

    other: (view or views, optional)
        An optional second view (or collection of views).
        If this is passed, then the energy between the
        views in this collections and others will be calculated.

    map: (dictionary, optional)
        An optional property map that will be used to find
        or map the properties used for this energy calculation.

    Returns
    -------

    sire.units.GeneralUnit (energy per quantity):
        Returns an energy, with attached components for the
        sub-components (if any) for this energy.
    """
    if hasattr(obj, "to_molecule_group"):
        mols = obj.to_molecule_group()
    else:
        from ..legacy.Mol import MoleculeGroup

        mols = MoleculeGroup("all")
        mols.add(obj)

    from ..system import calculate_energy

    if map is None:
        if other is None:
            return calculate_energy(mols.molecules())
        else:
            return calculate_energy(mols.molecules(), _to_molecules(other))
    elif other is None:
        return calculate_energy(mols.molecules(), map=map)
    else:
        return calculate_energy(
            mols.molecules(), _to_molecules(other), map=map
        )


Atom.energy = _atom_energy
Residue.energy = _energy
Chain.energy = _energy
Segment.energy = _energy
CutGroup.energy = _energy
Molecule.energy = _energy

SelectorMol.energy = _total_energy
Selector_Atom_.energy = _total_energy
Selector_Residue_.energy = _total_energy
Selector_Chain_.energy = _total_energy
Selector_Segment_.energy = _total_energy
Selector_CutGroup_.energy = _total_energy

SelectorM_Atom_.energy = _total_energy
SelectorM_Residue_.energy = _total_energy
SelectorM_Chain_.energy = _total_energy
SelectorM_Segment_.energy = _total_energy
SelectorM_CutGroup_.energy = _total_energy

SelectorMol.energies = _energies
Selector_Atom_.energies = _energies
Selector_Residue_.energies = _energies
Selector_Chain_.energies = _energies
Selector_Segment_.energies = _energies
Selector_CutGroup_.energies = _energies

SelectorM_Atom_.energies = _energies
SelectorM_Residue_.energies = _energies
SelectorM_Chain_.energies = _energies
SelectorM_Segment_.energies = _energies
SelectorM_CutGroup_.energies = _energies

MoleculeView.view = _viewfunc
MoleculeView.smiles = _to_smiles
MoleculeView.smarts = _to_smarts
MoleculeView.view2d = _view2d
SelectorMol.view = _viewfunc
SelectorMol.view2d = _selector_view2d
SelectorMol.smiles = _selector_to_smiles
SelectorMol.smarts = _selector_to_smarts
Selector_Atom_.view = _viewfunc
Selector_Residue_.view = _viewfunc
Selector_Chain_.view = _viewfunc
Selector_Segment_.view = _viewfunc
Selector_CutGroup_.view = _viewfunc
SelectorM_Atom_.view = _viewfunc
SelectorM_Atom_.view2d = _selector_view2d
SelectorM_Atom_.smiles = _selector_to_smiles
SelectorM_Atom_.smiles = _selector_to_smarts
SelectorM_Residue_.view = _viewfunc
SelectorM_Residue_.view2d = _selector_view2d
SelectorM_Residue_.smiles = _selector_to_smiles
SelectorM_Residue_.smiles = _selector_to_smarts
SelectorM_Chain_.view = _viewfunc
SelectorM_Chain_.view2d = _selector_view2d
SelectorM_Chain_.smiles = _selector_to_smiles
SelectorM_Chain_.smiles = _selector_to_smarts
SelectorM_Segment_.view = _viewfunc
SelectorM_Segment_.view2d = _selector_view2d
SelectorM_Segment_.smiles = _selector_to_smiles
SelectorM_Segment_.smiles = _selector_to_smarts
SelectorM_CutGroup_.view = _viewfunc
SelectorM_CutGroup_.view2d = _selector_view2d
SelectorM_CutGroup_.smiles = _selector_to_smiles
SelectorM_CutGroup_.smiles = _selector_to_smarts

if not hasattr(SelectorMol, "__orig__find__"):

    def __find__(obj, views):
        """
        Find the index(es) of the passed view(s) in this container.
        This returns a single index if a single view is passed,
        or a list of indexes if multiple views are passed
        (in the same order as the passed views).

        This raises an IndexError if any of the views are
        not in this container.
        """
        matches = obj.__orig__find__(views)

        if matches is None:
            raise IndexError("Cannot find the passed view in the container.")

        if views.is_selector():
            if len(matches) != len(views):
                raise IndexError(
                    "Could not find all of the passed views "
                    "in the container. We only found "
                    f"{len(matches)} of the required "
                    f"{len(views)}."
                )

            return matches
        else:
            if len(matches) != 1:
                raise IndexError(
                    "We matched too many indexes? "
                    f"{matches}. Only expected to match one?"
                )

            return matches[0]

    def __fix__find(CLASS):
        CLASS.__orig__find__ = CLASS.find
        CLASS.find = __find__

    for C in [
        Selector_Atom_,
        Selector_Residue_,
        Selector_Chain_,
        Selector_CutGroup_,
        Selector_Segment_,
        SelectorMol,
        SelectorM_Atom_,
        SelectorM_Residue_,
        SelectorM_Chain_,
        SelectorM_CutGroup_,
        SelectorM_Segment_,
    ]:
        __fix__find(C)


if not hasattr(AtomMapping, "__orig_find__"):

    def __mapping_find__(obj, atoms, container, find_all: bool = True):
        from ..system import System

        if System.is_system(container):
            container = container.atoms()

        return obj.__orig_find__(atoms, container, find_all=find_all)

    def __mapping_map__(obj, atoms, container, match_all: bool = True):
        from ..system import System

        if System.is_system(container):
            container = container.atoms()

        return obj.__orig_map__(atoms, container, find_all=match_all)

    AtomMapping.__orig_find__ = AtomMapping.find
    AtomMapping.__orig_map__ = AtomMapping.map

    AtomMapping.find = __mapping_find__
    AtomMapping.map = __mapping_map__


if not hasattr(Molecule, "perturbation"):

    def __molecule_perturbation(mol):
        """
        Return an interface to the perturbable properties of
        this molecule. Note that the molecule must be
        perturbable to call this function
        """
        from ..morph._perturbation import Perturbation

        return Perturbation(mol)

    def __molecule_is_perturbable(mol):
        """
        Return whether or not this molecule is perturbable
        (can be morphed with a lambda coordinate)
        """
        if mol.has_property("is_perturbable"):
            return mol.property("is_perturbable").as_boolean()
        else:
            return False

    Molecule.perturbation = __molecule_perturbation
    Molecule.is_perturbable = __molecule_is_perturbable


# Remove some temporary variables
del C
