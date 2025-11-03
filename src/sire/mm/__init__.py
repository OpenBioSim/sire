__all__ = [
    "AmberBond",
    "AmberAngle",
    "AmberDihPart",
    "AmberDihedral",
    "Angle",
    "AngleRestraint",
    "AngleRestraints",
    "Bond",
    "BondRestraint",
    "BondRestraints",
    "InverseBondRestraint",
    "InverseBondRestraints",
    "CMAPParameter",
    "Dihedral",
    "Improper",
    "PositionalRestraint",
    "PositionalRestraints",
    "DihedralRestraint",
    "DihedralRestraints",
    "MorsePotentialRestraint",
    "MorsePotentialRestraints",
    "RMSDRestraint",
    "RMSDRestraints",
    "SelectorAngle",
    "SelectorBond",
    "SelectorDihedral",
    "SelectorImproper",
    "SelectorMAngle",
    "SelectorMBond",
    "SelectorMDihedral",
    "SelectorMImproper",
    "LJParameter",
    "LJ1264Parameter",
]

from ..legacy import MM as _MM

from .. import use_new_api as _use_new_api

_use_new_api()

AngleRestraint = _MM.AngleRestraint
AngleRestraints = _MM.AngleRestraints

# It would be better if these were called "DistanceRestraints",
# but there is already a legacy Sire.MM class with this name
BondRestraint = _MM.BondRestraint
BondRestraints = _MM.BondRestraints

InverseBondRestraint = _MM.InverseBondRestraint
InverseBondRestraints = _MM.InverseBondRestraints

BoreschRestraint = _MM.BoreschRestraint
BoreschRestraints = _MM.BoreschRestraints

CMAPParameter = _MM.CMAPParameter

PositionalRestraint = _MM.PositionalRestraint
PositionalRestraints = _MM.PositionalRestraints

RMSDRestraint = _MM.RMSDRestraint
RMSDRestraints = _MM.RMSDRestraints

DihedralRestraint = _MM.DihedralRestraint
DihedralRestraints = _MM.DihedralRestraints

MorsePotentialRestraint = _MM.MorsePotentialRestraint
MorsePotentialRestraints = _MM.MorsePotentialRestraints

AmberBond = _MM.AmberBond
AmberAngle = _MM.AmberAngle
AmberDihPart = _MM.AmberDihPart
AmberDihedral = _MM.AmberDihedral

LJParameter = _MM.LJParameter
LJ1264Parameter = _MM.LJ1264Parameter

Bond = _MM.Bond
SelectorBond = _MM.SelectorBond
SelectorMBond = _MM.SelectorMBond

Angle = _MM.Angle
SelectorAngle = _MM.SelectorAngle
SelectorMAngle = _MM.SelectorMAngle

Dihedral = _MM.Dihedral
SelectorDihedral = _MM.SelectorDihedral
SelectorMDihedral = _MM.SelectorMDihedral

Improper = _MM.Improper
SelectorImproper = _MM.SelectorImproper
SelectorMImproper = _MM.SelectorMImproper

try:
    Bond.__len__ = Bond.nAtoms
except AttributeError:
    Bond.__len__ = Bond.num_atoms

try:
    Angle.__len__ = Angle.nAtoms
except AttributeError:
    Angle.__len__ = Angle.num_atoms

try:
    Dihedral.__len__ = Dihedral.nAtoms
except AttributeError:
    Dihedral.__len__ = Dihedral.num_atoms

try:
    Improper.__len__ = Improper.nAtoms
except AttributeError:
    Improper.__len__ = Improper.num_atoms


_have_fixed_siremm = False


def _fix_siremm():
    global _have_fixed_siremm
    if _have_fixed_siremm:
        return

    from ..mol import (
        __fix_getitem,
        _add_evals,
        _add_property_func,
        _add_apply_func,
        _cursor,
        _cursors,
        _cursorsm,
        _dynamics,
        _minimisation,
        _selector_to_smiles,
        _selector_to_smarts,
        _selector_view2d,
        _trajectory,
        _viewfunc,
    )

    _have_fixed_siremm = True

    for C in [
        Bond,
        SelectorBond,
        SelectorMBond,
        Angle,
        SelectorAngle,
        SelectorMAngle,
        Dihedral,
        SelectorDihedral,
        SelectorMDihedral,
        Improper,
        SelectorImproper,
        SelectorMImproper,
    ]:
        __fix_getitem(C)
        _add_evals(C)
        _add_property_func(C)
        _add_apply_func(C)

    Bond.cursor = _cursor
    SelectorBond.cursor = _cursors
    Angle.cursor = _cursor
    SelectorAngle.cursor = _cursors
    Dihedral.cursor = _cursor
    SelectorDihedral.cursor = _cursors
    Improper.cursor = _cursor
    SelectorImproper.cursor = _cursors

    SelectorMBond.cursor = _cursorsm
    SelectorMAngle.cursor = _cursorsm
    SelectorMDihedral.cursor = _cursorsm
    SelectorMImproper.cursor = _cursorsm

    SelectorMBond.dynamics = _dynamics
    SelectorMAngle.dynamics = _dynamics
    SelectorMDihedral.dynamics = _dynamics
    SelectorMImproper.dynamics = _dynamics

    SelectorMBond.minimisation = _minimisation
    SelectorMAngle.minimisation = _minimisation
    SelectorMDihedral.minimisation = _minimisation
    SelectorMImproper.minimisation = _minimisation

    SelectorMBond.smiles = _selector_to_smiles
    SelectorMAngle.smiles = _selector_to_smiles
    SelectorMDihedral.smiles = _selector_to_smiles
    SelectorMImproper.smiles = _selector_to_smiles
    SelectorMBond.smarts = _selector_to_smarts
    SelectorMAngle.smarts = _selector_to_smarts
    SelectorMDihedral.smarts = _selector_to_smarts
    SelectorMImproper.smarts = _selector_to_smarts

    SelectorMBond.trajectory = _trajectory
    SelectorMAngle.trajectory = _trajectory
    SelectorMDihedral.trajectory = _trajectory
    SelectorMImproper.trajectory = _trajectory

    SelectorBond.view = _viewfunc
    SelectorAngle.view = _viewfunc
    SelectorDihedral.view = _viewfunc
    SelectorImproper.view = _viewfunc

    SelectorMBond.view = _viewfunc
    SelectorMAngle.view = _viewfunc
    SelectorMDihedral.view = _viewfunc
    SelectorMImproper.view = _viewfunc

    SelectorMBond.view2d = _selector_view2d
    SelectorMAngle.view2d = _selector_view2d
    SelectorMDihedral.view2d = _selector_view2d
    SelectorMImproper.view2d = _selector_view2d


try:
    _fix_siremm()
except ImportError:
    pass
