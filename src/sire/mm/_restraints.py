__all__ = [
    "create_bond_restraints",
    "create_boresch_restraints",
    "create_positional_restraints",
]


def _to_atoms(mols, atoms):
    """
    Internal function used to convert `mols[atoms]` into a list
    of atoms
    """
    if hasattr(atoms, "atoms"):
        return atoms.atoms()

    from ..legacy.Base import NumberProperty, IntegerArrayProperty

    if type(atoms) is NumberProperty or type(atoms) is IntegerArrayProperty:
        return mols[atoms.value()].atoms()
    else:
        return mols[atoms].atoms()


def create_boresch_restraints(
    mols,
    receptor,
    ligand,
    kr,
    r0,
    ktheta,
    theta0,
    kphi,
    phi0,
    name: str = None,
    map=None,
):
    """
    Create a set of Boresch restraints that will hold the passed
    ligand in a its relative binding mode relative to the
    passed receptor. All of the atoms in both 'ligand' and
    'receptor' must be contained in 'mols'.

    The BoreschRestraint will be a set of six restraints between
    three identified ligand atoms, and three identified receptor
    atoms.

    1. A single distance restraint, with specified kr and r0 parameters
    2. Two angle restraints, with the specified two ktheta and theta0
       parameters
    3. Three torsion restraints, with the specified three kphi and phi0
       parameters

    This will create a single BoreschRestraint, which will be passed
    back in a BoreschRestraints object.
    """
    from . import BoreschRestraint, BoreschRestraints
    from .. import u
    from ..base import create_map

    map = create_map(map)

    receptor = _to_atoms(mols, receptor)
    ligand = _to_atoms(mols, ligand)

    if len(receptor) != 3 or len(ligand) != 3:
        # Eventually will choose the best atoms from the receptor
        # and ligand...
        raise ValueError(
            "You need to provide 3 receptor atoms and 3 ligand atoms"
        )

    mols = mols.atoms()

    b = BoreschRestraint(
        receptor=mols.find(receptor),
        ligand=mols.find(ligand),
        r0=r0,
        theta0=theta0,
        phi0=phi0,
        kr=kr,
        ktheta=ktheta,
        kphi=kphi,
    )

    if name is None:
        return BoreschRestraints(b)
    else:
        return BoreschRestraint(name, b)


def create_bond_restraints(
    mols, atoms0, atoms1, k, r0, name: str = None, map=None
):
    """
    Create a set of bond restraints from all of the atoms in 'atoms0'
    to all of the atoms in 'atoms1' where all atoms are
    contained in the container 'mols', using the
    passed values of the force constant 'k' and equilibrium
    bond length r0.

    These restraints will be per bond. If a list of k and/or r0
    values are passed, then different values could be used for
    different bonds (assuming the same number as the number of
    bodns). Otherwise, all bonds will use the same parameters.
    """
    from . import BondRestraint, BondRestraints
    from .. import u
    from ..base import create_map

    map = create_map(map)

    if type(k) is list:
        k = [u(x) for x in k]
    else:
        k = [u(k)]

    if type(r0) is list:
        r0 = [u(x) for x in r0]
    else:
        r0 = [u(r0)]

    atoms0 = _to_atoms(mols, atoms0)
    atoms1 = _to_atoms(mols, atoms1)

    if len(atoms0) != len(atoms1):
        raise ValueError(
            "Cannot create the set of bond restraints because the number "
            f"of atoms in the first group ({len(atoms0)}) is not equal to "
            f"the number of atoms in the second ({len(atoms1)})."
        )

    mols = mols.atoms()

    if name is None:
        restraints = BondRestraints()
    else:
        restraints = BondRestraints(name=name)

    for i, (atom0, atom1) in enumerate(zip(atoms0, atoms1)):
        idxs0 = mols.find(atom0)
        idxs1 = mols.find(atom1)

        if type(idxs0) is int:
            idxs0 = [idxs0]

        if type(idxs1) is int:
            idxs1 = [idxs1]

        if len(idxs0) == 0:
            raise KeyError(
                f"Could not find atom {atom0} in the molecules. Please ensure "
                "that 'mols' contains all of that atoms, or else we can't "
                "add the positional restraints."
            )

        if len(idxs1) == 0:
            raise KeyError(
                f"Could not find atom {atom1} in the molecules. Please ensure "
                "that 'mols' contains all of that atoms, or else we can't "
                "add the positional restraints."
            )

        if i < len(k):
            ik = k[i]
        else:
            ik = k[-1]

        if i < len(r0):
            ir0 = r0[i]
        else:
            ir0 = r0[-1]

        restraints.add(BondRestraint(idxs0[0], idxs1[0], ik, ir0))

    return restraints


def create_positional_restraints(
    mols, atoms, k=None, r0=None, position=None, name: str = None, map=None
):
    """
    Create a set of position restraints for the atoms specified in
    'atoms' that are contained in the container 'mols', using the
    passed values of the force constant 'k' and flat-bottom potential
    well-width 'r0' for the restraints.

    These restraints will be per atom. If a list of k and/or r0
    values are passed, then different values could be used for
    different atoms (assuming the same number as the number of
    atoms). Otherwise, all atoms will use the same parameters.

    If 'r0' is not specified, then a simple harmonic restraint
    is used.
    """
    from . import PositionalRestraint, PositionalRestraints
    from .. import u
    from ..base import create_map

    map = create_map(map)

    if k is None:
        k = [u("150 kcal mol-1 A-2")]
    elif type(k) is list:
        k = [u(x) for x in k]
    else:
        k = [u(k)]

    if r0 is None:
        r0 = [u("0")]
    elif type(r0) is list:
        r0 = [u(x) for x in r0]
    else:
        r0 = [u(r0)]

    atoms = _to_atoms(mols, atoms)

    mols = mols.atoms()

    if name is None:
        restraints = PositionalRestraints()
    else:
        restraints = PositionalRestraints(name=name)

    coords_prop = map["coordinates"]

    for i, atom in enumerate(atoms.atoms()):
        idxs = mols.find(atom)

        if type(idxs) is int:
            idxs = [idxs]

        elif len(idxs) == 0:
            raise KeyError(
                f"Could not find atom {atom} in the molecules. Please ensure "
                "that 'mols' contains all of that atoms, or else we can't "
                "add the positional restraints."
            )

        if i < len(k):
            ik = k[i]
        else:
            ik = k[-1]

        if i < len(r0):
            ir0 = r0[i]
        else:
            ir0 = r0[-1]

        if position is None:
            restraints.add(
                PositionalRestraint(
                    idxs[0], atom.property(coords_prop), ik, ir0
                )
            )
        else:
            restraints.add(PositionalRestraint(idxs[0], position, ik, ir0))

    return restraints
