__all__ = ["create_positional_restraints"]


def create_positional_restraints(mols, atoms, k=None, r0=None, map=None):
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

    if hasattr(atoms, "atoms"):
        atoms = atoms.atoms()
    else:
        atoms = mols[atoms].atoms()

    mols = mols.atoms()

    restraints = PositionalRestraints()

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

        restraints.add(
            PositionalRestraint(idxs[0], atom.property(coords_prop), ik, ir0)
        )

    return restraints
