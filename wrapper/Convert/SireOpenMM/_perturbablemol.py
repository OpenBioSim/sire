__all__ = [
    "_changed_atoms",
    "_changed_bonds",
    "_changed_angles",
    "_changed_torsions",
    "_changed_exceptions",
]


def _changed_atoms(obj, to_pandas: bool = True):
    """
    Return a list of the atoms that change parameters in this
    perturbation

    Parameters
    ----------

    to_pandas: bool, optional, default=True
        If True then the list of atoms will be returned as a pandas
        DataFrame
    """
    changed_atoms = []

    for atom, q0, s0, e0, a0, k0, q1, s1, e1, a1, k1 in zip(
        obj.atoms(),
        obj.get_charges0(),
        obj.get_sigmas0(),
        obj.get_epsilons0(),
        obj.get_alphas0(),
        obj.get_kappas0(),
        obj.get_charges1(),
        obj.get_sigmas1(),
        obj.get_epsilons1(),
        obj.get_alphas1(),
        obj.get_kappas1(),
    ):
        if q0 != q1 or s0 != s1 or e0 != e1 or a0 != a1 or k0 != k1:
            changed_atoms.append((atom, q0, s0, e0, a0, k0, q1, s1, e1, a1, k1))

    return changed_atoms


def _changed_bonds(obj, to_pandas: bool = True):
    """
    Return a list of the bonds that change parameters in this
    perturbation

    Parameters
    ----------

    to_pandas: bool, optional, default=True
        If True then the list of bonds will be returned as a pandas
        DataFrame
    """
    changed_bonds = []

    for bond, r0, k0, r1, k1 in zip(
        obj.bonds(),
        obj.get_bond_lengths0(),
        obj.get_bond_ks0(),
        obj.get_bond_lengths1(),
        obj.get_bond_ks1(),
    ):
        if r0 != r1 or k0 != k1:
            changed_bonds.append((bond, r0, k0, r1, k1))

    return changed_bonds


def _changed_angles(obj, to_pandas: bool = True):
    """
    Return a list of the angles that change parameters in this
    perturbation

    Parameters
    ----------

    to_pandas: bool, optional, default=True
        If True then the list of angles will be returned as a pandas
        DataFrame
    """
    changed_angles = []

    for angle, theta0, k0, theta1, k1 in zip(
        obj.angles(),
        obj.get_angle_sizes0(),
        obj.get_angle_ks0(),
        obj.get_angle_sizes1(),
        obj.get_angle_ks1(),
    ):
        if theta0 != theta1 or k0 != k1:
            changed_angles.append((angle, theta0, k0, theta1, k1))

    return changed_angles


def _changed_torsions(obj, to_pandas: bool = True):
    """
    Return a list of the torsions that change parameters in this
    perturbation. Note that this combines the dihedrals and improper
    dihedrals into a single list.

    Parameters
    ----------

    to_pandas: bool, optional, default=True
        If True then the list of torsions will be returned as a pandas
        DataFrame
    """
    changed_torsions = []

    for torsion, k0, p0, ph0, k1, p1, ph1 in zip(
        obj.torsions(),
        obj.get_torsion_ks0(),
        obj.get_torsion_periodicities0(),
        obj.get_torsion_phases0(),
        obj.get_torsion_ks1(),
        obj.get_torsion_periodicities1(),
        obj.get_torsion_phases1(),
    ):
        if k0 != k1 or p0 != p1 or ph0 != ph1:
            changed_torsions.append((torsion, k0, p0, ph0, k1, p1, ph1))

    return changed_torsions


def _changed_exceptions(obj, to_pandas: bool = True):
    """
    Return a list of the exceptions that change parameters in this
    perturbation

    Parameters
    ----------

    to_pandas: bool, optional, default=True
        If True then the list of exceptions will be returned as a pandas
        DataFrame
    """
    changed_exceptions = []

    return changed_exceptions
