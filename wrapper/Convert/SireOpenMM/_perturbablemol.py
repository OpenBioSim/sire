__all__ = [
    "_changed_atoms",
    "_changed_bonds",
    "_changed_angles",
    "_changed_torsions",
    "_changed_exceptions",
    "_changed_constraints",
    "_get_lever_values",
]


def _get_lever_values(
    obj,
    schedule=None,
    lambda_values=None,
    num_lambda: int = 101,
    to_pandas: bool = True,
):
    """
    Return the value of all of the parameters for this perturbable molecule
    at all of the specified values of lambda, given the passed
    lambda schedule. If no schedule is passed then a default morph
    will be used.
    """
    if lambda_values is None:
        import numpy as np

        lambda_values = np.linspace(0.0, 1.0, num_lambda)

    lambda_values = [float(x) for x in lambda_values]

    if schedule is None:
        from ...cas import LambdaSchedule

        schedule = LambdaSchedule.standard_morph()

    from . import LambdaLever

    lever = LambdaLever()
    lever.set_schedule(schedule)

    results = lever.get_lever_values(lambda_values=lambda_values, mol=obj)

    if to_pandas:
        import pandas as pd

        colnames = results[0]
        columns = results
        columns.pop_front()

        results = {}

        results["index"] = list(range(len(columns[0])))

        for i in range(len(colnames)):
            results[colnames[i]] = [x for x in columns[i]]

        results = pd.DataFrame(results)

    return results


def _name(atom):
    return f"{atom.name().value()}:{atom.number().value()}"


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
            if abs(q0) <= 1e-9:
                q0 = 0.0
            if abs(q1) <= 1e-9:
                q1 = 0.0
            if abs(s0) <= 1e-9:
                s0 = 0.0
            if abs(s1) <= 1e-9:
                s1 = 0.0
            if abs(e0) <= 1e-9:
                e0 = 0.0
            if abs(e1) <= 1e-9:
                e1 = 0.0
            if abs(a0) <= 1e-9:
                a0 = 0.0
            if abs(a1) <= 1e-9:
                a1 = 0.0
            if abs(k0) <= 1e-9:
                k0 = 0.0
            if abs(k1) <= 1e-9:
                k1 = 0.0

            if to_pandas:
                atom = _name(atom)

            changed_atoms.append((atom, q0, q1, s0, s1, e0, e1, a0, a1, k0, k1))

    if to_pandas:
        import pandas as pd

        changed_atoms = pd.DataFrame(
            changed_atoms,
            columns=[
                "atom",
                "charge0",
                "charge1",
                "sigma0",
                "sigma1",
                "epsilon0",
                "epsilon1",
                "alpha0",
                "alpha1",
                "kappa0",
                "kappa1",
            ],
        )

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
            if to_pandas:
                bond = f"{_name(bond[0])}-{_name(bond[1])}"

            changed_bonds.append((bond, r0, r1, k0, k1))

    if to_pandas:
        import pandas as pd

        changed_bonds = pd.DataFrame(
            changed_bonds, columns=["bond", "length0", "length1", "k0", "k1"]
        )

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
            if to_pandas:
                angle = f"{_name(angle[0])}-{_name(angle[1])}-{_name(angle[2])}"

            changed_angles.append((angle, theta0, theta1, k0, k1))

    if to_pandas:
        import pandas as pd

        changed_angles = pd.DataFrame(
            changed_angles, columns=["angle", "size0", "size1", "k0", "k1"]
        )

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
            if to_pandas:
                torsion = f"{_name(torsion[0])}-{_name(torsion[1])}-{_name(torsion[2])}-{_name(torsion[3])}"

            changed_torsions.append((torsion, k0, k1, p0, p1, ph0, ph1))

    if to_pandas:
        import pandas as pd

        changed_torsions = pd.DataFrame(
            changed_torsions,
            columns=[
                "torsion",
                "k0",
                "k1",
                "periodicity0",
                "periodicity1",
                "phase0",
                "phase1",
            ],
        )

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

    atoms = obj.atoms()

    for (atom0, atom1), q0, q1, lj0, lj1 in zip(
        obj.get_exception_atoms(),
        obj.get_charge_scales0(),
        obj.get_charge_scales1(),
        obj.get_lj_scales0(),
        obj.get_lj_scales1(),
    ):
        if q0 != q1 or lj0 != lj1:
            atom0 = atoms[atom0]
            atom1 = atoms[atom1]

            atompair = (atom0, atom1)

            if to_pandas:
                atompair = f"{_name(atom0)}-{_name(atom1)}"

            changed_exceptions.append((atompair, q0, q1, lj0, lj1))

    if to_pandas:
        import pandas as pd

        changed_exceptions = pd.DataFrame(
            changed_exceptions,
            columns=[
                "atompair",
                "charge_scale0",
                "charge_scale1",
                "lj_scale0",
                "lj_scale1",
            ],
        )

    return changed_exceptions


def _changed_constraints(obj, to_pandas: bool = True):
    """
    Return a list of the constraints that change parameters in this
    perturbation

    Parameters
    ----------

    to_pandas: bool, optional, default=True
        If True then the list of constraints will be returned as a pandas
        DataFrame
    """
    changed_constraints = []

    atoms = obj.atoms()

    for atom0, atom1, r0, r1 in obj.get_perturbable_constraints_with_atoms():
        if r0 != r1:
            atom0 = atoms[atom0]
            atom1 = atoms[atom1]

            atompair = (atom0, atom1)

            if to_pandas:
                atompair = f"{_name(atom0)}-{_name(atom1)}"

            changed_constraints.append((atompair, r0, r1))

    if to_pandas:
        import pandas as pd

        changed_constraints = pd.DataFrame(
            changed_constraints,
            columns=["atompair", "length0", "length1"],
        )

    return changed_constraints
