__all__ = [
    "boresch",
    "bond",
    "distance",
    "positional",
]

from .. import u


def _to_atoms(mols, atoms):
    """
    Internal function used to convert `mols[atoms]` into a list
    of atoms
    """
    from ..mol import selection_to_atoms

    return selection_to_atoms(mols, atoms)


def boresch(
    mols,
    receptor,
    ligand,
    kr=None,
    ktheta=None,
    kphi=None,
    r0=None,
    theta0=None,
    phi0=None,
    name=None,
    map=None,
    temperature=u("298 K"),
):
    """
    Create a set of Boresch restraints that will restrain the 6
    external degrees of freedom of the ligand relative to the receptor.
    All of the atoms in both 'ligand' and 'receptor' must be contained in
    'mols'.

    The BoreschRestraint will be a set of six restraints between
    three identified ligand atoms, and three identified receptor
    atoms:

    1. A single distance restraint, with specified force constant (kr)
       and equilibrium distance (r0) parameters.
    2. Two angle restraints, with specified force constants (ktheta)
       and equilibrium angles (theta0) parameters.
    3. Three torsion restraints, with specified force constants (kphi)
         and equilibrium angles (phi0) parameters.

    This will create a single BoreschRestraint, which will be passed
    back in a BoreschRestraints object.

    Parameters
    ----------
    mols : sire.system._system.System
        The system containing the receptor and ligand.

    receptor : SireMol::Selector<SireMol::Atom>
        The receptor atoms to restrain.

    ligand : SireMol::Selector<SireMol::Atom>
        The ligand atoms to restrain.

    kr : str or SireUnits::Dimension::GeneralUnit, optional
        The force constant for the distance restraint. If None, this will
        default to 10 kcal mol-1 A-2. Default is None.

    ktheta : str or SireUnits::Dimension::GeneralUnit or list of str or SireUnits::Dimension::GeneralUnit, optional
        The force constants for the angle restraints, in the order kthetaA,
        kthetaB If None, this will default to 100 kcal mol-1 rad-2 for
        both angle restraints.  If a list, then this should be a list of
        length 2 containing the force constants for the two angle
        restraints. If a single value, then this will be used for both
        angle restraints. Default is None.

    kphi : str or SireUnits::Dimension::GeneralUnit or list of str or SireUnits::Dimension::GeneralUnit, optional
        The force constants for the torsion restraints, in the order kthetaA,
        kthetaB, kthetaC. If None, this will default to 100 kcal mol-1 rad-2
        for all three torsion restraints.  If a list, then this should be a
        list of length 3 containing the force constants for the three
        torsion restraints. If a single value, then this will be used for
        all three torsion restraints. Default is None.

    r0 : str or SireUnits::Dimension::GeneralUnit, optional
        The equilibrium distance for the distance restraint. If None, this
        will be measured from the current coordinates of the atoms. Default
        is None.

    theta0 : list of str or SireUnits::Dimension::GeneralUnit, optional
        The equilibrium angles for the angle restraints. If None, these
        will be measured from the current coordinates of the atoms. If a
        list, then this should be a list of length 2 containing the
        equilibrium angles for the two angle restraints. Default is None.

    phi0 : list of str or SireUnits::Dimension::GeneralUnit, optional
        The equilibrium angles for the torsion restraints. If None, these
        will be measured from the current coordinates of the atoms. If a
        list, then this should be a list of length 3 containing the
        equilibrium angles for the three torsion restraints. Default is None.

    name : str, optional
        The name of the restraint. If None, then a default name will be
        used. Default is None.

    map : dict, optional
        A property map. Default is None.

    temperature : str or SireUnits::Dimension::GeneralUnit, optional
        The temperature to use when checking for unstable restraints. If
        None, then this will default to 298 K. Default is None.

    Returns
    -------
    BoreschRestraints : SireMM::BoreschRestraints
        The Boresch restraints.

    Examples
    --------
    Create a set of Boresch restraints for the ligand in the system
    'system', specifying all of the force constants and equilibrium
    values:
    ```
    my_boresch_restraint = boresch(
        system,
        receptor=system["protein"]["atomidx 1538 or atomidx 1518 or atomidx 1540"],
        ligand=system["resname LIG"]["atomidx 4 or atomidx 3 or atomidx 5"],
        kr="6.2012 kcal mol-1 A-2",
        ktheta=["28.7685 kcal mol-1 rad-2", "24.8204 kcal mol-1 rad-2"],
        kphi=["59.8626 kcal mol-1 rad-2", "0.7923 kcal mol-1 rad-2", "55.1775 kcal mol-1 rad-2"],
        r0="7.687 A",
        theta0=["1.3031 rad", "1.4777 rad"],
        phi0=["2.5569 rad", "2.9359 rad", "1.4147 rad"],
    )
    ```
    """
    from ..base import create_map
    from ..mm import BoreschRestraint, BoreschRestraints

    map = create_map(map)

    receptor = _to_atoms(mols, receptor)
    ligand = _to_atoms(mols, ligand)

    if len(receptor) != 3 or len(ligand) != 3:
        # Eventually will choose the best atoms from the receptor
        # and ligand...
        raise ValueError("You need to provide 3 receptor atoms and 3 ligand atoms")

    from .. import measure

    default_distance_k = u("10 kcal mol-1 A-2")
    default_angle_k = u("100 kcal mol-1 rad-2")

    # Use the user-specified equilibrium values if they are provided.
    distance = [[ligand[0], receptor[0]]]
    angles = [
        [receptor[1], receptor[0], ligand[0]],
        [receptor[0], ligand[0], ligand[1]],
    ]
    dihedrals = [
        [receptor[2], receptor[1], receptor[0], ligand[0]],
        [receptor[1], receptor[0], ligand[0], ligand[1]],
        [receptor[0], ligand[0], ligand[1], ligand[2]],
    ]

    restraint_components = {
        "distance": {
            "input_k": kr,
            "validated_k": [],
            "input_equil": r0,
            "measure": distance,
            "validated_equil": [],
        },
        "angle": {
            "input_k": ktheta,
            "validated_k": [],
            "input_equil": theta0,
            "measure": angles,
            "validated_equil": [],
        },
        "dihedral": {
            "input_k": kphi,
            "validated_k": [],
            "input_equil": phi0,
            "measure": dihedrals,
            "validated_equil": [],
        },
    }

    for restraint_component in restraint_components:
        n_measures = len(restraint_components[restraint_component]["measure"])

        # Get the force constants.
        if restraint_components[restraint_component]["input_k"] is None:
            restraint_components[restraint_component]["validated_k"] = n_measures * [
                default_distance_k
                if restraint_component == "distance"
                else default_angle_k
            ]
        elif type(restraint_components[restraint_component]["input_k"]) is not list:
            # Populate the list with the single specified value.
            restraint_components[restraint_component]["validated_k"] = n_measures * [
                u(restraint_components[restraint_component]["input_k"])
            ]
        else:
            if len(restraint_components[restraint_component]["input_k"]) == 0:
                # Empty list - populate with default values.
                restraint_components[restraint_component][
                    "validated_k"
                ] = n_measures * [
                    default_distance_k
                    if restraint_component == "distance"
                    else default_angle_k
                ]
            elif len(restraint_components[restraint_component]["input_k"]) == 1:
                # List of length 1 - populate with that value.
                restraint_components[restraint_component][
                    "validated_k"
                ] = n_measures * [
                    u(restraint_components[restraint_component]["input_k"][0])
                ]
            elif (
                len(restraint_components[restraint_component]["input_k"]) == n_measures
            ):
                # List of the correct length for this restraint component.
                restraint_components[restraint_component]["validated_k"] = [
                    u(x) for x in restraint_components[restraint_component]["input_k"]
                ]
            else:
                # List of the wrong length.
                raise ValueError(
                    f"Input force constants for {restraint_component} must be a single value or a list of length {n_measures}"
                )

        # Get the equilibrium values.
        if restraint_components[restraint_component]["input_equil"] is None:
            # Set all components to None - these will be measured from the structure later.
            restraint_components[restraint_component]["input_equil"] = [
                None for i in range(n_measures)
            ]

        if type(restraint_components[restraint_component]["input_equil"]) is not list:
            # Only allow this if we are dealing with the distance component, as this is the only one that can be a single value.
            if restraint_component == "distance":
                restraint_components[restraint_component][
                    "input_equil"
                ] = n_measures * [
                    u(restraint_components[restraint_component]["input_equil"])
                ]
            else:
                raise ValueError(
                    f"Input equilibrium values for {restraint_component} must be a list of length {n_measures} of values or Nones"
                )

        elif (
            len(restraint_components[restraint_component]["input_equil"]) != n_measures
        ):
            raise ValueError(
                f"If specified, equilibrium values for {restraint_component} must be a list of length {n_measures} of values or Nones"
            )

        # Now validate the input equilibrium values, replacing Nones with measured values.
        for i, equil in enumerate(
            restraint_components[restraint_component]["input_equil"]
        ):
            if equil is not None:
                restraint_components[restraint_component]["validated_equil"].append(
                    u(equil)
                )
            else:
                restraint_components[restraint_component]["validated_equil"].append(
                    measure(*restraint_components[restraint_component]["measure"][i])
                )

    # Warn the user if the restraint is likely to be unstable.
    _check_stability_boresch_restraint(restraint_components, temperature)

    mols = mols.atoms()

    b = BoreschRestraint(
        receptor=mols.find(receptor),
        ligand=mols.find(ligand),
        r0=restraint_components["distance"]["validated_equil"][0],
        theta0=restraint_components["angle"]["validated_equil"],
        phi0=restraint_components["dihedral"]["validated_equil"],
        kr=restraint_components["distance"]["validated_k"][0],
        ktheta=restraint_components["angle"]["validated_k"],
        kphi=restraint_components["dihedral"]["validated_k"],
    )

    if name is None:
        return BoreschRestraints(b)
    else:
        return BoreschRestraints(name, b)


def _check_stability_boresch_restraint(restraint_components, temperature=u("298 K")):
    """
    Internal function to check for likely unstable Boresch restraints.
    """
    import warnings as _warnings

    from .. import units

    # Check for unstable combinations of force constants.
    if restraint_components["distance"]["validated_k"][0].value() == 0:
        raise ValueError('"kr" cannot be zero')

    if restraint_components["angle"]["validated_k"][0].value() == 0:
        if (
            restraint_components["dihedral"]["validated_k"][0].value() != 0
            or restraint_components["dihedral"]["validated_k"][1].value() != 0
        ):
            raise ValueError(
                "Restraining phiA or phiB without restraining thetaA "
                "will produce unstable Boresch restraints."
            )

    if restraint_components["angle"]["validated_k"][1].value() == 0:
        if (
            restraint_components["dihedral"]["validated_k"][1].value() != 0
            or restraint_components["dihedral"]["validated_k"][2].value() != 0
        ):
            raise ValueError(
                "Restraining phiB or phiC without restraining thetaB "
                "will produce unstable Boresch restraints."
            )

    # Ensure that restrained angles are >= 10 kT from collinear.
    for equil_angle, k_angle in zip(
        restraint_components["angle"]["validated_equil"],
        restraint_components["angle"]["validated_k"],
    ):
        if k_angle.value() != 0:
            # Find the minimum stable angle "distance". We use the squared values as sire units don't support
            # taking the square root.
            min_stable_dist_sq = (2 * 10 * units.k_boltz * temperature) / k_angle
            min_dist_sq = min(
                [abs(equil_angle - 0), abs(equil_angle - 180 * units.degrees)]
            ).pow(2)
            if min_dist_sq < min_stable_dist_sq:
                _warnings.warn(
                    f"The equilibrium value angle value of {equil_angle} is within 10 kT of"
                    "collinearity, which may result in unstable Boresch restraints."
                    " Consider increasing the force constants or selecting equilibrium"
                    " values further from 0 or pi radians."
                )


def distance(mols, atoms0, atoms1, r0=None, k=None, name=None, map=None):
    """
    Create a set of distance restraints from all of the atoms in 'atoms0'
    to all of the atoms in 'atoms1' where all atoms are
    contained in the container 'mols', using the
    passed values of the force constant 'k' and equilibrium
    bond length r0.

    These restraints will be per atom-atom distance. If a list of k and/or r0
    values are passed, then different values could be used for
    different atom-atom distances (assuming the same number as the number of
    atom-atom distances). Otherwise, all atom-atom distances will use the
    same parameters.

    If r0 is None, then the current atom-atom distance for
    each atom-atom pair will be used as the equilibium value.

    If k is None, then a default value of 150 kcal mol-1 A-2 will be used
    """
    from .. import u
    from ..base import create_map
    from ..mm import BondRestraint, BondRestraints

    map = create_map(map)

    if k is None:
        k = [u("150 kcal mol-1 A-2")]
    elif type(k) is list:
        k = [u(x) for x in k]
    else:
        k = [u(k)]

    atoms0 = _to_atoms(mols, atoms0)
    atoms1 = _to_atoms(mols, atoms1)

    if atoms0.is_empty() or atoms1.is_empty():
        raise ValueError("We need at least one atom in each group")

    while len(atoms0) < len(atoms1):
        atoms0 += atoms0[-1]

    while len(atoms1) < len(atoms0):
        atoms1 += atoms1[-1]

    if r0 is None:
        # calculate all of the current distances
        from .. import measure

        r0 = []
        for atom0, atom1 in zip(atoms0, atoms1):
            r0.append(measure(atom0, atom1))
    elif type(r0) is list:
        r0 = [u(x) for x in r0]
    else:
        r0 = [u(r0)]

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


def bond(*args, **kwargs):
    """
    Synonym for distance(), as a bond restraint is treated the same
    as a distance restraint
    """
    return distance(*args, **kwargs)


def positional(mols, atoms, k=None, r0=None, position=None, name=None, map=None):
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

    If 'k' is not specified, then a default of 150 kcal mol-1 A-2
    will be used.
    """
    from .. import u
    from ..base import create_map
    from ..mm import PositionalRestraint, PositionalRestraints

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

    if position is not None:
        from ..maths import Vector

        if type(position) is not list:
            position = len(atoms) * [Vector.to_vector(position)]
        else:
            position = [Vector.to_vector(x) for x in position]

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
                PositionalRestraint(idxs[0], atom.property(coords_prop), ik, ir0)
            )
        else:
            restraints.add(PositionalRestraint(idxs[0], position[i], ik, ir0))

    return restraints
