__all__ = [
    "angle",
    "boresch",
    "bond",
    "dihedral",
    "distance",
    "inverse_bond",
    "inverse_distance",
    "morse_potential",
    "positional",
    "rmsd",
]

from .. import u


def _to_atoms(mols, atoms):
    """
    Internal function used to convert `mols[atoms]` into a list
    of atoms
    """
    from ..mol import selection_to_atoms

    return selection_to_atoms(mols, atoms)


def angle(mols, atoms, theta0=None, ktheta=None, use_pbc=None, name=None, map=None):
    """
    Create a set of angle restraints from all of the atoms in 'atoms'
    where all atoms are contained in the container 'mols', using the
    passed values of the force constant 'ktheta' and equilibrium
    angle value theta0.

    If theta0 is None, then the current angle for
    provided atoms will be used as the equilibium value.

    If ktheta is None, then a default value of 100 kcal mol-1 rad-2 will be used

    Parameters
    ----------
    mols : sire.system._system.System
        The system containing the atoms.

    atoms : SireMol::Selector<SireMol::Atom>
        The atoms to restrain.

    ktheta : str or SireUnits::Dimension::GeneralUnit or, optional
        The force constants for the angle restraints.
        If None, this will default to 100 kcal mol-1 rad-2.
        Default is None.

    theta0 : str or SireUnits::Dimension::GeneralUnit, optional
        The equilibrium angles for the angle restraints. If None, these
        will be measured from the current coordinates of the atoms.
        Default is None.

    use_pbc : bool, optional
        Whether to use periodic boundary conditions when calculating
        the angle. Default is None.

    Returns
    -------
    AngleRestraints : SireMM::AngleRestraints
        A container of angle restraints, where the first restraint is
        the AngleRestraint created. The angle restraint created can be
        extracted with AngleRestraints[0].
    """
    from .. import u
    from ..base import create_map
    from ..mm import AngleRestraint, AngleRestraints

    map = create_map(map)
    map_dict = map.to_dict()
    ktheta = ktheta if ktheta is not None else map_dict.get("ktheta", None)
    theta0 = theta0 if theta0 is not None else map_dict.get("theta0", None)
    use_pbc = use_pbc if use_pbc is not None else map_dict.get("use_pbc", None)
    name = name if name is not None else map_dict.get("name", None)

    if use_pbc is not None:
        if not isinstance(use_pbc, bool):
            raise ValueError("'use_pbc' must be of type 'bool'")
    else:
        use_pbc = False

    atoms = _to_atoms(mols, atoms)

    if len(atoms) != 3:
        raise ValueError(
            "You need to provide 3 atoms to create an angle restraint"
            f"whereas {len(atoms)} atoms were provided."
        )

    from .. import measure

    if ktheta is None:
        ktheta = u("100 kcal mol-1 rad-2")

    elif type(ktheta) is list:
        raise NotImplementedError(
            "Setup of multiple angle restraints simultaneously is not currently supported. "
            "Please set up each restraint individually and then combine them into multiple restraints."
        )

    if theta0 is None:
        from .. import measure

        theta0 = measure(atoms[0], atoms[1], atoms[2])

    elif type(theta0) is list:
        raise NotImplementedError(
            "Setup of multiple angle restraints simultaneously is not currently supported. "
            "Please set up each restraint individually and then combine them into multiple restraints."
        )
    else:
        theta0 = u(theta0)

    mols = mols.atoms()

    if name is None:
        restraints = AngleRestraints()
    else:
        restraints = AngleRestraints(name=name)

    restraints.add(AngleRestraint(mols.find(atoms), theta0, ktheta))

    # Set the use_pbc flag.
    restraints._use_pbc = use_pbc

    return restraints


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
    use_pbc=None,
    name=None,
    map=None,
    temperature=u("298 K"),
):
    """
    Create a set of Boresch restraints that will restrain the 6
    external degrees of freedom of the ligand relative to the receptor.
    All of the atoms in both 'ligand' and 'receptor' must be contained in
    'mols'. Note that restraint energies are defined as k*x**2 (so forces
    are defined as 2*k*x) and hence the 'kr', 'ktheta' and 'kphi' values are
    half the force constants for the distance, angle and torsion restraints.

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
        Half the force constant for the distance restraint. If None, this will
        default to 5 kcal mol-1 A-2. Default is None.

    ktheta : str or SireUnits::Dimension::GeneralUnit or list of str or SireUnits::Dimension::GeneralUnit, optional
        Half the force constants for the angle restraints, in the order kthetaA,
        kthetaB If None, this will default to 50 kcal mol-1 rad-2 for
        both angle restraints.  If a list, then this should be a list of
        length 2 containing the force constants for the two angle
        restraints. If a single value, then this will be used for both
        angle restraints. Default is None.

    kphi : str or SireUnits::Dimension::GeneralUnit or list of str or SireUnits::Dimension::GeneralUnit, optional
        Half the force constants for the torsion restraints, in the order kthetaA,
        kthetaB, kthetaC. If None, this will default to 50 kcal mol-1 rad-2
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

    use_pbc : bool, optional
        Whether to use periodic boundary conditions when calculating
        the distance, angles, and torsions. Default is None.

    name : str, optional
        The name of the restraint. If None, then a default name will be
        used. Default is None.

    map : dict, optional
        A dictionary of additional options. Note that any options
        set in this dictionary that are also specified via one of
        the arguments above will be overridden by the argument
        value

    temperature : str or SireUnits::Dimension::GeneralUnit, optional
        The temperature to use when checking for unstable restraints. If
        None, then this will default to 298 K. Default is None.

    Returns
    -------
    BoreschRestraints : SireMM::BoreschRestraints
        A container of Boresch restraints, where the first restraint is
        the BoreschRestraint created. The Boresch restraint created can be
        extracted with BoreschRestraints[0].

    Examples
    --------
    Create a set of Boresch restraints for the ligand in the system
    'system', specifying all of the force constants and equilibrium
    values:

    >>> my_boresch_restraints = boresch(
    >>>     system,
    >>>     receptor=system["protein"]["atomidx 1538 or atomidx 1518 or atomidx 1540"],
    >>>     ligand=system["resname LIG"]["atomidx 4 or atomidx 3 or atomidx 5"],
    >>>     kr="6.2012 kcal mol-1 A-2",
    >>>     ktheta=["28.7685 kcal mol-1 rad-2", "24.8204 kcal mol-1 rad-2"],
    >>>     kphi=["59.8626 kcal mol-1 rad-2", "0.7923 kcal mol-1 rad-2", "55.1775 kcal mol-1 rad-2"],
    >>>     r0="7.687 A",
    >>>     theta0=["1.3031 rad", "1.4777 rad"],
    >>>     phi0=["2.5569 rad", "2.9359 rad", "1.4147 rad"],
    >>> )
    >>> my_boresch_restraint = my_boresch_restraints[0]
    """
    from ..base import create_map
    from ..mm import BoreschRestraint, BoreschRestraints

    # If an argument is not specified in the function
    # arguments, but is specified in the map, update
    # the argument from the map.
    map = create_map(map)
    map_dict = map.to_dict()
    kr = kr if kr is not None else map_dict.get("kr", None)
    ktheta = ktheta if ktheta is not None else map_dict.get("ktheta", None)
    kphi = kphi if kphi is not None else map_dict.get("kphi", None)
    r0 = r0 if r0 is not None else map_dict.get("r0", None)
    theta0 = theta0 if theta0 is not None else map_dict.get("theta0", None)
    phi0 = phi0 if phi0 is not None else map_dict.get("phi0", None)
    use_pbc = use_pbc if use_pbc is not None else map_dict.get("use_pbc", None)
    name = name if name is not None else map_dict.get("name", None)
    temperature = (
        temperature if temperature is not None else map_dict.get("temperature", None)
    )

    receptor = _to_atoms(mols, receptor)
    ligand = _to_atoms(mols, ligand)

    if len(receptor) != 3 or len(ligand) != 3:
        # Eventually will choose the best atoms from the receptor
        # and ligand...
        raise ValueError(
            "You need to provide 3 receptor atoms and 3 ligand atoms"
            f"but only {len(receptor)} receptor atoms and {len(ligand)} "
            f"ligand atoms were provided."
        )

    if use_pbc is not None:
        if not isinstance(use_pbc, bool):
            raise ValueError("'use_pbc' must be of type 'bool'")
    else:
        use_pbc = False

    from .. import measure

    default_distance_k = u("5 kcal mol-1 A-2")
    default_angle_k = u("50 kcal mol-1 rad-2")

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
                (
                    default_distance_k
                    if restraint_component == "distance"
                    else default_angle_k
                )
            ]
        elif type(restraint_components[restraint_component]["input_k"]) is not list:
            # Populate the list with the single specified value.
            restraint_components[restraint_component]["validated_k"] = n_measures * [
                u(restraint_components[restraint_component]["input_k"])
            ]
        else:
            if len(restraint_components[restraint_component]["input_k"]) == 0:
                # Empty list - populate with default values.
                restraint_components[restraint_component]["validated_k"] = (
                    n_measures
                    * [
                        (
                            default_distance_k
                            if restraint_component == "distance"
                            else default_angle_k
                        )
                    ]
                )
            elif len(restraint_components[restraint_component]["input_k"]) == 1:
                # List of length 1 - populate with that value.
                restraint_components[restraint_component]["validated_k"] = (
                    n_measures
                    * [u(restraint_components[restraint_component]["input_k"][0])]
                )
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
                restraint_components[restraint_component]["input_equil"] = (
                    n_measures
                    * [u(restraint_components[restraint_component]["input_equil"])]
                )
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
        b = BoreschRestraints(b)
    else:
        b = BoreschRestraints(name, b)

    # Set the use_pbc flag.
    b._use_pbc = use_pbc

    return b


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


def dihedral(mols, atoms, phi0=None, kphi=None, use_pbc=None, name=None, map=None):
    """
    Create a set of dihedral restraints from all of the atoms in 'atoms'
    where all atoms are contained in the container 'mols', using the
    passed values of the force constant 'kphi' and equilibrium
    torsion angle phi0.

    If phi0 is None, then the current torsional angle for
    the provided atoms will be used as the equilibium value.

    If kphi is None, then a default value of 100 kcal mol-1 rad-2 will be used

    Parameters
    ----------
    mols : sire.system._system.System
        The system containing the atoms.

    atoms : SireMol::Selector<SireMol::Atom>
        The atoms to restrain.

    kphi : str or SireUnits::Dimension::GeneralUnit or, optional
        The force constants for the torsion restraints.
        If None, this will default to 100 kcal mol-1 rad-2.
        Default is None.

    phi0 : str or SireUnits::Dimension::GeneralUnit, optional
        The equilibrium torsional angle for restraints. If None, these
        will be measured from the current coordinates of the atoms.
        Default is None.

    use_pbc : bool, optional
        Whether to use periodic boundary conditions when calculating
        the dihedral. Default is None.

    Returns
    -------
    DihedralRestraints : SireMM::DihedralRestraints
        A container of Dihedral restraints, where the first restraint is
        the DihedralRestraint created. The Dihedral restraint created can be
        extracted with DihedralRestraints[0].
    """
    from .. import u
    from ..base import create_map
    from ..mm import DihedralRestraint, DihedralRestraints

    map = create_map(map)
    map_dict = map.to_dict()
    kphi = kphi if kphi is not None else map_dict.get("kphi", None)
    phi0 = phi0 if phi0 is not None else map_dict.get("phi0", None)
    use_pbc = use_pbc if use_pbc is not None else map_dict.get("use_pbc", None)
    name = name if name is not None else map_dict.get("name", None)

    atoms = _to_atoms(mols, atoms)

    if len(atoms) != 4:
        raise ValueError(
            "You need to provide 4 atoms to create a dihedral restraint"
            f"whereas {len(atoms)} atoms were provided."
        )

    if use_pbc is not None:
        if not isinstance(use_pbc, bool):
            raise ValueError("'use_pbc' must be of type 'bool'")
    else:
        use_pbc = False

    from .. import measure

    if kphi is None:
        kphi = u("100 kcal mol-1 rad-2")

    elif type(kphi) is list:
        raise NotImplementedError(
            "Setup of multiple dihedral restraints simultaneously is not currently supported. "
            "Please set up each restraint individually and then combine them into multiple restraints."
        )

    if phi0 is None:
        from .. import measure

        phi0 = measure(atoms[0], atoms[1], atoms[2], atoms[3])

    elif type(phi0) is list:
        raise NotImplementedError(
            "Setup of multiple dihedral restraints simultaneously is not currently supported. "
            "Please set up each restraint individually and then combine them into multiple restraints."
        )
    else:
        phi0 = u(phi0)

    mols = mols.atoms()

    if name is None:
        restraints = DihedralRestraints()
    else:
        restraints = DihedralRestraints(name=name)

    restraints.add(DihedralRestraint(mols.find(atoms), phi0, kphi))

    # Set the use_pbc flag.
    restraints._use_pbc = use_pbc

    return restraints


def distance(mols, atoms0, atoms1, r0=None, k=None, use_pbc=None, name=None, map=None):
    """
    Create a set of distance restraints from all of the atoms in 'atoms0'
    to all of the atoms in 'atoms1' where all atoms are
    contained in the container 'mols', using the
    passed values of 'k' and equilibrium bond length r0.
    Note that 'k' corresponds to half the force constant, because
    the restraint energy is defined as k*(r - r0)**2 (hence the force is
    defined as 2*k*(r-r0)).

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

    if use_pbc is not None:
        if not isinstance(use_pbc, bool):
            raise ValueError("'use_pbc' must be of type 'bool'")
    else:
        use_pbc = True

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

    # Set the use_pbc flag.
    restraints._use_pbc = use_pbc

    return restraints


def morse_potential(
    mols,
    atoms0=None,
    atoms1=None,
    r0=None,
    k=None,
    de=None,
    use_pbc=None,
    name=None,
    auto_parametrise=False,
    map=None,
):
    """
    Create a set of Morse restraints from all of the atoms in 'atoms'
    where all atoms are contained in the container 'mols', using the
    passed values of the force constant, 'k', equilibrium
    distance value, r0, and well depth, de.

    If r0 is None, then the current distance for
    provided atoms will be used as the equilibium value.

    The potential energy of the Morse potential is defined as:

    e_morse=de*(1-exp(-sqrt(k/(2*de))*delta))^2
    where, delta=(r-r0). Additionally, if alchemical Morse potential is used,
    the potential is scaled as:

    rho*e_morse
    where rho is the lambda scaling parameter.

    Parameters
    ----------
    mols : sire.system._system.System
        The system containing the atoms.

    atoms0 : SireMol::Selector<SireMol::Atom>
        The first atom involved in the Morse restraint.

    atoms1 : SireMol::Selector<SireMol::Atom>
        The second atom involved in the Morse restraint.

    k : str or SireUnits::Dimension::GeneralUnit,
        The force constants for the Morse restraints. Optional if
        auto_parametrise is True.
        Default is None.

    r0 : str or SireUnits::Dimension::GeneralUnit, optional
        The equilibrium distance for the Morse restraints. If None, this
        will be measured from the current coordinates of the atoms.
        Default is None.

    de : str or SireUnits::Dimension::GeneralUnit
        The well depth (dissociation energy) for the Morse potential.
        Default is 100 kcal mol-1.

    use_pbc : bool, optional
        Whether to use periodic boundary conditions when calculating
        the distance. Default is None.

    name : str, optional
        The name of the restraint.
        Default is None.

    auto_parametrise : bool, optional
        If True, will attempt to automatically parametrise the Morse potential
        from a perturbation that annihilates a bond. This requires that 'mols'
        contains exactly one molecule that is perturbable, and that this
        molecule contains exactly one bond that is annihilated at lambda=1.
        The atoms involved in the annihilated bond will be used as 'atoms0'
        and 'atoms1', the equilibrium distance r0 will be set to the original
        bond length, and the force constant k will be set to the force constant
        of the bond in the unperturbed state. Note that 'de' must still be provided.
        Default is False.

    Returns
    -------
    MorsePotentialRestraints : SireMM::MorsePotentialRestraints
        A container of Morse restraints, where the first restraint is
        the MorsePotentialRestraint created. The Morse restraint created can be
        extracted with MorsePotentialRestraints[0].
    """

    from .. import u
    from ..base import create_map
    from ..mm import MorsePotentialRestraint, MorsePotentialRestraints
    from ..morph import link_to_reference

    map = create_map(map)

    if use_pbc is not None:
        if not isinstance(use_pbc, bool):
            raise ValueError("'use_pbc' must be of type 'bool'")
    else:
        use_pbc = False

    if auto_parametrise is False:
        if atoms0 is None or atoms1 is None:
            raise ValueError(
                "If auto_parametrise is False, then atoms0 and atoms1 must be provided"
            )
        if k is None:
            raise ValueError(
                "If auto_parametrise is False, then the force constant k must be provided"
            )
        elif isinstance(k, list):
            k = [u(x) for x in k]
        else:
            k = [u(k)]

    else:
        mol = mols.molecules("molecule property is_perturbable")
        ref_mol = link_to_reference(mol)

        if len(ref_mol) != 1:
            raise ValueError(
                "We need exactly one molecule that is perturbable to automatically "
                "set up the Morse potential restraints"
            )
        perturbable_mol = ref_mol[0]
        pert = perturbable_mol.perturbation(map=map)
        pert_omm = pert.to_openmm()
        changed_bonds = pert_omm.changed_bonds(to_pandas=False)

        # Attempt to find the bond that is annihilated at lambda=1
        for bond in changed_bonds:
            bond_name, length0, length1, k0, k1 = bond
            if k1 == 0:

                atom0_idx = [bond_name.atom0().index().value()][0]
                atom1_idx = [bond_name.atom1().index().value()][0]

                length0 = u(f"{length0} nm")

                # Divide k0 by 2 to convert from force constant to sire half
                # force constant k
                if k is None:
                    k0 = k0 / 2.0
                    k0 = u(f"{k0} kJ mol-1 nm-2")
                    k = [k0]

                # User can still override the force constant if they want, but
                # we need to ensure it's a list of units
                elif isinstance(k, list):
                    k = [u(x) for x in k]
                else:
                    k = [u(k)]

                # Translate the atom numbers to the original system indexes
                atoms0 = mols[
                    f"molecule property is_perturbable and atomidx {atom0_idx}"
                ]
                atoms1 = mols[
                    f"molecule property is_perturbable and atomidx {atom1_idx}"
                ]
                break

    try:
        atoms0 = _to_atoms(mols, atoms0)
        atoms1 = _to_atoms(mols, atoms1)
    except:
        raise ValueError("Unable to find atoms0 or atoms1 in the provided system")

    if atoms0.is_empty() or atoms1.is_empty():
        raise ValueError("We need at least one atom in each group")

    if len(atoms0) != len(atoms1):
        raise ValueError("atoms0 and atoms1 must be the same length")

    if len(atoms0) > 1 or len(atoms1) > 1:
        if not auto_parametrise:
            raise ValueError(
                "Setting up multiple Morse potential restraints at once is not currently supported."
                "Please set up each restraint individually and then combine them into multiple restraints."
            )

    if r0 is None:
        if auto_parametrise:
            r0 = [length0]
        else:
            # calculate all of the current distances
            from .. import measure

            r0 = []
            for atom0, atom1 in zip(atoms0, atoms1):
                r0.append(measure(atom0, atom1))
    elif type(r0) is list:
        r0 = [u(x) for x in r0]
    else:
        try:
            r0 = [u(r0)]
        except:
            raise ValueError(f"Unable to parse 'r0' as a Sire GeneralUnit: {r0}")

    if de is None:
        de = u("100 kcal mol-1")
    else:
        try:
            de = u(de)
        except:
            raise ValueError(f"Unable to parse 'de' as a Sire GeneralUnit: {de}")

    mols = mols.atoms()

    if name is None:
        restraints = MorsePotentialRestraints()
    else:
        restraints = MorsePotentialRestraints(name=name)

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
                "add the morse potential restraints."
            )

        if len(idxs1) == 0:
            raise KeyError(
                f"Could not find atom {atom1} in the molecules. Please ensure "
                "that 'mols' contains all of that atoms, or else we can't "
                "add the morse potential restraints."
            )

        if i < len(k):
            ik = k[i]
        else:
            ik = k[-1]

        if i < len(r0):
            ir0 = r0[i]
        else:
            ir0 = r0[-1]

        restraints.add(MorsePotentialRestraint(idxs0[0], idxs1[0], ik, ir0, de))

    # Set the use_pbc flag.
    restraints._use_pbc = use_pbc

    return restraints


def bond(*args, use_pbc=False, **kwargs):
    """
    Synonym for distance(), as a bond restraint is treated the same
    as a distance restraint
    """
    return distance(*args, use_pbc=use_pbc, **kwargs)


def inverse_distance(
    mols, atoms0, atoms1, r0=None, k=None, use_pbc=None, name=None, map=None
):
    """
    Create a set of inverse distance restraints from all of the atoms in 'atoms0'
    to all of the atoms in 'atoms1' where all atoms are
    contained in the container 'mols', using the
    passed values of 'k' and radius r0.
    Note that 'k' corresponds to half the force constant, because
    the restraint energy is defined as k*(r - r0)**2 (hence the force is
    defined as 2*k*(r-r0)).

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
    from ..mm import InverseBondRestraint, InverseBondRestraints

    map = create_map(map)

    if k is None:
        k = [u("150 kcal mol-1 A-2")]
    elif type(k) is list:
        k = [u(x) for x in k]
    else:
        k = [u(k)]

    if use_pbc is not None:
        if not isinstance(use_pbc, bool):
            raise ValueError("'use_pbc' must be of type 'bool'")
    else:
        use_pbc = True

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
        restraints = InverseBondRestraints()
    else:
        restraints = InverseBondRestraints(name=name)

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

        restraints.add(InverseBondRestraint(idxs0[0], idxs1[0], ik, ir0))

    # Set the use_pbc flag.
    restraints._use_pbc = use_pbc

    return restraints


def inverse_bond(*args, use_pbc=False, **kwargs):
    """
    Synonym for distance(), as a bond restraint is treated the same
    as a distance restraint
    """
    return inverse_distance(*args, use_pbc=use_pbc, **kwargs)


def positional(
    mols, atoms, k=None, r0=None, position=None, use_pbc=None, name=None, map=None
):
    """
    Create a set of position restraints for the atoms specified in
    'atoms' that are contained in the container 'mols', using the
    passed values of 'k' and flat-bottom potential
    well-width 'r0' for the restraints. Note that 'k' values
    correspond to half the force constants for the harmonic
    restraints, because the harmonic restraint energy is defined as
    k*(r - r0)**2 (hence the force is defined as 2*(r - r0)).

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

    if use_pbc is not None:
        if not isinstance(use_pbc, bool):
            raise ValueError("'use_pbc' must be of type 'bool'")
    else:
        use_pbc = True

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

    # Set the use_pbc flag.
    restraints._use_pbc = use_pbc

    return restraints


def rmsd(mols, atoms, ref=None, k=None, r0=None, name=None, map=None):
    """
    Create a set of RMSD restraints for the atoms specified in
    'atoms' that are contained in the container 'mols', using the
    passed values of 'k' and flat-bottom potential
    well-width 'r0' for the restraints. Note that 'k' values
    correspond to half the force constants for the harmonic
    restraints, because the harmonic restraint energy is defined as
    k*(rmsd - r0)**2 (hence the force is defined as 2*(rmsd - r0)).

    The RMSD calculation is perfomed by default using the position
    of mols. Optionally, a different state of the system can be
    supplied as a reference by passing the 'ref' argument.

    If 'r0' is not specified, then a simple harmonic restraint
    is used.

    If 'k' is not specified, then a default of 150 kcal mol-1 A-2
    will be used.

    Parameters
    ----------
    mols : sire.system._system.System
        The system containing the atoms.

    atoms : SireMol::Selector<SireMol::Atom>
        The atoms to restrain.

    ref : sire.system._system.System
        The system from which the reference positions for the RMSD calculation
        are extracted from. If None, this will default to the current
        state of mols.

    k : str or SireUnits::Dimension::GeneralUnit or, optional
        The force constant for the RMSD restraints.
        If None, this will default to 150 kcal mol-1 A-2.
        Default is None.

    r0 : str or SireUnits::Dimension::GeneralUnit, optional
        The width of the flat bottom restraint. If None, this is zero
        and a simple harmonic restraint is used.
        Default is None.

    Returns
    -------
    RMSDRestraints : SireMM::RMSDRestraints
        A container of RMSD restraints, where the first restraint is
        the RMSDRestraint created. The RMSD restraint created can be
        extracted with RMSDRestraints[0].
    """
    from .. import u
    from ..base import create_map
    from ..mm import RMSDRestraint, RMSDRestraints

    map = create_map(map)

    if k is None:
        k = u("150 kcal mol-1 A-2")
    else:
        k = u(k)

    if r0 is None:
        r0 = u("0")
    else:
        r0 = u(r0)

    atoms = _to_atoms(mols, atoms)
    mols = mols.atoms()

    if name is None:
        restraints = RMSDRestraints()
    else:
        restraints = RMSDRestraints(name=name)

    # Set default reference positions to mols
    if ref is None:
        ref = mols
    else:
        try:
            ref = ref.atoms()
        except AttributeError:
            raise TypeError("The reference state must be a complete system.")

    # Generate list of all positions as reference for RMSD calculation
    ref_pos = ref.atoms().property("coordinates")

    restraints.add(RMSDRestraint(mols.find(atoms), ref_pos, k, r0))

    return restraints
