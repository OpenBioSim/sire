__all__ = ["repartition_hydrogen_masses"]


def repartition_hydrogen_masses(
    mols,
    mass_factor=1.5,
    ignore_water: bool = False,
    ignore_non_water: bool = False,
    include_end_states: bool = True,
    map=None,
):
    """
    Increase the mass of hydrogen atoms to hmass times * amu, and subtract the
    mass increase from the heavy atom the hydrogen is bonded to.

    (based heavily on the repartitionMasses function in
     Sire.Tools.OpenMMMD)

    Parameters
    ----------

    mol : sire.mol.Molecule or list of molecules, or System
        The molecule(s) whose hydrogen masses should be repartitioned

    mass_factor : float
        The factor to multiply the mass of hydrogen atoms by. Using
        the default of 1.5 is known to be a good value to use to
        achieve a 4 fs timestep with the (default) LangevinMiddle integrator

    ignore_water : bool
        Whether or not to ignore water molecules (default False)

    ignore_non_water : bool
        Whether or not to ignore non-water molecules (default False)

    include_end_states : bool
        Whether or not to repartition the masses of the end states
        of perturbable molecules (default True) (i.e. this will
        automatically repartition `mass0` and `mass1` in addition
        to `mass`)

    map : dict
        The property map used to identify molecular properties

    Returns
    -------

    sire.mol.Molecule, list of molecules or System
        The repartitioned molecule(s)
    """

    try:
        from ..legacy.IO import (
            repartitionHydrogenMass as _repartition_hydrogen_mass,
        )
    except Exception:
        from ..legacy.IO import (
            repartition_hydrogen_mass as _repartition_hydrogen_mass,
        )

    if ignore_water and ignore_non_water:
        # ignoring everything ;-)
        return mols

    from ..base import create_map

    map = create_map(map)

    # convert the flag into the integer defined in the BioSimSpace function
    if ignore_water:
        water = 0
    elif ignore_non_water:
        water = 2
    else:
        water = 1

    mols = mols.clone()

    for mol in mols.molecules():
        mol = _repartition_hydrogen_mass(mol, mass_factor, water, map)

        if include_end_states:
            mass0 = f"{map['mass'].source()}0"
            mass1 = f"{map['mass'].source()}1"

            if mol.has_property(mass0):
                map0 = map.clone()
                map0.set("mass", mass0)
                mol = _repartition_hydrogen_mass(mol, mass_factor, water, map0)

            if mol.has_property(mass1):
                map1 = map.clone()
                map1.set("mass", mass1)
                mol = _repartition_hydrogen_mass(mol, mass_factor, water, map1)

        mols.update(mol)

    return mols
