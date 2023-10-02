__all__ = ["repartition_hydrogen_masses"]


def repartition_hydrogen_masses(mol, mass_factor=4.0, map=None):
    """
    Increase the mass of hydrogen atoms to hmass times * amu, and subtract the
    mass increase from the heavy atom the hydrogen is bonded to.

    (based heavily on the repartitionMasses function in
     Sire.Tools.OpenMMMD)

    Parameters
    ----------

    mol : sire.mol.Molecule
        The molecule whose hydrogen masses should be repartitioned

    mass_factor : float
        The factor to multiply the mass of hydrogen atoms by. Using
        the default of 4 will ensure that the atoms. Note that this
        value can only be set to between 1 and 4.

    map : dict
        The property map used to identify molecular properties

    Returns
    -------

    sire.mol.Molecule
        The repartitioned molecule
    """

    # make sure that this is the whole molecule
    mol = mol.molecule()

    if mass_factor < 1 or mass_factor > 4:
        raise ValueError(
            f"The mass factor must be between 1 and 4, but was {mass_factor}"
        )

    from ..base import create_map
    from ..units import g_per_mol

    map = create_map(map)

    atoms = mol.atoms()

    if len(atoms) <= 1:
        # nothing to do
        return mol

    connectivity = mol.property(map["connectivity"])
    atom_masses = {}

    #
    # First pass. Initialise changes in atom_masses to effect
    #
    for atom in atoms:
        atom_masses[atom.index().value()] = 0 * g_per_mol

    total_delta = 0 * g_per_mol

    #
    # Second pass. Decide how the mass of each atom should change.
    #
    mass_property = map["mass"].source()

    for atom in atoms:
        atom_index = atom.index()
        atom_mass = atom.property(mass_property)

        # units are in g_per_mol
        if atom_mass.value() < 1.1:
            # Atoms with a mass < 1.1 g_per_mol are assumed to be hydrogen
            # atoms
            delta_mass = (atom_mass * mass_factor) - atom_mass
            # print("Increasing mass %s by %s  " % (at, deltamass))
            total_delta += delta_mass
            atom_masses[atom_index.value()] = delta_mass

            # Skip monoatomic systems without connectivity property
            if connectivity is None:
                continue

            bonds = connectivity.get_bonds(atom_index)

            # Get list of atoms that share one bond with this atom. Ignore
            # all atoms that have a mass < 1.1 g_mol in the ORIGINAL atoms
            # list.For each atom this was bonded to, substract
            # delta_mass / nbonded
            bonded_atoms = []
            for bond in bonds:
                atom0 = mol.atom(bond.atom0()).index()
                atom1 = mol.atom(bond.atom1()).index()
                if atom0 == atom_index:
                    heavy_atom_index = atom1
                else:
                    heavy_atom_index = atom0

                if heavy_atom_index in bonded_atoms:
                    continue

                heavy_atom = mol.atom(heavy_atom_index)
                heavy_atom_mass = heavy_atom.property(mass_property)

                # units are in g_per_mol
                if heavy_atom_mass.value() < 1.1:
                    continue

                bonded_atoms.append(heavy_atom_index)

            for bonded_atom in bonded_atoms:
                total_delta -= delta_mass
                atom_masses[bonded_atom.value()] -= delta_mass

    # Sanity check (g_per_mol)
    if total_delta.value() > 0.001:
        from ..utils import Console

        Console.warning(
            "WARNING! The mass repartitioning algorithm is not conserving "
            f"atomic masses for molecule {mol} (total delta is {total_delta})."
        )

    # Now that have worked out mass changes per atom, update the molecule
    c = mol.cursor()

    for atom in c.atoms():
        atom_index = atom.id().value()
        atom_mass = atom[mass_property]
        new_mass = atom_mass + atom_masses[atom_index]

        # Sanity check. Note this is likely to occur if mass_factor > 4
        if new_mass.value() < 0.0:
            raise ValueError(
                f"WARNING! The mass of {atom} is less than zero after "
                "hydrogen mass repartitioning. This should not happen! "
                "Decrease the hydrogen mass repartitioning factor "
                "and try again."
            )

        atom[mass_property] = new_mass

    return c.commit()
