def _configure_engine(engine, mols, qm_atoms, link_atoms, map):
    """
    Internal helper function to configure a QM engine ready for dynamics.
    """

    # Work out the indices of the QM atoms.
    try:
        idxs = mols.atoms().find(qm_atoms)
        engine.set_atoms(idxs)
    except:
        raise Exception("Unable to set atom indices for the QM region.")

    # Work out the atomic numbers of the QM atoms.
    try:
        elem_prop = map["element"]
        numbers = [atom.property(f"{elem_prop}").num_protons() for atom in qm_atoms]
        engine.set_numbers(numbers)
    except:
        raise Exception("Unable to set atomic numbers for the QM region.")

    # Work out the atomic charge for all atoms in the system.
    try:
        charge_prop = map["charge"]
        charges = [atom.property(f"{charge_prop}").value() for atom in mols.atoms()]
        engine.set_charges(charges)
    except:
        raise Exception("Unable to set atomic charges for the system.")

    # Set the link atom information.
    try:
        engine.set_link_atoms(link_atoms)
    except:
        raise Exception("Unable to set link atom information.")

    return engine


def _create_merged_mol(qm_mol, qm_atoms, map):
    """
    Internal helper function to create a merged molecule from the QM molecule.
    """

    from ..legacy import CAS as _CAS
    from ..legacy import Mol as _Mol
    from ..legacy import MM as _MM

    # Store the indices of the QM atoms.
    qm_idxs = [atom.index() for atom in qm_atoms]

    # Get the user defined names for the properties that we need to
    # merge.
    bond_prop = map["bond"]
    angle_prop = map["angle"]
    dihedral_prop = map["dihedral"]
    improper_prop = map["improper"]
    charge_prop = map["charge"]
    connectivity_prop = map["connectivity"]

    # Get the molecular info object.
    info = qm_mol.info()

    # Make an editable version of the molecule.
    edit_mol = qm_mol.edit()

    for prop in qm_mol.property_keys():
        # For all bonded properties we copy the MM terms to the lambda = 0
        # (MM) end state, create a null set of terms for the lambda = 1 (QM)
        # end state for any terms that only involve QM atoms,  then remove
        # the existing property. Charges also zeroed for the MM end state.
        # All other properties remain the same in both states.

        # Bonds.
        if prop == bond_prop:
            edit_mol = edit_mol.set_property(
                prop + "0", qm_mol.property(prop)
            ).molecule()

            bonds = _MM.TwoAtomFunctions(info)

            for bond in qm_mol.property(prop).potentials():
                atom0 = info.atom_idx(bond.atom0())
                atom1 = info.atom_idx(bond.atom1())

                if atom0 in qm_idxs and atom1 in qm_idxs:
                    r = _CAS.Symbol("r")
                    amber_bond = _MM.AmberBond(0, r)
                    bonds.set(atom0, atom1, amber_bond.to_expression(r))
                else:
                    bonds.set(atom0, atom1, bond.function())

            edit_mol = edit_mol.set_property(prop + "1", bonds).molecule()
            edit_mol = edit_mol.remove_property(prop).molecule()

        # Angles.
        elif prop == angle_prop:
            edit_mol = edit_mol.set_property(
                prop + "0", qm_mol.property(prop)
            ).molecule()

            angles = _MM.ThreeAtomFunctions(info)

            for angle in qm_mol.property(prop).potentials():
                atom0 = info.atom_idx(angle.atom0())
                atom1 = info.atom_idx(angle.atom1())
                atom2 = info.atom_idx(angle.atom2())

                if atom0 in qm_idxs and atom1 in qm_idxs and atom2 in qm_idxs:
                    theta = _CAS.Symbol("theta")
                    amber_angle = _MM.AmberAngle(0.0, theta)
                    angles.set(atom0, atom1, atom2, amber_angle.to_expression(theta))
                else:
                    angles.set(atom0, atom1, atom2, angle.function())

            edit_mol = edit_mol.set_property(prop + "1", angles).molecule()
            edit_mol = edit_mol.remove_property(prop).molecule()

        # Dihedrals.
        elif prop == dihedral_prop:
            edit_mol = edit_mol.set_property(
                prop + "0", qm_mol.property(prop)
            ).molecule()

            dihedrals = _MM.FourAtomFunctions(info)

            for dihedral in qm_mol.property(prop).potentials():
                atom0 = info.atom_idx(dihedral.atom0())
                atom1 = info.atom_idx(dihedral.atom1())
                atom2 = info.atom_idx(dihedral.atom2())
                atom3 = info.atom_idx(dihedral.atom3())

                if (
                    atom0 in qm_idxs
                    and atom1 in qm_idxs
                    and atom2 in qm_idxs
                    and atom3 in qm_idxs
                ):
                    phi = _CAS.Symbol("phi")
                    amber_dihedral = _MM.AmberDihedral(0, phi)
                    dihedrals.set(
                        atom0,
                        atom1,
                        atom2,
                        atom3,
                        amber_dihedral.to_expression(phi),
                    )
                else:
                    dihedrals.set(atom0, atom1, atom2, atom3, dihedral.function())

            edit_mol = edit_mol.set_property(prop + "1", dihedrals).molecule()
            edit_mol = edit_mol.remove_property(prop).molecule()

        # Impropers.
        elif prop == improper_prop:
            edit_mol = edit_mol.set_property(
                prop + "0", qm_mol.property(prop)
            ).molecule()

            impropers = _MM.FourAtomFunctions(info)

            for improper in qm_mol.property(prop).potentials():
                atom0 = info.atom_idx(improper.atom0())
                atom1 = info.atom_idx(improper.atom1())
                atom2 = info.atom_idx(improper.atom2())
                atom3 = info.atom_idx(improper.atom3())

                if (
                    atom0 in qm_idxs
                    and atom1 in qm_idxs
                    and atom2 in qm_idxs
                    and atom3 in qm_idxs
                ):
                    psi = _CAS.Symbol("psi")
                    amber_improper = _MM.AmberDihedral(0, psi)
                    impropers.set(
                        atom0,
                        atom1,
                        atom2,
                        atom3,
                        amber_improper.to_expression(psi),
                    )
                else:
                    impropers.set(atom0, atom1, atom2, atom3, improper.function())

            edit_mol = edit_mol.set_property(prop + "1", impropers).molecule()
            edit_mol = edit_mol.remove_property(prop).molecule()

        # Charge.
        elif prop == charge_prop:
            edit_mol = edit_mol.set_property(
                prop + "0", qm_mol.property(prop)
            ).molecule()

            charges = _Mol.AtomCharges(info)
            edit_mol = edit_mol.set_property(prop + "1", charges).molecule()
            edit_mol = edit_mol.remove_property(prop).molecule()

        # Connectivity.
        elif prop == connectivity_prop:
            pass

        # All other properties remain the same in both end states.
        else:
            edit_mol = edit_mol.set_property(
                prop + "0", qm_mol.property(prop)
            ).molecule()
            edit_mol = edit_mol.set_property(
                prop + "1", qm_mol.property(prop)
            ).molecule()

            edit_mol = edit_mol.remove_property(prop).molecule()

    # Flag the molecule as perturbable.
    edit_mol = edit_mol.set_property(map["is_perturbable"], True).molecule()

    # Commit the changes.
    qm_mol = edit_mol.commit()

    # Link to the perturbation to the reference state.
    qm_mol = qm_mol.perturbation().link_to_reference().commit()

    # Return the merged molecule.
    return qm_mol


def _get_link_atoms(mols, qm_mol, qm_atoms, map):
    """
    Internal helper function to get a dictionary with link atoms for each QM atom.
    """

    from ..legacy.Mol import Connectivity as _Connectivity
    from ..legacy.Mol import CovalentBondHunter as _CovalentBondHunter

    # Store the indices of the QM atoms.
    qm_idxs = [atom.index() for atom in qm_atoms]

    # Create a connectivity object.
    connectivity = _Connectivity(qm_mol, _CovalentBondHunter(), map)

    mm1_atoms = {}

    # Loop over the QM atoms and find any MM atoms that are bonded to them.
    for atom in qm_atoms:
        # Store the atom index.
        idx = atom.index()

        # Get the bonds for the atom.
        bonds = connectivity.get_bonds(idx)

        # A list to hold MM atoms involved in the bonds.
        mm_bonds = []

        # A flag to indicate if the atom has a bond to another QM atom.
        has_qm_bond = False

        # Loop over the bonds and find the MM atoms.
        for bond in bonds:
            # Get the indices of the two atoms in the bond.
            idx0 = bond.atom0()
            idx1 = bond.atom1()

            # Work out which atom isn't the current QM atom.
            if idx0 != idx:
                bond_idx = idx0
            else:
                bond_idx = idx1

            # If the atom is not in the QM region, add it to the list.
            if bond_idx not in qm_idxs:
                mm_bonds.append(bond_idx)
            else:
                has_qm_bond = True

        # If there are no QM bonds for this atom, raise an exception.
        if not has_qm_bond:
            raise ValueError(
                f"Atom {idx} in the QM region has no bonds to other QM atoms!"
            )

        # Store the list of MM atoms.
        if len(mm_bonds) > 0:
            mm1_atoms[idx] = mm_bonds

    # Now work out the MM atoms that are bonded to the link atoms.
    mm2_atoms = {}
    for k, v in mm1_atoms.items():
        for idx in v:
            if idx not in mm2_atoms:
                bonds = connectivity.get_bonds(idx)
                mm_bonds = []
                for bond in bonds:
                    idx0 = bond.atom0()
                    idx1 = bond.atom1()
                    if idx0 != idx:
                        bond_idx = idx0
                    else:
                        bond_idx = idx1
                    if bond_idx not in qm_idxs:
                        mm_bonds.append(bond_idx)
                mm2_atoms[idx] = mm_bonds

    # Convert MM1 atoms dictionary to absolute indices.
    abs_mm1_atoms = {}
    for k, v in mm1_atoms.items():
        qm_idx = mols.atoms().find(qm_mol.atoms()[k])
        link_idx = [mols.atoms().find(qm_mol.atoms()[x]) for x in v]
        abs_mm1_atoms[qm_idx] = link_idx

    # Convert MM2 atoms dictionary to absolute indices.
    abs_mm2_atoms = {}
    for k, v in mm2_atoms.items():
        link_idx = mols.atoms().find(qm_mol.atoms()[k])
        mm_idx = [mols.atoms().find(qm_mol.atoms()[x]) for x in v]
        abs_mm2_atoms[link_idx] = mm_idx

    return abs_mm1_atoms, abs_mm2_atoms
