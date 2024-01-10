def _configure_engine(engine, mols, qm_mol, map):
    """
    Internal helper function to configure a QM engine ready for dynamics.
    """

    # Work out the indices of the QM atoms.
    try:
        atoms_to_find = qm_mol.atoms()
        idxs = mols.atoms().find(atoms_to_find)
        engine.set_atoms(idxs)
    except:
        raise Exception("Unable to set atom indices for the QM region.")

    # Work out the atomic numbers of the QM atoms.
    try:
        elem_prop = map["element"]
        numbers = [
            atom.property(f"{elem_prop}").num_protons() for atom in atoms_to_find
        ]
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

    return engine


def _create_merged_mol(qm_mol, map):
    """
    Internal helper function to create a merged molecule from the QM molecule.
    """

    from ..legacy import CAS as _CAS
    from ..legacy import Mol as _Mol
    from ..legacy import MM as _MM

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
        # Bonds.
        if prop == bond_prop:
            # Copy the bonds to the lambda = 1 state.
            edit_mol = edit_mol.set_property(
                prop + "1", qm_mol.property(prop)
            ).molecule()

            # Create an equivalent set of bonds for the lambda = 0 state with
            # zeroed force constants.

            bonds = _MM.TwoAtomFunctions(info)

            for bond in qm_mol.property(prop).potentials():
                # Extract the atoms involved in the bond.:w
                atom0 = info.atom_idx(bond.atom0())
                atom1 = info.atom_idx(bond.atom1())

                # CAS variable for the bond.
                r = _CAS.Symbol("r")

                # Create a null AmberBond.
                amber_bond = _MM.AmberBond(0, r)

                # Set the new bond.
                bonds.set(atom0, atom1, amber_bond.to_expression(r))

            # Set the bonds for the lambda = 0 state.
            edit_mol = edit_mol.set_property(prop + "0", bonds).molecule()

            # Finally, delete the existing property.
            edit_mol = edit_mol.remove_property(prop).molecule()

        # Angles.
        elif prop == angle_prop:
            # Copy the angles to the lambda = 1 state.
            edit_mol = edit_mol.set_property(
                prop + "1", qm_mol.property(prop)
            ).molecule()

            # Create an equivalent set of angles for the lambda = 0 state with
            # zeroed force constants.

            angles = _MM.ThreeAtomFunctions(info)

            for angle in qm_mol.property(prop).potentials():
                # Extract the atoms involved in the angle.
                atom0 = info.atom_idx(angle.atom0())
                atom1 = info.atom_idx(angle.atom1())
                atom2 = info.atom_idx(angle.atom2())

                # CAS variable for the angle.
                theta = _CAS.Symbol("theta")

                # Create a null AmberAngle.
                amber_angle = _MM.AmberAngle(0.0, theta)

                # Set the new angle.
                angles.set(atom0, atom1, atom2, amber_angle.to_expression(theta))

            # Set the angles for the lambda = 0 state.
            edit_mol = edit_mol.set_property(prop + "0", angles).molecule()

            # Finally, delete the existing property.
            edit_mol = edit_mol.remove_property(prop).molecule()

        # Dihedrals.
        elif prop == dihedral_prop:
            # Copy the dihedrals to the lambda = 1 state.
            edit_mol = edit_mol.set_property(
                prop + "1", qm_mol.property(prop)
            ).molecule()

            # Create an equivalent set of dihedrals for the lambda = 0 state
            # with zeroed force constants.

            dihedrals = _MM.FourAtomFunctions(info)

            for dihedral in qm_mol.property(prop).potentials():
                # Extract the atoms involved in the dihedral.
                atom0 = info.atom_idx(dihedral.atom0())
                atom1 = info.atom_idx(dihedral.atom1())
                atom2 = info.atom_idx(dihedral.atom2())
                atom3 = info.atom_idx(dihedral.atom3())

                # CAS varialbe for the dihedral.
                phi = _CAS.Symbol("phi")

                # Create a hull AmberDihedral.
                amber_dihedral = _MM.AmberDihedral(0, phi)

                # Set the new dihedral.
                dihedrals.set(
                    atom0,
                    atom1,
                    atom2,
                    atom3,
                    amber_dihedral.to_expression(phi),
                )

            # Set the dihedrals for the lambda = 0 state.
            edit_mol = edit_mol.set_property(prop + "0", dihedrals).molecule()

            # Finally, delete the existing property.
            edit_mol = edit_mol.remove_property(prop).molecule()

        # Impropers.
        elif prop == improper_prop:
            # Copy the impropers to the lambda = 1 state.
            edit_mol = edit_mol.set_property(
                prop + "1", qm_mol.property(prop)
            ).molecule()

            # Create an equivalent set of impropers for the lambda = 0 state
            # with zeroed force constants.

            impropers = _MM.FourAtomFunctions(info)

            for improper in qm_mol.property(prop).potentials():
                # Extract the atoms involved in the improper.
                atom0 = info.atom_idx(improper.atom0())
                atom1 = info.atom_idx(improper.atom1())
                atom2 = info.atom_idx(improper.atom2())
                atom3 = info.atom_idx(improper.atom3())

                # CAS variable for the improper.
                psi = _CAS.Symbol("psi")

                # Create a null AmberDihedral.
                amber_improper = _MM.AmberDihedral(0, psi)

                # Set the new improper.
                impropers.set(
                    atom0,
                    atom1,
                    atom2,
                    atom3,
                    amber_improper.to_expression(psi),
                )

            # Set the impropers for the lambda = 0 state.
            edit_mol = edit_mol.set_property(prop + "0", impropers).molecule()

            # Finally, delete the existing property.
            edit_mol = edit_mol.remove_property(prop).molecule()

        # Charge.
        elif prop == charge_prop:
            # Copy the charges to the lambda = 1 state.
            edit_mol = edit_mol.set_property(
                prop + "1", qm_mol.property(prop)
            ).molecule()

            # Create a set of null charges for the lambda = 0 state.
            charges = _Mol.AtomCharges(info)
            edit_mol = edit_mol.set_property(prop + "0", charges).molecule()

            # Finally, delete the existing property.
            edit_mol = edit_mol.remove_property(prop).molecule()

        # Connectivity.
        elif prop == connectivity_prop:
            pass

        # Everything else.
        else:
            # All other properties remain the same in both end states.
            edit_mol = edit_mol.set_property(
                prop + "0", qm_mol.property(prop)
            ).molecule()
            edit_mol = edit_mol.set_property(
                prop + "1", qm_mol.property(prop)
            ).molecule()

            # Delete the existing property.
            edit_mol = edit_mol.remove_property(prop).molecule()

    # Flag the molecule as perturbable.
    edit_mol = edit_mol.set_property(map["is_perturbable"], True).molecule()

    # Commit the changes.
    qm_mol = edit_mol.commit()

    # Link to the perturbation to the reference state.
    qm_mol = qm_mol.perturbation().link_to_reference().commit()

    # Return the merged molecule.
    return qm_mol
