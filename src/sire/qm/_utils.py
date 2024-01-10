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
        # For all bonded properties we copy the MM terms to the lambda = 0
        # (MM) end state, create a null set of terms for the lambda = 1 (QM)
        # end state, then remove the existing property. Charges also zeroed
        # for the MM end state. All other properties remain the same in both
        # states.

        # Bonds.
        if prop == bond_prop:
            edit_mol = edit_mol.set_property(
                prop + "0", qm_mol.property(prop)
            ).molecule()

            bonds = _MM.TwoAtomFunctions(info)

            for bond in qm_mol.property(prop).potentials():
                atom0 = info.atom_idx(bond.atom0())
                atom1 = info.atom_idx(bond.atom1())

                r = _CAS.Symbol("r")
                amber_bond = _MM.AmberBond(0, r)

                bonds.set(atom0, atom1, amber_bond.to_expression(r))

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

                theta = _CAS.Symbol("theta")
                amber_angle = _MM.AmberAngle(0.0, theta)

                angles.set(atom0, atom1, atom2, amber_angle.to_expression(theta))

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

                phi = _CAS.Symbol("phi")
                amber_dihedral = _MM.AmberDihedral(0, phi)

                dihedrals.set(
                    atom0,
                    atom1,
                    atom2,
                    atom3,
                    amber_dihedral.to_expression(phi),
                )

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

                psi = _CAS.Symbol("psi")
                amber_improper = _MM.AmberDihedral(0, psi)

                impropers.set(
                    atom0,
                    atom1,
                    atom2,
                    atom3,
                    amber_improper.to_expression(psi),
                )

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
