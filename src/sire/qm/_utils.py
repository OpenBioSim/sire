def _check_charge(qm_atoms, map, tol=1e-6):
    """
    Internal helper function to check that the QM region has integer charge.

    Parameters
    ----------

    qm_atoms: [sire.legacy.Mol.AtomIdx]
        A list of QM atoms.

    map: sire.legacy.Base.PropertyMap
        The property map for the molecule.

    tol: float
        The tolerance for the charge check.

    Raises
    ------

    Exception
        If the charge of the QM region is not an integer.
    """

    import math as _math

    # Get the charge property.
    charge_prop = map["charge"]

    # Work out the charge of the QM atoms.
    qm_charge = 0
    for atom in qm_atoms:
        qm_charge += atom.property(charge_prop).value()

    # Check that the charge is an integer.
    if not _math.isclose(qm_charge, round(qm_charge), abs_tol=tol):
        raise Exception(f"Charge of the QM region ({qm_charge}) is not an integer!")


def _create_qm_mol_to_atoms(qm_atoms):
    """
    Internal helper function to create a mapping between molecule numbers and
    a list of QM atoms.

    Parameters
    ----------

    qm_atoms: [sire.legacy.Mol.AtomIdx]
        A list of QM atoms.

    Returns
    -------

    qm_mol_to_atoms: {int: [sire.legacy.Mol.AtomIdx]}
        A dictionary with molecule numbers as keys and a list of QM atoms as
        values.
    """
    qm_mol_to_atoms = {}
    for atom in qm_atoms:
        mol_num = atom.molecule().number()
        if mol_num not in qm_mol_to_atoms:
            qm_mol_to_atoms[mol_num] = [atom]
        else:
            qm_mol_to_atoms[mol_num].append(atom)

    return qm_mol_to_atoms


def _check_qm_atom_bonds(mol, atom, qm_idxs, map):
    """
    Internal helper function to check the bonding for QM atoms.

    Parameters
    ----------

    mol: sire.legacy.Mol.Molecule
        The molecule containing the QM atoms.

    atom: sire.legacy.Mol.Atom
        The QM atom to check the bonding for.

    qm_idxs: [sire.legacy.Mol.AtomIdx]
        The indices of the QM atoms.

    map: sire.legacy.Base.PropertyMap
        The property map for the molecule.

    Returns
    -------

    mm_atoms: [sire.legacy.Mol.AtomIdx]
        A list of MM atoms that are bonded to the QM atom.

    has_qm_bond: bool
        A flag to indicate if the QM atom has a bond to another QM atom.
    """

    # Get the bonds for the molecule.
    bonds = mol.property(map["bond"]).potentials()

    # Store the info for the molecule.
    info = mol.info()

    # Store the cut-group atom index pair.
    cg_atom_idx = info.cg_atom_idx(atom.index())

    # Initialise a list to store the MM atoms.
    mm_atoms = []

    # A flag to indicate if the atom has a bond to another QM atom.
    has_qm_bond = False

    # Loop over all of the bonds.
    for bond in bonds:
        # Get the indices of the atoms in the bond.
        idx0 = bond.atom0()
        idx1 = bond.atom1()

        # Work out which is the other atom in the bond.
        if idx0 == cg_atom_idx:
            idx = idx1
        elif idx1 == cg_atom_idx:
            idx = idx0
        else:
            continue

        # Convert to an atom index.
        idx = info.atom_idx(idx)

        # The atom is not in the QM region.
        if idx not in qm_idxs:
            mm_atoms.append(idx)
        else:
            has_qm_bond = True

    return mm_atoms, has_qm_bond


def _get_link_atoms(mols, qm_mol_to_atoms, map):
    """
    Internal helper function to get a dictionary with link atoms for each QM atom.

    Parameters
    ----------

    mols: sire.legacy.System.System
        The Sire system containing the QM atoms.

    qm_mol_to_atoms: {sire.legacy.Mol.MolNum: [sire.legacy.Mol.AtomIdx]}
        A dictionary with molecule numbers as keys and a list of QM atoms as
        values.

    map: sire.legacy.Base.PropertyMap
        The property map for the system.

    Returns

    mm1_to_qm: {sire.legacy.Mol.AtomIdx: sire.legacy.Mol.AtomIdx}
        A dictionary with link atoms as keys and QM atoms as values.

    mm1_to_mm2: {sire.legacy.Mol.AtomIdx: [sire.legacy.Mol.AtomIdx]}
        A dictionary with link atoms as keys and a list of MM atoms as values.

    bond_scale_factors: {sire.legacy.Mol.AtomIdx: float}
        A dictionary with link atoms as keys and bond scale factors as values.

    mm1_indices: [[sire.legacy.Mol.AtomIdx]]
        A list of lists of MM1 atom indices.
    """

    import warnings as _warnings

    from ..legacy import CAS as _CAS
    from ..legacy import Mol as _Mol
    from ..legacy import MM as _MM

    # Initialise the dictionaries.

    # List of MM1 atom indices as sire.legacy.Mol.AtomIdx objects.
    mm1_indices = []

    # Link atoms to QM atoms.
    mm1_to_qm = {}

    # Link atoms to MM atoms.
    mm1_to_mm2 = {}

    # QM to link atom bond scale factors. These are the ratios of the equilibrium
    # bond lengths for the QM-L bonds and the QM-MM1 bonds, taken from the MM
    # bond potentials, i.e. R0(QM-L) / R0(QM-MM1).
    bond_scale_factors = {}

    # Store carbon and hydrogen elements.
    carbon = _Mol.Element("C")
    hydrogen = _Mol.Element("H")

    # Store the element property.
    elem_prop = map["element"]

    # Loop over all molecules containing QM atoms.
    for mol_num, qm_atoms in qm_mol_to_atoms.items():
        # Get the molecule.
        qm_mol = qm_atoms[0].molecule()

        # Store the indices of the QM atoms.
        qm_idxs = [atom.index() for atom in qm_atoms]

        # Create a connectivity object.
        connectivity = _Mol.Connectivity(qm_mol, _Mol.CovalentBondHunter(), map)

        # Dictionary to store the MM1 atoms.
        mm1_atoms = {}

        # Loop over the QM atoms and find any MM atoms that are bonded to them.
        for atom in qm_atoms:
            # Store the atom index.
            idx = atom.index()

            # Store the element of the atom.
            elem = atom.property(elem_prop)

            # Check the bonding for this atom.
            mm_atoms, has_qm_bond = _check_qm_atom_bonds(qm_mol, atom, qm_idxs, map)

            # If there are no QM bonds for this atom, raise an exception.
            if not has_qm_bond:
                raise ValueError(
                    f"Atom {idx} in the QM region has no bonds to other QM atoms!"
                )

            # Store the list of MM atoms.
            if len(mm_atoms) > 0:
                if len(mm_atoms) > 1:
                    raise Exception(f"QM atom {idx} has more than one MM bond!")
                else:
                    # Get the element of the cut atom.
                    link_elem = qm_mol[mm_atoms[0]].property(elem_prop)

                    # If the element is hydrogen, raise an exception.
                    if elem == hydrogen:
                        abs_idx = mols.atoms().find(atom)
                        raise Exception(
                            "Attempting replace a hydrogen with a link atom "
                            f"(atom index {abs_idx})!"
                        )

                    # Warn if the link atom is not for a carbon-carbon bond.
                    if elem != carbon or link_elem != carbon:
                        abs_idx = mols.atoms().find(atom)
                        _warnings.warn(
                            "Attempting to add a link atom for a non carbon-carbon "
                            f"bond (atom index {abs_idx})!"
                        )

                    # Store the link (MM1) atom.
                    mm1_atoms[idx] = mm_atoms[0]

        # Now work out the MM atoms that are bonded to the link atoms. (MM2 atoms.)
        mm2_atoms = {}
        for qm_idx, mm1_idx in mm1_atoms.items():
            if mm1_idx not in mm2_atoms:
                bonds = connectivity.get_bonds(mm1_idx)
                mm_bonds = []
                for bond in bonds:
                    idx0 = bond.atom0()
                    idx1 = bond.atom1()
                    if idx0 != mm1_idx:
                        bond_idx = idx0
                    else:
                        bond_idx = idx1
                    if bond_idx not in qm_idxs:
                        mm_bonds.append(bond_idx)
                mm2_atoms[mm1_idx] = mm_bonds

        # Convert MM1 to QM atom dictionary to absolute indices.
        mm1_to_qm_local = {}
        for k, v in mm1_atoms.items():
            qm_idx = mols.atoms().find(qm_mol.atoms()[k])
            link_idx = mols.atoms().find(qm_mol.atoms()[v])
            # Make sure that we haven't assigned this link atom already.
            if link_idx in mm1_to_qm_local or link_idx in mm1_to_qm:
                raise Exception(
                    f"Cannot substitue the same MM atom (index {link_idx}) "
                    "for more than one link atom!"
                )
            mm1_to_qm_local[link_idx] = qm_idx

        # Convert MM1 to MM2 atom dictionary to absolute indices.
        mm1_to_mm2_local = {}
        for k, v in mm2_atoms.items():
            link_idx = mols.atoms().find(qm_mol.atoms()[k])
            mm_idx = [mols.atoms().find(qm_mol.atoms()[x]) for x in v]
            mm1_to_mm2_local[link_idx] = mm_idx

        # Store the MM1 atom indices.
        mm1_indices.append(list(mm1_atoms.values()))

        # Now work out the QM-MM1 bond distances based on the equilibrium
        # MM bond lengths.

        # A dictionary to store the bond lengths. Here 'bond_lengths' are the
        # equilibrium bond lengths for the QM-MM1 bonds, and 'link_bond_lengths',
        # are the equilibrium bond lengths for the QM-L (QM-link) bonds. Both of
        # these are evaluated from the MM bond potentials.
        bond_lengths = {}
        link_bond_lengths = {}

        # Get the MM bond potentials.
        bonds = qm_mol.property(map["bond"]).potentials()

        # Store the info for the QM molecule.
        info = qm_mol.info()

        # Store the bond potential symbol.
        r = _CAS.Symbol("r")

        # Loop over the link atoms.
        for qm_idx, mm1_idx in mm1_atoms.items():
            # Convert to cg_atom_idx objects.
            cg_qm_idx = info.cg_atom_idx(qm_idx)
            cg_mm1_idx = info.cg_atom_idx(mm1_idx)

            # Store the element of the QM atom.
            qm_elem = qm_mol[cg_qm_idx].element()
            hydrogen = _Mol.Element("H")

            qm_m1_bond_found = False
            qm_link_bond_found = False

            # Loop over the bonds.
            for bond in bonds:
                # Get the indices of the atoms in the bond.
                bond_idx0 = bond.atom0()
                bond_idx1 = bond.atom1()

                # If the bond is between the QM atom and the MM atom, store the
                # bond length.
                if (
                    not qm_m1_bond_found
                    and cg_qm_idx == bond_idx0
                    and cg_mm1_idx == bond_idx1
                    or cg_qm_idx == bond_idx1
                    and cg_mm1_idx == bond_idx0
                ):
                    # Cast as an AmberBond.
                    ab = _MM.AmberBond(bond.function(), r)
                    bond_lengths[mm1_idx] = ab.r0()
                    qm_m1_bond_found = True
                    if qm_link_bond_found:
                        break
                else:
                    elem0 = qm_mol[bond_idx0].element()
                    elem1 = qm_mol[bond_idx1].element()

                    # Is this bond is between a hydrogen and and the same element
                    # as the QM atom? If so, store the bond length.
                    if (
                        not qm_link_bond_found
                        and elem0 == hydrogen
                        and elem1 == qm_elem
                        or elem0 == qm_elem
                        and elem1 == hydrogen
                    ):
                        # Cast as an AmberBond.
                        ab = _MM.AmberBond(bond.function(), r)
                        link_bond_lengths[mm1_idx] = ab.r0()
                        qm_link_bond_found = True
                        if qm_m1_bond_found:
                            break

        # Work out the bond scale factors: R0(QM-L) / R0(QM-MM1)
        try:
            bond_scale_factors_local = {}
            for idx in bond_lengths:
                abs_idx = mols.atoms().find(qm_mol.atoms()[idx])
                bond_scale_factors_local[abs_idx] = (
                    link_bond_lengths[idx] / bond_lengths[idx]
                )
        except:
            raise Exception(
                f"Unable to compute the scaled the QM-MM1 bond lengths for MM1 atom {idx}!"
            )

        # Update the dictionaries.
        mm1_to_qm.update(mm1_to_qm_local)
        mm1_to_mm2.update(mm1_to_mm2_local)
        bond_scale_factors.update(bond_scale_factors_local)

    return mm1_to_qm, mm1_to_mm2, bond_scale_factors, mm1_indices


def _create_merged_mols(qm_mol_to_atoms, mm1_indices, map):
    """
    Internal helper function to create a merged molecule from the QM molecule.

    Parameters
    ----------

    qm_mol_to_atoms: {sire.legacy.Mol.MolNum: [sire.legacy.Mol.AtomIdx]}
        A dictionary with molecule numbers as keys and a list of QM atoms as
        values.

    mm1_indices: [[sire.legacy.Mol.AtomIdx]]
        A list of lists of MM1 atom indices.

    map: sire.legacy.Base.PropertyMap
        The property map for the system.

    Returns
    -------

    qm_mols: [sire.legacy.Mol.Molecule]
        A list of merged molecules.
    """

    from ..legacy import CAS as _CAS
    from ..legacy import Mol as _Mol
    from ..legacy import MM as _MM
    from ..morph import link_to_reference as _link_to_reference

    # Initialise a list to store the merged molecules.
    qm_mols = []

    # Loop over all molecules containing QM atoms.
    for (mol_num, qm_atoms), mm1_idxs in zip(qm_mol_to_atoms.items(), mm1_indices):
        # Get the molecule.
        qm_mol = qm_atoms[0].molecule()

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
                        angles.set(
                            atom0, atom1, atom2, amber_angle.to_expression(theta)
                        )
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

                # Set the charge for all non-QM and non-MM1 atoms to the MM value.
                for atom in qm_mol.atoms():
                    idx = info.atom_idx(atom.index())
                    if idx not in qm_idxs and idx not in mm1_idxs:
                        idx = info.cg_atom_idx(idx)
                        charges.set(idx, atom.property(prop))

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
        qm_mol = _link_to_reference(qm_mol)

        # Add the molecule to the list.
        qm_mols.append(qm_mol)

    # Return the merged molecule.
    return qm_mols


def _configure_engine(engine, mols, qm_atoms, mm1_to_qm, mm1_to_mm2, bond_lengths, map):
    """
    Internal helper function to configure a QM engine ready for dynamics.

    Parameters
    ----------

    engine: sire.legacy.QM.Engine
        The QM engine to configure.

    mols: sire.legacy.System.System
        The Sire system containing the QM atoms.

    qm_atoms: [sire.legacy.Mol.AtomIdx]
        A list of QM atoms.

    mm1_to_qm: {sire.legacy.Mol.AtomIdx: sire.legacy.Mol.AtomIdx}
        A dictionary with link atoms as keys and QM atoms as values.

    mm1_to_mm2: {sire.legacy.Mol.AtomIdx: [sire.legacy.Mol.AtomIdx]}
        A dictionary with link atoms as keys and a list of MM atoms as values.

    bond_lengths: {sire.legacy.Mol.AtomIdx: float}
        A dictionary with link atoms as keys and bond lengths as values.

    map: sire.legacy.Base.PropertyMap
        The property map for the system.
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
        engine.set_link_atoms(mm1_to_qm, mm1_to_mm2, bond_lengths)
    except:
        raise Exception("Unable to set link atom information.")

    return mols, engine
