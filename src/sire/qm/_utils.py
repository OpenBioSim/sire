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
