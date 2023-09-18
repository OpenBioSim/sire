__all__ = ["shrink_ghost_atoms"]


def shrink_ghost_atoms(mols, length=None, map=None):
    """
    Update all of the molecules (or single molecule) in 'mols'
    so that the bond lengths involving ghost atoms at one or other
    end state (but not both) are shrunk to `length`. This is
    0.6 A by default (this seems to be a value that causes
    fewest NaNs).
    """
    from ..base import create_map
    from ..cas import Symbol
    from ..mm import AmberBond

    if length is None:
        length = 0.6
    else:
        from ..units import angstrom, u

        length = u(length).to(angstrom)

    map = create_map(map)

    perturbable_mols = mols.molecules("property is_perturbable")

    props = [
        "LJ",
        "bond",
        "charge",
    ]

    map0 = map.add_suffix("0", props)
    map1 = map.add_suffix("1", props)

    # identify all of the ghost atoms
    from_ghosts = []
    to_ghosts = []

    for mol in perturbable_mols:
        for atom in mol.atoms():
            chg0 = atom.property(map0["charge"])
            chg1 = atom.property(map1["charge"])

            lj0 = atom.property(map0["LJ0"])
            lj1 = atom.property(map1["LJ1"])

            is_ghost0 = chg0.is_zero() and lj0.is_dummy()
            is_ghost1 = chg1.is_zero() and lj1.is_dummy()

            if is_ghost0 and not is_ghost1:
                from_ghosts.append(atom.index())
            elif is_ghost1:
                to_ghosts.append(atom.index())

        r = Symbol("r")

        bonds0 = mol.property(map0["bond"])
        bonds1 = mol.property(map1["bond"])

        changed0 = False
        changed1 = False

        for bond in mol.bonds():
            # are either of the atoms in the bond in from_ghosts or to_ghosts
            is_from_ghost = (
                bond[0].index() in from_ghosts
                or bond[1].index() in from_ghosts
            )

            is_to_ghost = (
                bond[0].index() in to_ghosts or bond[1].index() in to_ghosts
            )

            if is_from_ghost and not is_to_ghost:
                # this is a bond that is turning into a ghost
                r0_0 = length
                pot0 = AmberBond(bond.potential(map0), r)

                if r0_0 < pot0.r0():
                    # only change it if the bond is not already shorter
                    bonds0.set(
                        bond.id(),
                        AmberBond(pot0.k(), r0_0).to_expression(r),
                    )
                    changed0 = True
            elif is_to_ghost:
                r0_1 = length
                pot1 = AmberBond(bond.potential(map1), r)

                if r0_1 < pot1.r0():
                    # only change it if the bond is not already shorter
                    bonds1.set(
                        bond.id(),
                        AmberBond(pot1.k(), r0_1).to_expression(r),
                    )
                    changed1 = True

        if changed0:
            mol = (
                mol.edit().set_property(map0["bond"].source(), bonds0).commit()
            )

        if changed1:
            mol = (
                mol.edit().set_property(map1["bond"].source(), bonds1).commit()
            )

        if changed0 or changed1:
            mols.update(mol)
