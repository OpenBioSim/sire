__all__ = ["Perturbation"]


class Perturbation:
    """
    This class provides lots of convenience functions that make it
    easier to work with an visualise perturbations
    """

    def __init__(self, mol, shrink_ghosts: bool = True, map=None):
        """
        Construct the Perturbation object from the passed molecule.
        Note that this molecule must be perturbable
        """
        from ..base import create_map

        map = create_map(map)

        if not mol.has_property(map["is_perturbable"]):
            raise ValueError(
                "You can only create a `Perturbation` from a "
                "perturbable molecule!"
            )

        if not mol.property(map["is_perturbable"]):
            raise ValueError(
                "You can only create a `Perturbation` from a "
                "perturbable molecule!"
            )

        # construct the perturbation objects that can move the
        # coordinates between the end states
        from ..legacy.Mol import (
            BondPerturbation,
            AnglePerturbation,
            GeometryPerturbations,
        )

        from ..legacy.MM import AmberBond, AmberAngle
        from ..cas import Symbol
        from ..units import angstrom, radian

        self._perturbations = GeometryPerturbations()

        props = [
            "LJ",
            "ambertype",
            "angle",
            "atomtype",
            "bond",
            "charge",
            "coordinates",
            "dihedral",
            "element",
            "forcefield",
            "gb_radii",
            "gb_screening",
            "improper",
            "intrascale",
            "mass",
            "name",
            "parameters",
            "treechain",
        ]

        map0 = map.add_suffix("0", props)
        map1 = map.add_suffix("1", props)

        # identify all of the ghost atoms
        from_ghosts = []
        to_ghosts = []

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
        theta = Symbol("theta")

        if shrink_ghosts:
            bonds0 = mol.property(map0["bond"])
            bonds1 = mol.property(map1["bond"])

            changed0 = False
            changed1 = False

        connectivity = mol.property(map0["connectivity"])

        for angle in mol.angles():
            # get the bond lengths desired by the forcefield at the
            # two end states
            in_ring0 = connectivity.in_ring(angle.atom0().index())
            in_ring1 = connectivity.in_ring(angle.atom1().index())
            in_ring2 = connectivity.in_ring(angle.atom2().index())

            if not (in_ring0 and in_ring1 and in_ring2):
                pot0 = AmberAngle(angle.potential(map0), theta)
                pot1 = AmberAngle(angle.potential(map1), theta)

                theta0_0 = pot0.theta0()
                theta0_1 = pot1.theta0()

                self._perturbations.append(
                    AnglePerturbation(
                        angle=angle.id(),
                        start=theta0_0 * radian,
                        end=theta0_1 * radian,
                        map=map,
                    )
                )

        for bond in mol.bonds():
            # get the bond lengths desired by the forcefield at the
            # two end states
            pot0 = AmberBond(bond.potential(map0), r)
            pot1 = AmberBond(bond.potential(map1), r)

            r0_0 = pot0.r0()
            r0_1 = pot1.r0()

            if shrink_ghosts:
                is_from_ghost = (
                    bond[0].index() in from_ghosts
                    or bond[1].index() in from_ghosts
                )

                is_to_ghost = (
                    bond[0].index() in to_ghosts
                    or bond[1].index() in to_ghosts
                )

                if is_from_ghost and not is_to_ghost:
                    r0_0 = 0.2
                    bonds0.set(
                        bond.id(),
                        AmberBond(pot0.k(), r0_0).to_expression(r),
                    )
                    changed0 = True
                elif is_to_ghost:
                    r0_1 = 0.2
                    bonds1.set(
                        bond.id(),
                        AmberBond(pot1.k(), r0_1).to_expression(r),
                    )
                    changed1 = True

            self._perturbations.append(
                BondPerturbation(
                    bond=bond.id(),
                    start=r0_0 * angstrom,
                    end=r0_1 * angstrom,
                    map=map,
                )
            )

        if shrink_ghosts:
            if changed0:
                mol = (
                    mol.edit()
                    .set_property(map0["bond"].source(), bonds0)
                    .commit()
                )

            if changed1:
                mol = (
                    mol.edit()
                    .set_property(map1["bond"].source(), bonds1)
                    .commit()
                )

        self._mol = (
            mol.edit()
            .set_property(
                map["coordinates"].source(), mol.property(map0["coordinates"])
            )
            .commit()
        )

    def __str__(self):
        return f"Perturbation( {self._mol} )"

    def __repr__(self):
        return self.__str__()

    def set_lambda(self, lam_val: float):
        """
        Set the lambda value to the passed value
        """
        from ..cas import Symbol
        from ..legacy.CAS import Values

        vals = Values({Symbol("lambda"): lam_val})
        self._mol = self._perturbations.perturb(self._mol, vals)

    def view(self, *args, **kwargs):
        """
        View the perturbation
        """
        from ..cas import Symbol
        from ..legacy.CAS import Values

        mol = self._mol.clone()
        mol.delete_all_frames()

        for lam in range(0, 11):
            vals = Values({Symbol("lambda"): 0.1 * lam})
            mol = self._perturbations.perturb(mol, vals)
            mol.save_frame()

        for lam in range(9, 0, -1):
            vals = Values({Symbol("lambda"): 0.1 * lam})
            mol = self._perturbations.perturb(mol, vals)
            mol.save_frame()

        return mol.view(*args, **kwargs)
