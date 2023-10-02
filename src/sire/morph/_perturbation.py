__all__ = ["Perturbation"]


class Perturbation:
    """
    This class provides lots of convenience functions that make it
    easier to work with an visualise perturbations
    """

    def __init__(self, mol, map=None):
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

        self._map0 = map.add_suffix("0", props)
        self._map1 = map.add_suffix("1", props)

        # identify all of the ghost atoms
        from_ghosts = []
        to_ghosts = []

        for atom in mol.atoms():
            chg0 = atom.property(self._map0["charge"])
            chg1 = atom.property(self._map1["charge"])

            lj0 = atom.property(self._map0["LJ0"])
            lj1 = atom.property(self._map1["LJ1"])

            is_ghost0 = chg0.is_zero() and lj0.is_dummy()
            is_ghost1 = chg1.is_zero() and lj1.is_dummy()

            if is_ghost0 and not is_ghost1:
                from_ghosts.append(atom.index())
            elif is_ghost1:
                to_ghosts.append(atom.index())

        r = Symbol("r")
        theta = Symbol("theta")

        connectivity = mol.property(self._map0["connectivity"])

        for angle in mol.angles():
            # get the bond lengths desired by the forcefield at the
            # two end states
            in_ring0 = connectivity.in_ring(angle.atom0().index())
            in_ring1 = connectivity.in_ring(angle.atom1().index())
            in_ring2 = connectivity.in_ring(angle.atom2().index())

            if not (in_ring0 and in_ring1 and in_ring2):
                pot0 = AmberAngle(angle.potential(self._map0), theta)
                pot1 = AmberAngle(angle.potential(self._map1), theta)

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
            pot0 = AmberBond(bond.potential(self._map0), r)
            pot1 = AmberBond(bond.potential(self._map1), r)

            r0_0 = pot0.r0()
            r0_1 = pot1.r0()

            self._perturbations.append(
                BondPerturbation(
                    bond=bond.id(),
                    start=r0_0 * angstrom,
                    end=r0_1 * angstrom,
                    map=map,
                )
            )

        self._mol = mol.clone()

    def __str__(self):
        return f"Perturbation( {self._mol} )"

    def __repr__(self):
        return self.__str__()

    def link_to_reference(self, properties: list[str] = None):
        """
        Link all of the properties of the molecule to their values in the
        reference molecule (lambda=0).

        If a list of properties is passed then only those properties will
        be linked to the reference molecule

        Parameters
        ----------

        properties: list[str]
            The list of properties to link to the reference molecule
        """
        if properties is None:
            properties = []

            for key in self._mol.property_keys():
                if key.endswith("0"):
                    properties.append(key[:-1])

        elif type(properties) is str:
            properties = [properties]

        mol = self._mol.molecule().edit()

        for key in properties:
            if mol.has_property(key):
                mol.remove_property(key)

            mol.add_link(key, f"{key}0")

        self._mol.update(mol.commit())
        return self

    def link_to_perturbed(self, properties: list[str] = None):
        """
        Link all of the properties of the molecule to their values in the
        perturbed molecule (lambda=1).

        If a list of properties is passed then only those properties will
        be linked to the perturbed molecule

        Parameters
        ----------

        properties: list[str]
            The list of properties to link to the perturbed molecule
        """
        if properties is None:
            properties = []

            for key in self._mol.property_keys():
                if key.endswith("1"):
                    properties.append(key[:-1])

        elif type(properties) is str:
            properties = [properties]

        mol = self._mol.molecule().edit()

        for key in properties:
            if mol.has_property(key):
                mol.remove_property(key)

            mol.add_link(key, f"{key}1")

        self._mol.update(mol.commit())
        return self

    def set_lambda(self, lam_val: float):
        """
        Set the lambda value to the passed value
        """
        from ..cas import Symbol
        from ..base import create_map
        from ..legacy.CAS import Values

        mol = self._mol.molecule()

        map = create_map({})

        if mol.is_link(map["coordinates"]):
            # make sure we remove any link to the 'coordinates' property
            # as we will be replacing it with the calculated perturbed
            # coordinates
            mol = mol.edit().remove_link(map["coordinates"].source()).commit()

        if not mol.has_property(map["coordinates"]):
            mol = (
                mol.edit()
                .set_property(
                    map["coordinates"].source(),
                    mol.property(self._map0["coordinates"]),
                )
                .commit()
            )

        vals = Values({Symbol("lambda"): lam_val})
        self._mol.update(self._perturbations.perturb(mol, vals))
        return self

    def commit(self):
        """
        Return the modified molecule
        """
        return self._mol.clone()

    def view(self, *args, **kwargs):
        """
        View the perturbation
        """
        from ..cas import Symbol
        from ..base import create_map
        from ..legacy.CAS import Values

        map = create_map({})

        mol = self._mol.clone()
        mol.delete_all_frames(map=map)

        if mol.is_link(map["coordinates"]):
            # make sure we remove any link to the 'coordinates' property
            # as we will be replacing it with the calculated perturbed
            # coordinates
            mol.update(
                mol.molecule()
                .edit()
                .remove_link(map["coordinates"].source())
                .commit()
            )

        if not mol.has_property(map["coordinates"]):
            mol.update(
                mol.molecule()
                .edit()
                .set_property(
                    map["coordinates"].source(),
                    mol.property(self._map0["coordinates"]),
                )
                .commit()
            )

        for lam in range(0, 11):
            vals = Values({Symbol("lambda"): 0.1 * lam})
            mol = self._perturbations.perturb(mol, vals)
            mol.save_frame()

        for lam in range(9, 0, -1):
            vals = Values({Symbol("lambda"): 0.1 * lam})
            mol = self._perturbations.perturb(mol, vals)
            mol.save_frame()

        return mol.view(*args, **kwargs)
