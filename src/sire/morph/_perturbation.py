__all__ = [
    "Perturbation",
    "link_to_reference",
    "link_to_perturbed",
    "extract_reference",
    "extract_perturbed",
    "zero_ghost_torsions",
]


def link_to_reference(mols, map=None):
    """
    Return the passed molecule(s), where any perturbable molecules
    are linked to their reference (lambda=0) states
    """
    mols = mols.clone()

    for mol in mols.molecules("property is_perturbable"):
        mols.update(mol.perturbation(map=map).link_to_reference(auto_commit=True))

    return mols


def link_to_perturbed(mols, map=None):
    """
    Return the passed molecule(s), where any perturbable molecules
    are linked to their perturbed (lambda=1) states
    """
    mols = mols.clone()

    for mol in mols.molecules("property is_perturbable"):
        mols.update(mol.perturbation(map=map).link_to_perturbed(auto_commit=True))

    return mols


def extract_reference(mols, remove_ghosts: bool = True, map=None):
    """
    Return the passed molecule(s) where any perturbable molecules
    have been extracted to their reference (lambda=0) state (i.e. deleting
    their perturbed state)
    """
    mols = mols.clone()

    for mol in mols.molecules("property is_perturbable"):
        mols.update(
            mol.perturbation(map=map).extract_reference(remove_ghosts=remove_ghosts)
        )

    return mols


def extract_perturbed(mols, remove_ghosts: bool = True, map=None):
    """
    Return the passed molecule(s) where any perturbable molecules
    have been extracted to their perturbed (lambda=1) state (i.e. deleting
    their reference state)
    """
    mols = mols.clone()

    for mol in mols.molecules("property is_perturbable"):
        mols.update(
            mol.perturbation(map=map).extract_perturbed(remove_ghosts=remove_ghosts)
        )

    return mols


def zero_ghost_torsions(mols, map=None):
    """
    Zero any torsions (dihedrals or impropers) in the reference or perturbed
    states where any of the atoms in those states is a ghost (dummy) atom
    """
    mols = mols.clone()

    for mol in mols.molecules("property is_perturbable"):
        mols.update(mol.perturbation(map=map).zero_ghost_torsions(auto_commit=True))

    return mols


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
                "You can only create a `Perturbation` from a " "perturbable molecule!"
            )

        if not mol.property(map["is_perturbable"]):
            raise ValueError(
                "You can only create a `Perturbation` from a " "perturbable molecule!"
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

        self._map = map
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

    def extract_reference(
        self,
        properties: list[str] = None,
        remove_ghosts: bool = True,
    ):
        """
        Extract the reference state properties of the molecule
        and remove the perturbed state.

        If a list of properties is passed then only those properties will
        be extracted to the reference molecule

        Parameters
        ----------

        properties: list[str]
            The list of properties to extract to the reference molecule

        remove_ghosts: bool
            If True then any ghost atoms will be removed from the molecule
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
            ref_prop = f"{key}0"
            pert_prop = f"{key}1"

            if mol.has_property(ref_prop):
                if mol.has_property(key):
                    mol.remove_property(key)

                mol.set_property(key, mol.property(ref_prop))
                mol.remove_property(ref_prop)

                if mol.has_property(pert_prop):
                    mol.remove_property(pert_prop)

        if mol.has_property("is_perturbable"):
            mol.remove_property("is_perturbable")

        mol = mol.commit().molecule()

        if remove_ghosts:
            mol = mol["not element Xx"].extract(to_same_molecule=True)

        return mol

    def extract_perturbed(
        self,
        properties: list[str] = None,
        remove_ghosts: bool = True,
    ):
        """
        Extract the perturbed state properties of the molecule
        and remove the reference state.

        If a list of properties is passed then only those properties will
        be extracted to the reference molecule

        Parameters
        ----------

        properties: list[str]
            The list of properties to extract to the perturbed molecule

        remove_ghosts: bool
            If True then any ghost atoms will be removed from the molecule
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
            ref_prop = f"{key}0"
            pert_prop = f"{key}1"

            if mol.has_property(pert_prop):
                if mol.has_property(key):
                    mol.remove_property(key)

                mol.set_property(key, mol.property(pert_prop))
                mol.remove_property(pert_prop)

                if mol.has_property(ref_prop):
                    mol.remove_property(ref_prop)

        if mol.has_property("is_perturbable"):
            mol.remove_property("is_perturbable")

        mol = mol.commit().molecule()

        if remove_ghosts:
            mol = mol["not element Xx"].extract(to_same_molecule=True)

        return mol

    def link_to_reference(self, properties: list[str] = None, auto_commit: bool = True):
        """
        Link all of the properties of the molecule to their values in the
        reference molecule (lambda=0).

        If a list of properties is passed then only those properties will
        be linked to the reference molecule

        Parameters
        ----------

        properties: list[str]
            The list of properties to link to the reference molecule

        auto_commit: bool
            If True then the molecule will be committed and returned
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

        if auto_commit:
            return self.commit()
        else:
            return self

    def link_to_perturbed(self, properties: list[str] = None, auto_commit: bool = True):
        """
        Link all of the properties of the molecule to their values in the
        perturbed molecule (lambda=1).

        If a list of properties is passed then only those properties will
        be linked to the perturbed molecule

        Parameters
        ----------

        properties: list[str]
            The list of properties to link to the perturbed molecule

        auto_commit: bool
            If True then the molecule will be committed and returned
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

        if auto_commit:
            return self.commit()
        else:
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

    def view(self, *args, state="perturbed", **kwargs):
        """
        View the perturbation
        """
        from ..cas import Symbol
        from ..base import create_map
        from ..legacy.CAS import Values

        map = create_map({})

        if str(state).lower() == "perturbed":
            mol = self._mol.clone().perturbation(map=self._map).link_to_perturbed()
        else:
            mol = self._mol.clone().perturbation(map=self._map).link_to_reference()

        mol.delete_all_frames(map=map)

        if mol.is_link(map["coordinates"]):
            # make sure we remove any link to the 'coordinates' property
            # as we will be replacing it with the calculated perturbed
            # coordinates
            mol.update(
                mol.molecule().edit().remove_link(map["coordinates"].source()).commit()
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

        return mol["not element Xx"].view(*args, **kwargs)

    def view_reference(self, *args, **kwargs):
        """
        View the reference state of the perturbation
        """
        return self.view(*args, state="reference", **kwargs)

    def view_perturbed(self, *args, **kwargs):
        """
        View the perturbed state of the perturbation
        """
        return self.view(*args, state="perturbed", **kwargs)

    def zero_ghost_torsions(self, auto_commit: bool = True):
        """
        Zero the torsions (dihedrals and impropers) in the reference and
        perturbed states where any of the atoms in those states is a ghost atom
        """
        p = self.to_openmm()

        zero_in_ref = []
        zero_in_pert = []

        from_ghosts = []
        to_ghosts = []

        atoms = p.atoms()

        for idx in p.get_to_ghost_idxs():
            to_ghosts.append(atoms[idx].index())

        for idx in p.get_from_ghost_idxs():
            from_ghosts.append(atoms[idx].index())

        for torsion in p.torsions():
            if (
                torsion.atom0().index() in from_ghosts
                or torsion.atom1().index() in from_ghosts
                or torsion.atom2().index() in from_ghosts
                or torsion.atom3().index() in from_ghosts
            ):
                zero_in_ref.append(torsion)

            if (
                torsion.atom0().index() in to_ghosts
                or torsion.atom1().index() in to_ghosts
                or torsion.atom2().index() in to_ghosts
                or torsion.atom3().index() in to_ghosts
            ):
                zero_in_pert.append(torsion)

        if len(zero_in_ref) == 0 and len(zero_in_pert) == 0:
            # nothing to do
            return self

        mol = self._mol.molecule().edit()

        if len(zero_in_ref) > 0:
            dihs = self._mol.property(self._map["dihedral0"])
            imps = self._mol.property(self._map["improper0"])

            for torsion in zero_in_ref:
                dihs.clear(torsion.id())
                imps.clear(torsion.id())

            mol.set_property(self._map["dihedral0"].source(), dihs)
            mol.set_property(self._map["improper0"].source(), imps)

        if len(zero_in_pert) > 0:
            dihs = self._mol.property(self._map["dihedral1"])
            imps = self._mol.property(self._map["improper1"])

            for torsion in zero_in_pert:
                dihs.clear(torsion.id())
                imps.clear(torsion.id())

            mol.set_property(self._map["dihedral1"].source(), dihs)
            mol.set_property(self._map["improper1"].source(), imps)

        self._mol.update(mol.commit())

        if auto_commit:
            return self.commit()
        else:
            return self

    def to_openmm(
        self,
        constraint: str = None,
        perturbable_constraint: str = None,
        swap_end_states: bool = None,
        include_constrained_energies: bool = None,
        map: dict = None,
    ):
        """
        Return the PerturbableOpenMMMolecule that corresponds to
        this perturbation in the OpenMM simulation. The arguments
        to this function have the same meaning as the equivalents
        in the dynamics() and minimisation() functions of a molecule.

        Parameters
        ----------

        constraint: str
            The constraint algorithm to use

        perturbable_constraint: str
            The constraint algorithm to use for perturbable atoms

        swap_end_states: bool
            If True then the end states will be swapped

        include_constrained_energies: bool
            If True then the constrained energies will be included

        map: dict
            The property map to use


        Returns
        -------
        PerturbableOpenMMMolecule
            The perturbable OpenMM molecule
        """
        from ..base import create_map

        map = create_map(self._map, map)

        if constraint is not None:
            map.set("constraint", str(constraint))

        if perturbable_constraint is not None:
            map.set("perturbable_constraint", str(perturbable_constraint))

        if swap_end_states is not None:
            map.set("swap_end_states", bool(swap_end_states))

        if include_constrained_energies is not None:
            map.set("include_constrained_energies", bool(include_constrained_energies))

        from ..convert.openmm import PerturbableOpenMMMolecule

        try:
            return PerturbableOpenMMMolecule(self._mol.molecule(), map=map)
        except KeyError:
            # probably need to choose an end-state - default to reference
            return PerturbableOpenMMMolecule(
                self._mol.perturbation().link_to_reference(auto_commit=True).molecule(),
                map=map,
            )
