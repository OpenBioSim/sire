__all__ = ["SOMMContext"]

from openmm import Context as _Context


class SOMMContext(_Context):
    """
    This is a specialised version of an OpenMM Context
    class that provides extra information and functions
    to support the integration with sire.

    Extra information includes an atom index (so
    we can map atom metadata to atom index in the
    openmm context) and a LambdaLever that can be used
    to morph a context between end states. There are
    extra functions that use this LambdaLever to
    set the lambda value of the context
    """

    def __init__(
        self,
        system=None,
        integrator=None,
        platform=None,
        metadata=None,
        map=None,
    ):
        """
        Construct from a passed OpenMM Context, the
        atom index, and the lambda lever
        """
        if system is not None:
            from ...base import create_map
            from ._SireOpenMM import set_openmm_coordinates_and_velocities

            map = create_map(map)

            self._atom_index = metadata.index()
            self._lambda_lever = metadata.lambda_lever()

            # we need to update the constraints in the system
            # to match the current value of lambda, before we
            # turn this system into a context
            if map.specified("lambda"):
                lambda_value = map["lambda"].value().as_double()
            elif map.specified("lambda_value"):
                lambda_value = map["lambda_value"].value().as_double()
            elif map.specified("qm_engine"):
                # Default to full QM.
                lambda_value = 1.0
            else:
                lambda_value = 0.0

            # we have space here, if we need it, to update
            # the system based on the requested value of lambda
            # There are some things that are unchangeable once
            # the context has been created (e.g. constraints)
            #
            # Note the initialisation of the constraints is slow,
            # so if there are many constraints then this constructor
            # can take several seconds to complete. However, the
            # increase in dynamics performance through having a larger
            # timestep if we constraint bonds and angles makes this
            # more than worth it. Note that this is why minimisations
            # start running more quickly than dynamics jobs. There
            # are no constraints in minimisations
            super().__init__(system, integrator, platform)

            # place the coordinates and velocities into the context
            set_openmm_coordinates_and_velocities(self, metadata)

            # Check for a REST2 scaling factor.
            if map.specified("rest2_scale"):
                try:
                    rest2_scale = map["rest2_scale"].value().as_double()
                except:
                    raise ValueError("'rest2_scale' must be of type 'float'")
                if rest2_scale < 1.0:
                    raise ValueError("'rest2_scale' must be >= 1.0")
                self._rest2_scale = rest2_scale

                if map.specified("rest2_selection"):
                    self._has_rest2_selection = True
                else:
                    self._has_rest2_selection = False
            else:
                self._rest2_scale = 1.0
                self._has_rest2_selection = False

            self._lambda_value = self._lambda_lever.set_lambda(
                self,
                lambda_value=lambda_value,
                rest2_scale=self._rest2_scale,
                update_constraints=True,
            )

            self._map = map
        else:
            self._atom_index = None
            self._lambda_lever = None
            self._lambda_value = 0.0
            self._map = None

    def __str__(self):
        p = self.getPlatform()
        s = self.getSystem()
        i = self.getIntegrator()
        return (
            f"openmm::Context( num_atoms={s.getNumParticles()} "
            f"integrator={i.__class__.__name__} "
            f"timestep={i.getStepSize()._value * 1000} fs "
            f"platform={p.getName()} )"
        )

    def __repr__(self):
        return self.__str__()

    def get_constraint(self):
        """
        Return the constraint applied to the system
        """
        if self._map.specified("constraint"):
            return self._map["constraint"].source()
        else:
            return None

    def get_perturbable_constraint(self):
        """
        Return the perturbable constraint applied to the system
        """
        if self._map.specified("perturbable_constraint"):
            return self._map["perturbable_constraint"].source()
        else:
            return None

    def get_platform(self):
        """
        Return the platform used for this simulation
        """
        return self.getPlatform()

    def get_integrator(self):
        """
        Return the integrator used for the simulation
        """
        return self.getIntegrator()

    def get_platform_properties(self):
        """
        Return all of the properties (and their values) of the
        platform used to run the simulation.
        """
        props = {}

        platform = self.getPlatform()

        for key in platform.getPropertyNames():
            props[key] = platform.getPropertyValue(self, key)

        return props

    def get_platform_property(self, key):
        """
        Return the value of the specified platform property
        """
        platform = self.getPlatform()

        keys = platform.getPropertyNames()

        if key not in keys:
            keys = ", ".join(keys)

            raise KeyError(
                f"There is no platform property called {key} in the "
                f"platform {platform.getName()}. Available properties "
                f"are [ {keys} ]"
            )

        return platform.getPropertyValue(self, key)

    def set_platform_property(self, key, value):
        """
        Set the value of the platform property 'key' to 'value'

        Note that this probably doesn't do anything as it looks
        like platform properties cannot be changed after
        construction. This function may be removed.
        """
        value = str(value)

        platform = self.getPlatform()

        keys = platform.getPropertyNames()

        if key not in keys:
            keys = ", ".join(keys)

            raise KeyError(
                f"There is no platform property called {key} in the "
                f"platform {platform.getName()}. Available properties "
                f"are [ {keys} ]"
            )

        from ._SireOpenMM import set_context_platform_property

        set_context_platform_property(self, key, value)

    def get_atom_index(self):
        """
        Return the index mapping from atom metadata to its index
        in the openmm context
        """
        return self._atom_index

    def get_lambda(self):
        """
        Return the current value of lambda for this context
        """
        return self._lambda_value

    def get_lambda_lever(self):
        """
        Return the LambdaLever used to change the lambda value
        of the parameters
        """
        return self._lambda_lever

    def get_lambda_schedule(self):
        """
        Return the LambdaSchedule used to control how different forcefield
        parameters will be changed with lambda
        """
        if self._lambda_lever is None:
            return None

        return self._lambda_lever.get_schedule()

    def set_lambda_schedule(self, schedule):
        """
        Set the LambdaSchedule used to control how different forcefield
        parameters will be changed with lambda
        """
        if self._lambda_lever is None:
            raise ValueError("Cannot set the schedule in a null context!")

        self._lambda_lever.set_schedule(schedule)

    def set_lambda(
        self,
        lambda_value: float,
        rest2_scale: float = None,
        update_constraints: bool = True,
    ):
        """
        Update the parameters in the context to set the lambda value
        to 'lambda_value'. The 'rest2_scale' defines the temperature of
        the REST2 region relative to the rest of the system. If
        'update_constraints' is True then the constraints will be updated
        to match the new value of lambda.
        """
        if self._lambda_lever is None:
            return

        # If not provided, use the REST2 scaling factor used to initalise
        # the context.
        if rest2_scale is None:
            rest2_scale = self._rest2_scale
        else:
            if rest2_scale < 1.0:
                raise ValueError("'rest2_scale' must be >= 1.0")

        if (lambda_value == self._lambda_value) and (rest2_scale == self._rest2_scale):
            # Nothing to do.
            return

        self._lambda_value = self._lambda_lever.set_lambda(
            self,
            lambda_value=lambda_value,
            rest2_scale=rest2_scale,
            update_constraints=update_constraints,
        )

        # Update any additional parameters in the REST2 region.
        if self._has_rest2_selection and rest2_scale != self._rest2_scale:
            self._update_rest2(lambda_value, rest2_scale)
            self._rest2_scale = rest2_scale

    def get_rest2_scale(self):
        """
        Return the temperature scale factor for the REST2 region.
        """
        return self._rest2_scale

    def set_rest2_scale(self, rest2_scale):
        """
        Set the temperature scale factor for the REST2 region.
        """
        self._set_lambda(self._lambda_value, rest2_scale=rest2_scale)

    def set_temperature(self, temperature, rescale_velocities=True):
        """
        Set the target temperature for the dynamics. If
        rescale_velocities is True then the velocities will
        be rescaled to the new temperature
        """
        raise NotImplementedError("We can't yet set the temperature")

    def set_pressure(self, pressure):
        """
        Set the target pressure for the dynamics.
        """
        raise NotImplementedError("We can't yet set the pressure")

    def get_potential_energy(self, to_sire_units: bool = True):
        """
        Calculate and return the potential energy of the system
        """
        s = self.getState(getEnergy=True)
        nrg = s.getPotentialEnergy()

        if to_sire_units:
            import openmm
            from ...units import kcal_per_mol

            return nrg.value_in_unit(openmm.unit.kilocalorie_per_mole) * kcal_per_mol
        else:
            return nrg

    def get_energy(self, to_sire_units: bool = True):
        """
        Synonym for self.get_potential_energy()
        """
        return self.get_potential_energy(to_sire_units=to_sire_units)

    def get_constraints(self):
        """
        Return all pairs of atoms that are constrained, together with
        the constraint distance
        """
        s = self.getSystem()

        num_constraints = s.getNumConstraints()

        constraints = []

        import openmm
        from ...units import nanometer

        for i in range(num_constraints):
            a1, a2, dist = s.getConstraintParameters(i)

            constraints.append(
                (
                    self._atom_index[a1] + self._atom_index[a2],
                    dist.value_in_unit(openmm.unit.nanometer) * nanometer,
                )
            )

        return constraints

    def to_xml(self, f=None):
        """
        Save the current state of the dynamics to XML.
        This is mostly used for debugging. This will return the
        XML string if 'f' is None. Otherwise it will write the
        XML to 'f' (either a filename, or a FILE object)
        """
        from openmm.openmm import XmlSerializer as _XmlSerializer

        if f is None:
            return _XmlSerializer.serialize(self.getSystem())
        elif isinstance(f, str):
            with open(f, "w") as handle:
                handle.write(_XmlSerializer.serialize(self.getSystem()))
        else:
            f.write(_XmlSerializer.serialize(self.getSystem()))

    def _prepare_rest2(self, system, atoms):
        """
        Internal method to prepare the REST2 data structures.
        """

        # Adapted from code in meld: https://github.com/maccallumlab/meld

        import openmm
        from ..Mol import AtomIdx, Connectivity

        # Work out the molecules to which the atoms belong.
        mols = []
        for atom in atoms:
            mol = atom.molecule()
            # Perturbable molecules are handled separately.
            if mol not in mols and not mol.has_property("is_perturbable"):
                mols.append(mol)

        # Store the OpenMM system.
        omm_system = self.getSystem()

        # Get the NonBonded force.
        nonbonded_force = None
        for force in omm_system.getForces():
            if isinstance(force, openmm.NonbondedForce):
                nonbonded_force = force
                break
        if nonbonded_force is None:
            raise ValueError("No NonbondedForce found in the OpenMM system.")
        self._nonbonded_force = nonbonded_force

        # Get the PeriodicTorsionForce.
        periodic_torsion_force = None
        for force in omm_system.getForces():
            if isinstance(force, openmm.PeriodicTorsionForce):
                periodic_torsion_force = force
                break
        if periodic_torsion_force is None:
            raise ValueError("No PeriodicTorsionForce found in the OpenMM system.")
        self._periodic_torsion_force = periodic_torsion_force

        # Initialise the parameter dictionaries.
        self._nonbonded_params = {}
        self._exception_params = {}
        self._torsion_params = {}

        # Store the molecules in the system.
        system_mols = system.molecules()

        # Process each of the molecules.
        for mol in mols:
            # Create the connectivity object for the molecule.
            connectivity = Connectivity(mol.info()).edit()

            # Loop over the bonds in the molecule and connect the atoms.
            for bond in mol.bonds():
                connectivity.connect(bond.atom0().index(), bond.atom1().index())
            connectivity = connectivity.commit()

            # Find the index of the molecule in the system.
            mol_idx = system._system.mol_nums().index(mol.number())

            # Work out the offset to apply to the atom indices to convert to system indices.
            num_atoms = 0
            for i in range(mol_idx):
                num_atoms += system_mols[i].num_atoms()

            # Create a list of atom indices.
            atom_idxs = [atom.index().value() + num_atoms for atom in atoms]

            # Gather the nonbonded parameters for the atoms in the selection.
            for idx in atom_idxs:
                self._nonbonded_params[idx] = nonbonded_force.getParticleParameters(idx)

            # Store the exception parameters.
            for param_index in range(nonbonded_force.getNumExceptions()):
                params = nonbonded_force.getExceptionParameters(param_index)
                if params[0] in atom_idxs and params[1] in atom_idxs:
                    self._exception_params[param_index] = params

            # Gather the torsion parameters for the atoms in the selection.
            for param_index in range(periodic_torsion_force.getNumTorsions()):
                params = periodic_torsion_force.getTorsionParameters(param_index)
                i, j, k, l, _, _, _ = params

                # Don't modify non-REST2 torsions.
                if (
                    i not in atom_idxs
                    or j not in atom_idxs
                    or k not in atom_idxs
                    or l not in atom_idxs
                ):
                    continue

                # Convert to AtomIdx objects.
                idx_i = AtomIdx(i - num_atoms)
                idx_l = AtomIdx(l - num_atoms)

                if connectivity.are_dihedraled(idx_i, idx_l):
                    self._torsion_params[param_index] = params

    def _update_rest2(self, lambda_value, rest2_scale):
        """
        Internal method to update the REST2 parameters.
        """

        from math import sqrt

        # This is the temperature scale factor, so we need to invert to get the energy
        # scale factor.
        scale = 1.0 / rest2_scale

        # Store the REST2 charge scaling factor for non-bonded interactions.
        sqrt_scale = sqrt(scale)

        # Update the non-bonded parameters.
        for index, params in self._nonbonded_params.items():
            q, sigma, epsilon = params
            self._nonbonded_force.setParticleParameters(
                index, q * sqrt_scale, sigma, epsilon * scale
            )

        # Update the exception parameters.
        for index, params in self._exception_params.items():
            i, j, q, sigma, epsilon = params
            self._nonbonded_force.setExceptionParameters(
                index, i, j, q * scale, sigma, epsilon * scale
            )

        # Update the parameters in the context.
        self._nonbonded_force.updateParametersInContext(self)

        # Update the torsion parameters.
        for index, params in self._torsion_params.items():
            i, j, k, l, periodicity, phase, fc = params
            self._periodic_torsion_force.setTorsionParameters(
                index, i, j, k, l, periodicity, phase, fc * scale
            )

        # Update the parameters in the context.
        self._periodic_torsion_force.updateParametersInContext(self)
