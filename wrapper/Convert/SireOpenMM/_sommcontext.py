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
            from ._SireOpenMM import _set_openmm_coordinates_and_velocities

            map = create_map(map)

            self._atom_index = metadata.index()
            self._lambda_lever = metadata.lambdaLever()

            # we need to update the constraints in the system
            # to match the current value of lambda, before we
            # turn this system into a context
            if map.specified("lambda"):
                lambda_value = map["lambda"].value().as_double()
            elif map.specified("lambda_value"):
                lambda_value = map["lambda_value"].value().as_double()
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
            _set_openmm_coordinates_and_velocities(self, metadata)

            self._lambda_value = self._lambda_lever.set_lambda(
                self, lambda_value
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
            return "none"

    def get_perturbable_constraint(self):
        """
        Return the perturbable constraint applied to the system
        """
        if self._map.specified("perturbable_constraint"):
            return self._map["perturbable_constraint"].source()
        else:
            return "none"

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

        from ._SireOpenMM import _openmm_set_context_platform_property

        _openmm_set_context_platform_property(self, key, value)

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

        return self._lambda_lever.schedule()

    def set_lambda_schedule(self, schedule):
        """
        Set the LambdaSchedule used to control how different forcefield
        parameters will be changed with lambda
        """
        if self._lambda_lever is None:
            raise ValueError("Cannot set the schedule in a null context!")

        self._lambda_lever.set_schedule(schedule)

    def set_lambda(self, lambda_value: float):
        """
        Update the parameters in the context to set the lambda value
        to 'lamval'
        """
        if self._lambda_lever is None:
            return

        self._lambda_value = self._lambda_lever.set_lambda(self, lambda_value)

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

            return (
                nrg.value_in_unit(openmm.unit.kilocalorie_per_mole)
                * kcal_per_mol
            )
        else:
            return nrg

    def get_energy(self, to_sire_units: bool = True):
        """
        Synonym for self.get_potential_energy()
        """
        return self.get_potential_energy(to_sire_units=to_sire_units)
