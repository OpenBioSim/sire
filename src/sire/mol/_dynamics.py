__all__ = ["Dynamics"]


class DynamicsData:
    """
    Internal class that is designed to only be used by the Dynamics
    class. This holds the shared state for dynamics on a set
    of molecule(s).
    """

    def __init__(self, mols=None, map=None, **kwargs):
        from ..base import create_map

        map = create_map(map, kwargs)

        # Save the original view, so that we can return an object
        # of the same type in .commit()
        self._orig_mols = mols

        if mols is not None:
            self._map = map  # this is already a PropertyMap

            # see if there are any fixed atoms
            if map.specified("fixed"):
                if map["fixed"].has_value():
                    fixed_atoms = map["fixed"].value()
                else:
                    fixed_atoms = map["fixed"].source()

                from . import selection_to_atoms

                # turn the fixed atoms into a list of atoms
                self._map.set(
                    "fixed",
                    mols.atoms().find(selection_to_atoms(mols, fixed_atoms)),
                )

            # see if there is a QM/MM engine
            if map.specified("qm_engine"):
                qm_engine = map["qm_engine"].value()

                from ..legacy.Convert import QMEngine
                from warnings import warn

                if qm_engine and not isinstance(qm_engine, QMEngine):
                    raise ValueError(
                        "'qm_engine' must be an instance of 'sire.legacy.Convert.QMEngine'"
                    )

                # If a QM/MM engine is specified, then we need to check that there is a
                # perturbable molecule.
                try:
                    pert_mol = mols["property is_perturbable"]
                except:
                    raise ValueError(
                        "You are trying to run QM/MM dynamics for a system without "
                        "a QM/MM enabled molecule!"
                    )

                # Check the constraints and raise a warning if the perturbable_constraint
                # is not "none".

                if map.specified("perturbable_constraint"):
                    perturbable_constraint = map["perturbable_constraint"].source()
                    if perturbable_constraint.lower() != "none":
                        warn(
                            "Running a QM/MM simulation with constraints on the QM "
                            "region is not recommended."
                        )
                else:
                    # The perturbable constraint is unset, so will follow the constraint.
                    # Make sure this is "none".
                    if map.specified("constraint"):
                        constraint = map["constraint"].source()
                        if constraint.lower() != "none":
                            warn(
                                "Running a QM/MM simulation with constraints on the QM "
                                "region is not recommended."
                            )
                    # Constraints will be automatically applied, so we can't guarantee that
                    # the constraint is "none".
                    else:
                        warn(
                            "Running a QM/MM simulation with constraints on the QM "
                            "region is not recommended."
                        )

            if map.specified("cutoff"):
                cutoff = map["cutoff"]

                if cutoff.has_source():
                    cutoff = cutoff.source()

                    if cutoff.lower() == "none" or cutoff.lower().startswith("infinit"):
                        self._map.set("cutoff_type", "NONE")
                        self._map.unset("cutoff")
                    elif cutoff.lower() == "auto":
                        self._map.unset("cutoff")
                    elif cutoff != "cutoff":
                        from .. import u

                        self._map.set("cutoff", u(cutoff))

            # get the forcefield info from the passed parameters
            # and from whatever we can glean from the molecules
            from ..system import ForceFieldInfo, System

            self._ffinfo = ForceFieldInfo(mols, map=self._map)

            # We want to store the molecules as a System so that
            # we can more easily track the space and time properties
            if type(mols) is System:
                # work on our own copy of the system
                self._sire_mols = mols.clone()
            else:
                # create a system to work on
                self._sire_mols = System()
                self._sire_mols._system.add(mols.molecules().to_molecule_group())
                self._sire_mols._system.set_property("space", self._ffinfo.space())

            # see if this is an interpolation simulation
            if map.specified("lambda_interpolate"):
                if map["lambda_interpolate"].has_value():
                    lambda_interpolate = map["lambda_interpolate"].value()
                else:
                    lambda_interpolate = map["lambda_interpolate"].source()

                # Single lambda value.
                try:
                    lambda_interpolate = float(lambda_interpolate)
                    map.set("lambda_value", lambda_interpolate)
                # Two lambda values.
                except:
                    try:
                        if not len(lambda_interpolate) == 2:
                            raise
                        lambda_interpolate = [float(x) for x in lambda_interpolate]
                        map.set("lambda_value", lambda_interpolate[0])
                    except:
                        raise ValueError(
                            "'lambda_interpolate' must be a float or a list of two floats"
                        )

                from ..units import kcal_per_mol

                self._is_interpolate = True
                self._lambda_interpolate = lambda_interpolate
                self._work = 0 * kcal_per_mol
                self._nrg_prev = 0 * kcal_per_mol

            else:
                self._is_interpolate = False

            # find the existing energy trajectory - we will build on this
            self._energy_trajectory = self._sire_mols.energy_trajectory(
                to_pandas=False, map=self._map
            )

            # make sure that the ensemble is recorded in the trajectory
            self._energy_trajectory.set_property("ensemble", self.ensemble())

            self._num_atoms = len(self._sire_mols.atoms())

            from ..units import nanosecond

            try:
                current_time = (
                    self._sire_mols.property(self._map["time"]).to(nanosecond)
                    * nanosecond
                )
            except Exception:
                current_time = 0 * nanosecond
                self._sire_mols.set_property("time", current_time)

            self._current_time = current_time
            self._current_step = 0
            self._prev_step = 0
            self._elapsed_time = 0 * nanosecond
            self._prev_current_time = 0 * nanosecond
            self._prev_elapsed_time = 0 * nanosecond
            self._walltime = 0 * nanosecond
            self._is_running = False
            self._schedule_changed = False

            # Initialise the GCMC sampler. This will be updated externally.
            # if the dynamics object is coupled to a sampler.
            self._gcmc_sampler = None

            # Check for a REST2 scaling factor.
            if map.specified("rest2_scale"):
                try:
                    rest2_scale = map["rest2_scale"].value().as_double()
                except:
                    raise ValueError("'rest2_scale' must be of type 'float'")
                if rest2_scale < 1.0:
                    raise ValueError("'rest2_scale' must be >= 1.0")
            else:
                rest2_scale = 1.0

            # Check for an additional REST2 selection.
            if map.specified("rest2_selection"):
                try:
                    rest2_selection = str(map["rest2_selection"])
                except:
                    raise ValueError("'rest2_selection' must be of type 'str'")

                try:
                    from . import selection_to_atoms

                    # Try to find the REST2 selection.
                    atoms = selection_to_atoms(mols, rest2_selection)
                except:
                    raise ValueError(
                        "Invalid 'rest2_selection' string. Please check the selection syntax."
                    )

                # Store all the perturbable molecules associated with the selection
                # and remove perturbable atoms from the selection. Remove alchemical ions
                # from the selection.
                pert_mols = {}
                non_pert_atoms = atoms.to_list()
                for atom in atoms:
                    mol = atom.molecule()
                    if mol.has_property("is_alchemical_ion"):
                        non_pert_atoms.remove(atom)
                    elif mol.has_property("is_perturbable"):
                        non_pert_atoms.remove(atom)
                        if mol.number() not in pert_mols:
                            pert_mols[mol.number()] = [atom]
                        else:
                            pert_mols[mol.number()].append(atom)

                # Now create a boolean is_rest2 mask for the atoms in the perturbable molecules.
                # Only do this if there are perturbable atoms in the selection.
                if len(non_pert_atoms) != len(atoms):
                    for num in pert_mols:
                        mol = self._sire_mols[num]
                        is_rest2 = [False] * mol.num_atoms()
                        for atom in pert_mols[num]:
                            is_rest2[atom.index().value()] = True

                        # Set the is_rest2 property for each perturbable molecule.
                        mol = (
                            mol.edit()
                            .set_property("is_rest2", is_rest2)
                            .molecule()
                            .commit()
                        )

                        # Update the system.
                        self._sire_mols.update(mol)

            # Search for alchemical ions and exclude them via a REST2 mask.
            try:
                for mol in self._sire_mols.molecules("property is_alchemical_ion"):
                    is_rest2 = [False] * mol.num_atoms()
                    mol = (
                        mol.edit()
                        .set_property("is_rest2", is_rest2)
                        .molecule()
                        .commit()
                    )
                    self._sire_mols.update(mol)
            except:
                pass

            from ..convert import to

            self._omm_mols = to(self._sire_mols, "openmm", map=self._map)
            self._clear_state()

            if self._ffinfo.space().is_periodic():
                self._enforce_periodic_box = True

                if self._map.specified("enforce_periodic_box"):
                    self._enforce_periodic_box = (
                        self._map["enforce_periodic_box"].value().as_boolean()
                    )
            else:
                self._enforce_periodic_box = False

            # Prepare the OpenMM REST2 data structures.
            if map.specified("rest2_selection"):
                if len(non_pert_atoms) > 0:
                    self._omm_mols._prepare_rest2(self._sire_mols, non_pert_atoms)
        else:
            self._sire_mols = None
            self._energy_trajectory = None

        # Store the pressure value needed for the contribution to the reduced potential.
        if self._map.specified("pressure"):
            try:
                from sire import u
                from sire.units import mole

                NA = 6.02214076e23 / mole
                pressure = u(self._map["pressure"].source())
                self._pressure = (pressure * NA).value()
            except Exception as e:
                raise ValueError(
                    "'Unable to computed reduced pressure from 'pressure' "
                    f"value: {self._map['pressure']}'"
                ) from e
        else:
            self._pressure = None

    def is_null(self):
        return self._sire_mols is None

    def _clear_state(self):
        self._omm_state = None
        self._omm_state_has_cv = (False, False)

    def _update_from(self, state, state_has_cv, nsteps_completed):
        if self.is_null():
            return

        self._current_step = nsteps_completed

        if not state_has_cv[0]:
            # there is no information to update
            return

        from ..legacy.Convert import (
            openmm_extract_coordinates_and_velocities,
            openmm_extract_coordinates,
            openmm_extract_space,
        )

        if self._sire_mols.num_atoms() == self._omm_mols.get_atom_index().count():
            # all of the atoms in all molecules are in the context,
            # and we can assume they are in atom index order
            mols_to_update = self._sire_mols.molecules()
        else:
            # some of the atoms aren't in the context, and they may be
            # in a different order
            mols_to_update = self._omm_mols.get_atom_index().atoms()
            mols_to_update.update(self._sire_mols.molecules())

        if state_has_cv[1]:
            # get velocities too
            mols_to_update = openmm_extract_coordinates_and_velocities(
                state,
                mols_to_update,
                # black auto-formats this to a long line
                perturbable_maps=self._omm_mols.get_lambda_lever().get_perturbable_molecule_maps(),  # noqa: E501
                map=self._map,
            )
        else:
            mols_to_update = openmm_extract_coordinates(
                state,
                mols_to_update,
                # black auto-formats this to a long line
                perturbable_maps=self._omm_mols.get_lambda_lever().get_perturbable_molecule_maps(),  # noqa: E501
                map=self._map,
            )

        self._sire_mols.update(mols_to_update.molecules())

        if self._ffinfo.space().is_periodic():
            # don't change the space if it is infinite - this
            # cannot change during the simulation and OpenMM
            # likes giving back fake spaces
            space = openmm_extract_space(state)
            self._sire_mols.set_property("space", space)
            self._ffinfo.set_space(space)

        self._sire_mols.set_property("time", self._current_time)

    def _enter_dynamics_block(self):
        if self._is_running:
            raise SystemError("Cannot start dynamics while it is already running!")

        self._is_running = True
        self._clear_state()

    def _exit_dynamics_block(
        self,
        save_frame: bool = False,
        save_energy: bool = False,
        lambda_windows=[],
        rest2_scale_factors=[],
        save_velocities: bool = False,
        delta_lambda: float = None,
        num_energy_neighbours: int = None,
        null_energy: float = None,
        excess_chemical_potential: float = None,
        num_waters: int = None,
    ):
        if not self._is_running:
            raise SystemError("Cannot stop dynamics that is not running!")

        import openmm
        from ..units import nanosecond, kcal_per_mol

        if save_frame:
            self._omm_state = self._omm_mols.getState(
                getEnergy=True,
                getPositions=True,
                getVelocities=save_velocities,
                enforcePeriodicBox=self._enforce_periodic_box,
            )

            self._omm_state_has_cv = (True, save_velocities)
        else:
            self._omm_state = self._omm_mols.getState(getEnergy=True)
            self._omm_state_has_cv = (False, False)

        current_time = (
            self._omm_state.getTime().value_in_unit(openmm.unit.nanosecond) * nanosecond
        )

        delta = current_time - self._elapsed_time

        self._elapsed_time = current_time
        self._current_time += delta

        # store the number of lambda windows
        if lambda_windows is not None:
            num_lambda_windows = len(lambda_windows)

            # compute energies for all windows
            if num_energy_neighbours is None:
                num_energy_neighbours = num_lambda_windows

        if save_energy:
            # should save energy here
            nrgs = {}

            nrgs["kinetic"] = (
                self._omm_state.getKineticEnergy().value_in_unit(
                    openmm.unit.kilocalorie_per_mole
                )
                * kcal_per_mol
            )

            nrgs["potential"] = (
                self._omm_state.getPotentialEnergy().value_in_unit(
                    openmm.unit.kilocalorie_per_mole
                )
                * kcal_per_mol
            )

            sim_lambda_value = self._omm_mols.get_lambda()
            sim_rest2_scale = self._omm_mols.get_rest2_scale()

            # Get the current volume.
            if self._pressure is not None:
                volume = self._omm_state.getPeriodicBoxVolume().value_in_unit(
                    openmm.unit.angstrom**3
                )

            # store the potential energy and accumulated non-equilibrium work
            if self._is_interpolate:
                nrg = nrgs["potential"]

                if (
                    len(self._lambda_interpolate) == 2
                    and sim_lambda_value != self._lambda_interpolate[0]
                ):
                    self._work += delta_lambda * (nrg - self._nrg_prev)
                self._nrg_prev = nrg
                nrgs["work"] = self._work
            else:
                nrg = (nrgs["potential"] / kcal_per_mol).value()
                if self._pressure is not None:
                    nrg += self._pressure * volume
                if excess_chemical_potential is not None:
                    nrg += excess_chemical_potential * num_waters
                nrgs[str(sim_lambda_value)] = nrg * kcal_per_mol

                if lambda_windows is not None:
                    # get the index of the simulation lambda value in the
                    # lambda windows list
                    try:
                        lambda_index = lambda_windows.index(sim_lambda_value)
                        has_lambda_index = True
                    except:
                        has_lambda_index = False

                    for i, (lambda_value, rest2_scale) in enumerate(
                        zip(lambda_windows, rest2_scale_factors)
                    ):
                        if lambda_value != sim_lambda_value:
                            if (
                                not has_lambda_index
                                or abs(lambda_index - i) <= num_energy_neighbours
                            ):
                                self._omm_mols.set_lambda(
                                    lambda_value,
                                    rest2_scale=rest2_scale,
                                    update_constraints=False,
                                )
                                nrg = self._omm_mols.get_potential_energy(
                                    to_sire_units=False
                                ).value_in_unit(openmm.unit.kilocalorie_per_mole)
                                if self._pressure is not None:
                                    nrg += self._pressure * volume
                                if excess_chemical_potential is not None:
                                    nrg += excess_chemical_potential * num_waters
                                nrgs[str(lambda_value)] = nrg * kcal_per_mol
                            else:
                                nrgs[str(lambda_value)] = null_energy * kcal_per_mol

                self._omm_mols.set_lambda(
                    sim_lambda_value,
                    rest2_scale=sim_rest2_scale,
                    update_constraints=False,
                )

            if self._is_interpolate:
                self._energy_trajectory.set(
                    self._current_time, nrgs, {"lambda": f"{sim_lambda_value:.8f}"}
                )
            else:
                self._energy_trajectory.set(
                    self._current_time, nrgs, {"lambda": str(sim_lambda_value)}
                )

            # update the interpolation lambda value
            if self._is_interpolate:
                if delta_lambda:
                    sim_lambda_value += delta_lambda
                    # clamp to [0, 1]
                    if sim_lambda_value < 0.0:
                        sim_lambda_value = 0.0
                    elif sim_lambda_value > 1.0:
                        sim_lambda_value = 1.0
                self._omm_mols.set_lambda(sim_lambda_value)

        self._is_running = False

        return (self._omm_state, self._omm_state_has_cv)

    def _get_current_state(
        self, include_coords: bool = False, include_velocities: bool = False
    ):
        if self._omm_state is not None:
            if self._omm_state_has_cv == (include_coords, include_velocities):
                return self._omm_state

        # we need to get this again as either it doesn't exist,
        # or we now want velocities as well
        if include_coords or include_velocities:
            self._omm_state = self._omm_mols.getState(
                getEnergy=True,
                getPositions=include_coords,
                getVelocities=include_velocities,
                enforcePeriodicBox=self._enforce_periodic_box,
            )
        else:
            self._omm_state = self._omm_mols.getState(getEnergy=True)

        self._omm_state_has_cv = (include_coords, include_velocities)

        return self._omm_state

    def timestep(self):
        if self.is_null():
            return 0
        else:
            from openmm.unit import picosecond as _omm_ps
            from ..units import picosecond as _sire_ps

            t = self._omm_mols.getIntegrator().getStepSize()

            return t.value_in_unit(_omm_ps) * _sire_ps

    def ensemble(self):
        if self.is_null():
            return None
        else:
            from ..move import Ensemble

            return Ensemble(map=self._map)

    def randomise_velocities(self, temperature=None, random_seed: int = None):
        """
        Set the velocities to random values, drawn from the Boltzmann
        distribution for the current temperature.

        Parameters
        ----------

        - temperature (temperature): The temperature to use. If None, then
          the current temperature will be used
        - random_seed (int): The random seed to use. If None, then
          a random seed will be generated
        """
        if self.is_null():
            return

        if temperature is not None:
            self.set_temperature(temperature, rescale_velocities=False)

        from ..units import kelvin

        if random_seed is None:
            self._omm_mols.setVelocitiesToTemperature(
                self.ensemble().temperature().to(kelvin)
            )
        else:
            self._omm_mols.setVelocitiesToTemperature(
                self.ensemble().temperature().to(kelvin), random_seed
            )

    def set_ensemble(self, ensemble, rescale_velocities: bool = True):
        """
        Set the ensemble for the dynamics. Note that this will
        only let you change the temperature and/or pressure of the
        ensemble. You can't change its fundemental nature.

        If rescale_velocities is True, then the velocities will
        be rescaled to the new temperature.
        """
        if self.is_null():
            return

        if ensemble.name() != self.ensemble().name():
            raise ValueError(
                "You cannot change the ensemble of the dynamics. "
                f"Currently the ensemble is {self.ensemble().name()}, "
                f"but you tried to set it to {ensemble.name()}."
            )

        if ensemble.temperature() != self.ensemble().temperature():
            self._map["temperature"] = ensemble.temperature()
            self._omm_mols.set_temperature(
                ensemble.temperature(), rescale_velocities=rescale_velocities
            )

        if ensemble.is_constant_pressure():
            if ensemble.pressure() != self.ensemble().pressure():
                self._map["pressure"] = ensemble.pressure()
                self._omm_mols.set_pressure(ensemble.pressure())

                from sire import u
                from sire.units import mole

                NA = 6.02214076e23 / mole
                pressure = u(self._map["pressure"].source())
                self._pressure = (pressure * NA).value()

    def set_temperature(self, temperature, rescale_velocities: bool = True):
        """
        Set the temperature for the dynamics. Note that this will only
        let you change the temperature of the ensemble.
        You can't change its fundemental nature.

        If rescale_velocities is True, then the velocities will be
        rescaled to the new temperature.
        """
        if self.is_null():
            return

        ensemble = self.ensemble()

        ensemble.set_temperature(temperature)

        self.set_ensemble(ensemble=ensemble, rescale_velocities=rescale_velocities)

    def set_pressure(self, pressure):
        """
        Set the pressure for the dynamics. Note that this will only
        let you change the pressure of the ensemble.
        You can't change its fundemental nature.
        """
        if self.is_null():
            return

        ensemble = self.ensemble()

        ensemble.set_pressure(pressure)

        self.set_ensemble(ensemble=ensemble)

    def constraint(self):
        if self.is_null():
            return None
        else:
            return self._omm_mols.get_constraint()

    def perturbable_constraint(self):
        if self.is_null():
            return None
        else:
            return self._omm_mols.get_perturbable_constraint()

    def get_constraints(self):
        """
        Return the actual list of constraints that have been applied
        to this system. This is two lists of atoms, plus a list of
        distances. The constraint is atom0[i]::atom1[i] with distance[i]
        """
        if self.is_null():
            return None
        else:
            return self._omm_mols.get_constraints()

    def get_schedule(self):
        if self.is_null():
            return None
        else:
            return self._omm_mols.get_lambda_schedule()

    def set_schedule(self, schedule):
        if not self.is_null():
            self._omm_mols.set_lambda_schedule(schedule)
            self._schedule_changed = True
            self.set_lambda(
                self._omm_mols.get_lambda(),
                rest2_scale=self._omm_mols.get_rest2_scale(),
            )

    def get_lambda(self):
        if self.is_null():
            return None
        else:
            return self._omm_mols.get_lambda()

    def set_lambda(
        self,
        lambda_value: float,
        rest2_scale: float = 1.0,
        update_constraints: bool = True,
    ):
        if not self.is_null():
            s = self.get_schedule()

            if s is None:
                return

            lambda_value = s.clamp(lambda_value)

            if (
                (not self._schedule_changed)
                and (lambda_value == self._omm_mols.get_lambda())
                and (rest2_scale == self._omm_mols.get_rest2_scale())
            ):
                # nothing to do
                return

            self._omm_mols.set_lambda(
                lambda_value=lambda_value,
                rest2_scale=rest2_scale,
                update_constraints=update_constraints,
            )
            self._schedule_changed = False
            self._clear_state()

    def integrator(self):
        if self.is_null():
            return None
        else:
            return self._omm_mols.getIntegrator()

    def info(self):
        if self.is_null():
            return None
        else:
            return self._ffinfo

    def platform(self):
        if self.is_null():
            return None
        else:
            return self._omm_mols.getPlatform().getName()

    def current_step(self):
        if self.is_null():
            return 0
        else:
            return self._current_step

    def current_time(self):
        if self.is_null():
            return 0
        else:
            return self._current_time

    def current_space(self):
        if self.is_null():
            return None
        else:
            return self._ffinfo.space()

    def current_walltime(self):
        if self.is_null():
            return 0
        else:
            return self._walltime

    def elapsed_time(self):
        if self.is_null():
            return 0
        else:
            return self._elapsed_time

    def walltime(self):
        if self.is_null():
            return 0
        else:
            return self._walltime

    def current_energy(self):
        if self.is_null():
            return 0
        else:
            from openmm.unit import kilocalorie_per_mole as _omm_kcal_mol
            from ..units import kcal_per_mol as _sire_kcal_mol

            state = self._get_current_state()

            nrg = state.getKineticEnergy() + state.getPotentialEnergy()

            return nrg.value_in_unit(_omm_kcal_mol) * _sire_kcal_mol

    def current_potential_energy(self):
        if self.is_null():
            return 0
        else:
            from openmm.unit import kilocalorie_per_mole as _omm_kcal_mol
            from ..units import kcal_per_mol as _sire_kcal_mol

            state = self._get_current_state()

            nrg = state.getPotentialEnergy()

            return nrg.value_in_unit(_omm_kcal_mol) * _sire_kcal_mol

    def current_kinetic_energy(self):
        if self.is_null():
            return 0
        else:
            from openmm.unit import kilocalorie_per_mole as _omm_kcal_mol
            from ..units import kcal_per_mol as _sire_kcal_mol

            state = self._get_current_state()

            nrg = state.getKineticEnergy()

            return nrg.value_in_unit(_omm_kcal_mol) * _sire_kcal_mol

    def energy_trajectory(self):
        return self._energy_trajectory.clone()

    def step(self, num_steps: int = 1):
        """
        Just perform 'num_steps' steps of dynamics, without saving
        anything or running anything in a background thread. This is
        designed for times when we want a minimial overhead, e.g.
        when we want to run a small number of steps quickly.
        """
        if self._is_running:
            raise SystemError("Cannot step dynamics while it is already running!")

        self._clear_state()

        self._omm_mols.getIntegrator().step(num_steps)

    def get_minimisation_log(self):
        """
        Return the log from the last minimisation
        """
        if self.is_null():
            return None
        else:
            try:
                return self._minimisation_log
            except Exception:
                return None

    def run_minimisation(
        self,
        max_iterations: int = 10000,
        tolerance: float = 10.0,
        max_restarts: int = 10,
        max_ratchets: int = 20,
        ratchet_frequency: int = 500,
        starting_k: float = 400.0,
        ratchet_scale: float = 10.0,
        max_constraint_error: float = 0.001,
        timeout: str = "300s",
    ):
        """
        Internal method that runs minimisation on the molecules.

        If the system is constrained, then a ratcheting algorithm is used.
        The constraints are replaced by harmonic restraints with an
        force constant based on `tolerance` and `starting_k`. Minimisation
        is performed, with the actual constrained bond lengths checked
        whenever minimisation converges, or when ratchet_frequency steps
        have completed (whichever is sooner). The force constant of
        the restraints is ratcheted up by `ratchet_scale`, and minimisation
        continues until there is no large change in energy or the maximum
        number of ratchets has been reached. In addition, at each ratchet,
        the actual bond lengths of constrained bonds are compared against
        the constrained values. If these have drifted too far away from
        the constrained values, then the minimisation is restarted,
        going back to the starting conformation and starting minimisation
        at one higher ratchet level. This will repeat a maximum of
        `max_restarts` times.

        If a stable structure cannot be reached, then an exception
        will be raised.

        Parameters:

        - max_iterations (int): The maximum number of iterations to run
        - tolerance (float): The tolerance to use for the minimisation
        - max_restarts (int): The maximum number of restarts before giving up
        - max_ratchets (int): The maximum number of ratchets before giving up
        - ratchet_frequency (int): The maximum number of steps between ratchets
        - starting_k (float): The starting value of k for the minimisation
        - ratchet_scale (float): The amount to scale k at each ratchet
        - max_constraint_error (float): The maximum error in the constraint in nm
        - timeout (float): The maximum time to run the minimisation for in seconds.
                           A value of <=0 will disable the timeout.
        """
        from ..legacy.Convert import minimise_openmm_context

        if max_iterations <= 0:
            max_iterations = 0

        try:
            from ..units import second
            from .. import u

            timeout = u(timeout)
            if not timeout.has_same_units(second):
                raise ValueError("'timeout' must have units of time")
        except:
            raise ValueError("Unable to parse 'timeout' as a time")

        self._clear_state()

        self._minimisation_log = minimise_openmm_context(
            self._omm_mols,
            tolerance=tolerance,
            max_iterations=max_iterations,
            max_restarts=max_restarts,
            max_ratchets=max_ratchets,
            ratchet_frequency=ratchet_frequency,
            starting_k=starting_k,
            ratchet_scale=ratchet_scale,
            max_constraint_error=max_constraint_error,
            timeout=timeout.to(second),
        )

    def _rebuild_and_minimise(self):
        if self.is_null():
            return

        from ..utils import Console

        Console.warning(
            "Something went wrong when running dynamics. The system will be "
            "minimised, and then dynamics will be attempted again. If an "
            "error still occurs, then it is likely that the step size is too "
            "large, the molecules are over-constrained, or there is something "
            "more fundemental going wrong..."
        )

        # rebuild the molecules
        from ..convert import to

        self._omm_mols = to(self._sire_mols, "openmm", map=self._map)

        # reset the water state
        if self._gcmc_sampler is not None:
            self._gcmc_sampler.push()
            self._gcmc_sampler._set_water_state(self._omm_mols)
            self._gcmc_sampler.pop()

        if self._save_crash_report:
            import openmm
            import numpy as np
            from copy import deepcopy
            from uuid import uuid4

            # Create a unique identifier for this crash report.
            crash_id = str(uuid4())[:8]

            # Get the current context and system.
            context = self._omm_mols
            system = deepcopy(context.getSystem())

            # Add each force to a unique group.
            for i, f in enumerate(system.getForces()):
                f.setForceGroup(i)

            # Create a new context.
            new_context = openmm.Context(system, deepcopy(context.getIntegrator()))
            new_context.setPositions(context.getState(getPositions=True).getPositions())

            # Write the  energies for each force group.
            with open(f"crash_{crash_id}.log", "w") as f:
                f.write(f"Current lambda: {str(self.get_lambda())}\n")
                for i, force in enumerate(system.getForces()):
                    state = new_context.getState(getEnergy=True, groups={i})
                    f.write(f"{force.getName()}, {state.getPotentialEnergy()}\n")

            # Save the serialised system.
            with open(f"system_{crash_id}.xml", "w") as f:
                f.write(openmm.XmlSerializer.serialize(system))

            # Save the positions.
            positions = (
                new_context.getState(getPositions=True).getPositions(asNumpy=True)
                / openmm.unit.nanometer
            )
            np.savetxt(f"positions_{crash_id}.txt", positions)

        self.run_minimisation()

    def run(
        self,
        time,
        save_frequency=None,
        frame_frequency=None,
        energy_frequency=None,
        lambda_windows=None,
        rest2_scale_factors=None,
        save_velocities: bool = None,
        save_frame_on_exit: bool = False,
        save_energy_on_exit: bool = False,
        auto_fix_minimise: bool = True,
        num_energy_neighbours: int = None,
        null_energy: str = None,
        excess_chemical_potential: float = None,
        num_waters: int = None,
        save_crash_report: bool = False,
    ):
        if self.is_null():
            return

        orig_args = {
            "time": time,
            "save_frequency": save_frequency,
            "frame_frequency": frame_frequency,
            "energy_frequency": energy_frequency,
            "lambda_windows": lambda_windows,
            "rest2_scale_factors": rest2_scale_factors,
            "save_velocities": save_velocities,
            "save_frame_on_exit": save_frame_on_exit,
            "save_energy_on_exit": save_energy_on_exit,
            "auto_fix_minimise": auto_fix_minimise,
            "num_energy_neighbours": num_energy_neighbours,
            "null_energy": null_energy,
            "excess_chemical_potential": excess_chemical_potential,
            "num_waters": num_waters,
            "save_crash_report": save_crash_report,
        }

        from concurrent.futures import ThreadPoolExecutor
        import openmm

        from ..units import picosecond
        from .. import u

        time = u(time)

        if save_frequency is not None:
            save_frequency = u(save_frequency)

        if frame_frequency is not None:
            frame_frequency = u(frame_frequency)

        if energy_frequency is not None:
            energy_frequency = u(energy_frequency)

        if null_energy is not None:
            null_energy = u(null_energy)
        else:
            null_energy = u("1e6 kcal/mol")

        if num_energy_neighbours is not None:
            try:
                num_energy_neighbours = int(num_energy_neighbours)
            except:
                num_energy_neighbours = len(lambda_windows)

        if excess_chemical_potential is not None:
            excess_chemical_potential = u(excess_chemical_potential).value()

        if num_waters is not None:
            try:
                num_waters = int(num_waters)
            except:
                raise ValueError("'num_waters' must be an integer")

        if save_crash_report is not None:
            if not isinstance(save_crash_report, bool):
                raise ValueError("'save_crash_report' must be True or False")
            self._save_crash_report = save_crash_report
        else:
            self._save_crash_report = False

        try:
            steps_to_run = int(time.to(picosecond) / self.timestep().to(picosecond))
        except Exception:
            # passed in the number of steps instead
            steps_to_run = int(time)

        if steps_to_run < 1:
            steps_to_run = 1

        if steps_to_run <= 0:
            return

        # get the energy and frame save frequencies
        if save_frequency != 0:
            if save_frequency is None:
                if self._map.specified("save_frequency"):
                    save_frequency = self._map["save_frequency"].value().to(picosecond)
                else:
                    save_frequency = 25
            else:
                save_frequency = save_frequency.to(picosecond)
            no_save = False
        else:
            no_save = True
            save_frequency = steps_to_run + 1

        if energy_frequency != 0:
            no_save_energy = False
            if energy_frequency is None:
                if self._map.specified("energy_frequency"):
                    energy_frequency = (
                        self._map["energy_frequency"].value().to(picosecond)
                    )
                else:
                    energy_frequency = save_frequency
                    no_save_energy = no_save
            else:
                energy_frequency = energy_frequency.to(picosecond)
        else:
            energy_frequency = steps_to_run + 1
            no_save_energy = True

        if frame_frequency != 0:
            no_save_frame = False
            if frame_frequency is None:
                if self._map.specified("frame_frequency"):
                    frame_frequency = (
                        self._map["frame_frequency"].value().to(picosecond)
                    )
                else:
                    frame_frequency = save_frequency
                    no_save_frame = no_save
            else:
                frame_frequency = frame_frequency.to(picosecond)
        else:
            frame_frequency = steps_to_run + 1
            no_save_frame = True

        completed = 0

        frame_frequency_steps = int(frame_frequency / self.timestep().to(picosecond))

        energy_frequency_steps = int(energy_frequency / self.timestep().to(picosecond))

        # If performing QM/MM lambda interpolation, then we just compute energies
        # for the pure MM (0.0) and QM (1.0) potentials.
        if self._is_interpolate:
            lambda_windows = [0.0, 1.0]
            # Work out the lambda increment.
            if isinstance(self._lambda_interpolate, list):
                divisor = (steps_to_run / energy_frequency_steps) - 1.0
                delta_lambda = (
                    self._lambda_interpolate[1] - self._lambda_interpolate[0]
                ) / divisor
            # Fixed lambda value.
            else:
                delta_lambda = None

            # Set the REST2 scaling factors.
            rest2_scale_factors = [1.0, 1.0]
        else:
            delta_lambda = None
            if lambda_windows is not None:
                if not isinstance(lambda_windows, list):
                    lambda_windows = [lambda_windows]
            else:
                if self._map.specified("lambda_windows"):
                    lambda_windows = self._map["lambda_windows"].value()

            if rest2_scale_factors is not None:
                if len(rest2_scale_factors) != len(lambda_windows):
                    raise ValueError(
                        "len(rest2_scale_factors) must be equal to len(lambda_windows)"
                    )
            else:
                if lambda_windows is not None:
                    rest2_scale_factors = [1.0] * len(lambda_windows)

        def runfunc(num_steps):
            try:
                integrator = self._omm_mols.getIntegrator()
                integrator.step(num_steps)
                return 0
            except Exception as e:
                return e

        def process_block(state, state_has_cv, nsteps_completed):
            self._update_from(state, state_has_cv, nsteps_completed)

            if state_has_cv[0] or state_has_cv[1]:
                self._sire_mols.save_frame(
                    map={
                        "space": self.current_space(),
                        "time": self.current_time(),
                    }
                )

        state = None
        state_has_cv = (False, False)
        saved_last_frame = False

        # whether the energy or frame were saved after the current block
        have_saved_frame = False
        have_saved_energy = False

        class NeedsMinimiseError(Exception):
            pass

        # store the previous time and number of steps to allow us to reset
        # when a crash occurs
        self._prev_current_time = self._current_time
        self._prev_elapsed_time = self._elapsed_time
        self._prev_step = self._current_step

        nsteps_before_run = self._current_step

        # if this is the first call, then set the save frequencies
        if nsteps_before_run == 0:
            self._next_save_frame = frame_frequency_steps
            self._next_save_energy = energy_frequency_steps
            self._prev_frame_frequency_steps = frame_frequency_steps
            self._prev_energy_frequency_steps = energy_frequency_steps
            self._prev_no_frame = no_save_frame
            self._prev_no_energy = no_save_energy
        # handle adjustments to the save frequencies
        else:
            if frame_frequency_steps != self._prev_frame_frequency_steps:
                if self._prev_no_frame:
                    self._next_save_frame = nsteps_before_run + frame_frequency_steps
                else:
                    self._next_save_frame = (
                        self._next_save_frame
                        + frame_frequency_steps
                        - self._prev_frame_frequency_steps
                    )
            if energy_frequency_steps != self._prev_energy_frequency_steps:
                if self._prev_no_energy:
                    self._next_save_energy = nsteps_before_run + energy_frequency_steps
                else:
                    self._next_save_energy = (
                        self._next_save_energy
                        + energy_frequency_steps
                        - self._prev_energy_frequency_steps
                    )
            self._prev_no_frame = no_save_frame
            self._prev_frame_frequency_steps = frame_frequency_steps
            self._prev_no_energy = no_save_energy
            self._prev_energy_frequency_steps = energy_frequency_steps

        from ..base import ProgressBar
        from ..units import second
        from datetime import datetime
        from math import isnan

        try:
            with ProgressBar(total=steps_to_run, text="dynamics") as progress:
                progress.set_speed_unit("steps / s")

                start_time = datetime.now()

                with ThreadPoolExecutor() as pool:
                    while completed < steps_to_run:
                        block_size = 50

                        steps_till_frame = self._next_save_frame - (
                            completed + nsteps_before_run
                        )
                        if steps_till_frame <= 0 or (
                            steps_till_frame <= block_size
                            and steps_till_frame <= steps_to_run - completed
                        ):
                            save_frame = True
                            self._next_save_frame += frame_frequency_steps
                            if frame_frequency_steps < block_size:
                                block_size = frame_frequency_steps
                        else:
                            save_frame = False

                        steps_till_energy = self._next_save_energy - (
                            completed + nsteps_before_run
                        )
                        if steps_till_energy <= 0 or (
                            steps_till_energy <= block_size
                            and steps_till_energy <= steps_to_run - completed
                        ):
                            save_energy = True
                            self._next_save_energy += energy_frequency_steps
                            if energy_frequency_steps < block_size:
                                block_size = energy_frequency_steps
                        else:
                            save_energy = False

                        # save the last frame if we're about to exit and the user
                        # has requested it
                        if (
                            save_frame_on_exit
                            and completed + block_size >= steps_to_run
                        ):
                            save_frame = True

                        # save the last energy if we're about to exit and the user
                        # has requested it
                        if (
                            save_energy_on_exit
                            and completed + block_size >= steps_to_run
                        ):
                            save_energy = True

                        self._enter_dynamics_block()

                        # process the last block in the foreground
                        if state is not None:
                            process_promise = pool.submit(
                                process_block,
                                state,
                                state_has_cv,
                                nsteps_before_run + completed,
                            )
                        else:
                            process_promise = None

                        nrun = block_size

                        # this block will exceed the run time so reduce the size
                        if nrun > steps_to_run - completed:
                            nrun = steps_to_run - completed

                        # run the current block in the background
                        run_promise = pool.submit(runfunc, nrun)

                        result = None

                        while not run_promise.done():
                            try:
                                result = run_promise.result(timeout=1.0)
                            except Exception as e:
                                if (
                                    "NaN" in str(e)
                                    and not have_saved_frame
                                    and not have_saved_energy
                                    and auto_fix_minimise
                                ):
                                    raise NeedsMinimiseError()

                        # catch rare edge case where the promise timed
                        # out, but then completed before the .done()
                        # test in the next loop iteration
                        if result is None:
                            result = run_promise.result()

                        if result == 0:
                            completed += nrun
                            progress.set_progress(completed)
                            run_promise = None
                        else:
                            # make sure we finish processing the last block
                            if process_promise is not None:
                                try:
                                    process_promise.result()
                                except Exception as e:
                                    if (
                                        "NaN" in str(e)
                                        and not have_saved_frame
                                        and not have_saved_energy
                                        and auto_fix_minimise
                                    ):
                                        raise NeedsMinimiseError()

                            if (
                                not have_saved_frame
                                and not have_saved_energy
                                and auto_fix_minimise
                            ):
                                raise NeedsMinimiseError()

                            # something went wrong - re-raise the exception
                            raise result

                        # make sure we've finished processing the last block
                        if process_promise is not None:
                            process_promise.result()
                            process_promise = None
                            saved_last_frame = True

                        # get the state, including coordinates and velocities
                        state, state_has_cv = self._exit_dynamics_block(
                            save_frame=save_frame,
                            save_energy=save_energy,
                            lambda_windows=lambda_windows,
                            rest2_scale_factors=rest2_scale_factors,
                            save_velocities=save_velocities,
                            delta_lambda=delta_lambda,
                            num_energy_neighbours=num_energy_neighbours,
                            null_energy=null_energy.value(),
                            excess_chemical_potential=excess_chemical_potential,
                            num_waters=num_waters,
                        )

                        saved_last_frame = False

                        have_saved_frame = save_frame
                        have_saved_energy = save_energy

                        kinetic_energy = state.getKineticEnergy().value_in_unit(
                            openmm.unit.kilocalorie_per_mole
                        )

                        ke_per_atom = kinetic_energy / self._num_atoms

                        if isnan(ke_per_atom) or ke_per_atom > 1000:
                            # The system has blown up!
                            state = None
                            saved_last_frame = True

                            if (
                                not have_saved_frame
                                and not have_saved_energy
                                and auto_fix_minimise
                            ):
                                raise NeedsMinimiseError()

                            raise RuntimeError(
                                "The kinetic energy has exceeded 1000 kcal "
                                f"mol-1 per atom (it is {ke_per_atom} kcal "
                                f"mol-1 atom-1, and {kinetic_energy} kcal "
                                "mol-1 total). This suggests that the "
                                "simulation has become unstable. Try reducing "
                                "the timestep and/or minimising the system "
                                "and run again."
                            )

                self._walltime += (datetime.now() - start_time).total_seconds() * second

            if state is not None and not saved_last_frame:
                # we can process the last block in the main thread
                process_block(
                    state=state,
                    state_has_cv=state_has_cv,
                    nsteps_completed=nsteps_before_run + completed,
                )

        except NeedsMinimiseError:
            from openmm.unit import picosecond

            # try to fix this problem by minimising,
            # then running again
            self._is_running = False
            self._clear_state()
            self._rebuild_and_minimise()
            orig_args["auto_fix_minimise"] = False
            self._current_step = self._prev_step
            self._current_time = self._prev_current_time
            self._elapsed_time = self._prev_elapsed_time
            self._omm_mols.setTime(
                self._prev_elapsed_time.to("picosecond") * picosecond
            )

            self.run(**orig_args)
            return

    def to_xml(self, f=None):
        """
        Save the current state of the dynamics to XML.
        This is mostly used for debugging. This will return the
        XML string if 'f' is None. Otherwise it will write the
        XML to 'f' (either a filename, or a FILE object)
        """
        return self._omm_mols.to_xml(f=f)

    def commit(self, return_as_system: bool = False):
        if self.is_null():
            return

        self._update_from(
            state=self._get_current_state(include_coords=True, include_velocities=True),
            state_has_cv=(True, True),
            nsteps_completed=self._current_step,
        )

        self._sire_mols.set_energy_trajectory(self._energy_trajectory, map=self._map)

        self._sire_mols.set_ensemble(self.ensemble())

        if return_as_system:
            return self._sire_mols.clone()

        elif self._orig_mols is not None:
            from ..system import System

            if System.is_system(self._orig_mols):
                return self._sire_mols.clone()
            else:
                r = self._orig_mols.clone()
                r.update(self._sire_mols.molecules())
                return r
        else:
            return self._sire_mols.clone()


def _add_extra(extras, key, value):
    if value is not None:
        extras[key] = value


class Dynamics:
    """
    Class that runs dynamics on the contained molecule(s). Note that
    this class is not designed to be constructed directly. You should only
    use this class by calling `.dynamics()` on the molecules(s)
    you want to simulate
    """

    def __init__(
        self,
        mols=None,
        map=None,
        cutoff=None,
        cutoff_type=None,
        timestep=None,
        constraint=None,
        perturbable_constraint=None,
        schedule=None,
        lambda_value=None,
        swap_end_states=None,
        ignore_perturbations=None,
        shift_delta=None,
        shift_coulomb=None,
        coulomb_power=None,
        restraints=None,
        fixed=None,
        qm_engine=None,
        lambda_interpolate=None,
    ):
        from ..base import create_map
        from .. import u

        extras = {}

        _add_extra(extras, "cutoff", cutoff)
        _add_extra(extras, "cutoff_type", cutoff_type)

        if timestep is not None:
            _add_extra(extras, "timestep", u(timestep))

        _add_extra(extras, "constraint", constraint)
        _add_extra(extras, "perturbable_constraint", perturbable_constraint)
        _add_extra(extras, "schedule", schedule)
        _add_extra(extras, "lambda", lambda_value)
        _add_extra(extras, "swap_end_states", swap_end_states)
        _add_extra(extras, "ignore_perturbations", ignore_perturbations)

        if shift_delta is not None:
            _add_extra(extras, "shift_delta", u(shift_delta))

        if shift_coulomb is not None:
            _add_extra(extras, "shift_coulomb", u(shift_coulomb))

        _add_extra(extras, "coulomb_power", coulomb_power)
        _add_extra(extras, "restraints", restraints)
        _add_extra(extras, "fixed", fixed)
        _add_extra(extras, "qm_engine", qm_engine)
        _add_extra(extras, "lambda_interpolate", lambda_interpolate)

        map = create_map(map, extras)

        self._d = DynamicsData(mols=mols, map=map)

        # Save the original view too, as this is the view that
        # will be returned to the user
        self._view = mols

    def __str__(self):
        speed = self.time_speed()

        if speed == 0:
            if self.current_step() > 0:
                return (
                    f"Dynamics(completed={self.current_time()}, "
                    f"energy={self.current_energy()}, speed=FAST ns day-1)"
                )
            else:
                return "Dynamics(completed=0)"
        else:
            return (
                f"Dynamics(completed={self.current_time()}, "
                f"energy={self.current_energy()}, speed={speed:.1f} ns day-1)"
            )

    def __repr__(self):
        return self.__str__()

    def minimise(
        self,
        max_iterations: int = 10000,
        tolerance: float = 10.0,
        max_restarts: int = 10,
        max_ratchets: int = 20,
        ratchet_frequency: int = 500,
        starting_k: float = 400.0,
        ratchet_scale: float = 10.0,
        max_constraint_error: float = 0.001,
        timeout: str = "300s",
    ):
        """
        Internal method that runs minimisation on the molecules.

        If the system is constrained, then a ratcheting algorithm is used.
        The constraints are replaced by harmonic restraints with an
        force constant based on `tolerance` and `starting_k`. Minimisation
        is performed, with the actual constrained bond lengths checked
        whenever minimisation converges, or when ratchet_frequency steps
        have completed (whichever is sooner). The force constant of
        the restraints is ratcheted up by `ratchet_scale`, and minimisation
        continues until there is no large change in energy or the maximum
        number of ratchets has been reached. In addition, at each ratchet,
        the actual bond lengths of constrained bonds are compared against
        the constrained values. If these have drifted too far away from
        the constrained values, then the minimisation is restarted,
        going back to the starting conformation and starting minimisation
        at one higher ratchet level. This will repeat a maximum of
        `max_restarts` times.

        If a stable structure cannot be reached, then an exception
        will be raised.

        Parameters:

        - max_iterations (int): The maximum number of iterations to run
        - tolerance (float): The tolerance to use for the minimisation
        - max_restarts (int): The maximum number of restarts before giving up
        - max_ratchets (int): The maximum number of ratchets before giving up
        - ratchet_frequency (int): The maximum number of steps between ratchets
        - starting_k (float): The starting value of k for the minimisation
        - ratchet_scale (float): The amount to scale k at each ratchet
        - max_constraint_error (float): The maximum error in the constraints in nm
        - timeout (float): The maximum time to run the minimisation for in seconds.
                           A value of <=0 will disable the timeout.
        """
        if not self._d.is_null():
            self._d.run_minimisation(
                max_iterations=max_iterations,
                tolerance=tolerance,
                max_restarts=max_restarts,
                max_ratchets=max_ratchets,
                ratchet_frequency=ratchet_frequency,
                starting_k=starting_k,
                ratchet_scale=ratchet_scale,
                max_constraint_error=max_constraint_error,
                timeout=timeout,
            )

        return self

    def step(self, num_steps: int = 1):
        """
        Simple function that performs `num_steps` steps of dynamics.
        This does not save any frames or energies - it is designed for
        times when you want to run a small number of steps quickly
        with minimal overhead.
        """
        if not self._d.is_null():
            self._d.step(num_steps=num_steps)

        return self

    def run(
        self,
        time,
        save_frequency=None,
        frame_frequency=None,
        energy_frequency=None,
        lambda_windows=None,
        rest2_scale_factors=None,
        save_velocities: bool = None,
        save_frame_on_exit: bool = False,
        save_energy_on_exit: bool = False,
        auto_fix_minimise: bool = True,
        num_energy_neighbours: int = None,
        null_energy: str = None,
        excess_chemical_potential: str = None,
        num_waters: int = None,
        save_crash_report: bool = False,
    ):
        """
        Perform dynamics on the molecules.

        Parameters

        time: Time
            The amount of simulation time to run, e.g.
            dynamics.run(sr.u("5 ps")) would perform
            5 picoseconds of dynamics. The number of steps
            is determined automatically based on the current
            timestep (e.g. if the timestep was 1 femtosecond,
            then 5 picoseconds would involve running 5000 steps)

        save_frequency: Time
            The amount of simulation time between saving frames
            (coordinates, velocities) and energies from the trajectory. The
            number of timesteps between saves will depend on the
            timestep. For example, if save_frequency was 0.1 picoseconds
            and the timestep was 2 femtoseconds, then the coordinates
            would be saved every 50 steps of dynamics. Note that
            `frame_frequency` or `energy_frequency` can be used
            to override the frequency of saving frames or energies,
            if you want them to be saved with different frequencies.
            Specifying both will mean that the value of
            `save_frequency` will be ignored.

        frame_frequency: Time
            The amount of simulation time between saving
            frames (coordinates, velocities) from the trajectory.
            The number of timesteps between saves will depend on the
            timestep. For example, if save_frequency was 0.1 picoseconds
            and the timestep was 2 femtoseconds, then the coordinates
            would be saved every 50 steps of dynamics. The energies
            will be saved into this object and are accessible via the
            `energy_trajectory` function.

        energy_frequency: Time
            The amount of simulation time between saving
            energies (kinetic and potential) from the trajectory.
            The number of timesteps between saves will depend on the
            timestep. For example, if save_frequency was 0.1 picoseconds
            and the timestep was 2 femtoseconds, then the coordinates
            would be saved every 50 steps of dynamics. The energies
            will be saved into this object and are accessible via the
            `energy_trajectory` function.

        lambda_windows: list[float]
            The values of lambda for which the potential energy will be
            evaluated at every save. If this is None (the default) then
            only the current energy will be saved every `energy_frequency`
            time. If this is not None, then the potential energy for
            each of the lambda values in this list will be saved. Note that
            we always save the potential energy of the simulated lambda
            value, even if it is not in the list of lambda windows.

        rest2_scale_factors: list[float]
            The scaling factors for the REST2 region for each lambda
            window.

        save_velocities: bool
            Whether or not to save the velocities when running dynamics.
            By default this is False. Set this to True if you are
            interested in saving the velocities.

        save_frame_on_exit: bool
            Whether to save a trajectory frame on exit, regardless of
            whether the frame frequency has been reached.

        save_energy_on_exit: bool
            Whether to save the energy on exit, regardless of whether
            the energy frequency has been reached.

        save_crash_report: bool
            Whether to save a crash report if the dynamics fails due to an
            instability. This will save a named log file containing the energy
            for each force, an XML file containing the OpenMM system at the
            start of the dynamics block, and a NumPy text file containing
            the atomic positions at the start of the dynamics block. This
            option is only used when auto_fix_minimise is True.

        auto_fix_minimise: bool
            Whether or not to automatically run minimisation if the
            trajectory exits with an error in the first few steps.
            Such failures often indicate that the system needs
            minimsing. This automatically runs the minimisation
            in these cases, and then runs the requested dynamics.

        num_energy_neighbours: int
            The number of neighbouring windows to use when computing
            the energy trajectory for the simulation lambda value.
            This can be used to compute energies over a subset of the
            values in 'lambda_windows', hence reducing the cost of
            computing the energy trajectory. Note that the simulation
            lambda value must be contained in 'lambda_windows', so it
            is recommended that the values are rounded. A value of
            'null_energy' will be added to the energy trajectory for the
            lambda windows that are omitted. Note that a similar result
            can be achieved by simply removing any lambda values from
            'lambda_windows' that you don't want, but this will result
            in an energy trajectory that only contains results for
            the specified lambda values. If None, then all lambda windows
            will be used.

        null_energy: str
            The energy value to use for lambda windows that are not
            being computed as part of the energy trajectory, i.e. when
            'num_energy_neighbours' is less than len(lambda_windows).
            By default, a value of '10000 kcal mol-1' is used.

        excess_chemical_potential: str
            The excess chemical potential of water. This is required when
            running dynamics as part of a Grand Canonical water sampling
            simulation.

        num_waters: int
            The current number of water molecules in the simulation box.
            This is used to compute the Grand Canonical contribution to
            the potential energy.
        """
        if not self._d.is_null():
            if save_velocities is None:
                if self._d._map.specified("save_velocities"):
                    save_velocities = self._d._map["save_velocities"].value().as_bool()
                else:
                    save_velocities = False

            if save_crash_report is None:
                if self._d._map.specified("save_crash_report"):
                    save_crash_report = (
                        self._d._map["save_crash_report"].value().as_bool()
                    )
                else:
                    save_crash_report = False

            self._d.run(
                time=time,
                save_frequency=save_frequency,
                frame_frequency=frame_frequency,
                energy_frequency=energy_frequency,
                lambda_windows=lambda_windows,
                rest2_scale_factors=rest2_scale_factors,
                save_velocities=save_velocities,
                save_frame_on_exit=save_frame_on_exit,
                save_energy_on_exit=save_energy_on_exit,
                save_crash_report=save_crash_report,
                auto_fix_minimise=auto_fix_minimise,
                num_energy_neighbours=num_energy_neighbours,
                null_energy=null_energy,
                excess_chemical_potential=excess_chemical_potential,
                num_waters=num_waters,
            )

        return self

    def get_schedule(self):
        """
        Return the LambdaSchedule that shows how lambda changes the
        underlying forcefield parameters in the system.
        Returns None if this isn't a perturbable system.
        """
        return self._d.get_schedule()

    def set_schedule(self, schedule):
        """
        Set the LambdaSchedule that will be used to control how
        lambda changes the underlying forcefield parameters
        in the system. This does nothing if this isn't
        a perturbable system
        """
        self._d.set_schedule(schedule)

    def get_lambda(self):
        """
        Return the current value of lambda for this system. This
        does nothing if this isn't a perturbable system
        """
        return self._d.get_lambda()

    def set_lambda(
        self,
        lambda_value: float,
        rest2_scale: float = 1.0,
        update_constraints: bool = True,
    ):
        """
        Set the current value of lambda for this system. This will
        update the forcefield parameters in the context according
        to the data in the LambdaSchedule. This does nothing if
        this isn't a perturbable system.

        The `rest2_scale` parameter specifies the temperature of the
        REST2 region relative to the rest of the system at the specified
        lambda value.

        If `update_constraints` is True, then this will also update
        the constraint length of any constrained perturbable bonds.
        These will be set to the r0 value for that bond at this
        value of lambda. If `update_constraints` is False, then
        the constraint will not be changed.
        """
        self._d.set_lambda(
            lambda_value=lambda_value,
            rest2_scale=rest2_scale,
            update_constraints=update_constraints,
        )

    def get_rest2_scale(self):
        """
        Return the current REST2 scaling factor.
        """
        if self.is_null():
            return None
        return self._d.get_rest2_scale()

    def set_rest2_scale(self, rest2_scale: float):
        """
        Set the current REST2 scaling factor.
        """
        self._d.set_rest2_scale(rest2_scale=rest2_scale)

    def ensemble(self):
        """
        Return the ensemble in which the simulation is being performed
        """
        return self._d.ensemble()

    def set_ensemble(self, ensemble):
        """
        Set the ensemble that should be used for this simulation. Note
        that you can only use this function to change temperature and/or
        pressure values. You can't change the fundemental ensemble
        of the simulation.
        """
        self._d.set_ensemble(ensemble)

    def set_temperature(self, temperature, rescale_velocities: bool = True):
        """
        Set the temperature for the dynamics. Note that this will only
        let you change the temperature of the ensemble.
        You can't change its fundemental nature.

        If rescale_velocities is True, then the velocities will be
        rescaled to the new temperature.
        """
        self._d.set_temperature(
            temperature=temperature, rescale_velocities=rescale_velocities
        )

    def set_pressure(self, pressure):
        """
        Set the pressure for the dynamics. Note that this will only
        let you change the pressure of the ensemble.
        You can't change its fundemental nature.
        """
        self._d.set_pressure(pressure=pressure)

    def randomise_velocities(self, temperature=None, random_seed: int = None):
        """
        Set the velocities to random values, drawn from the Boltzmann
        distribution for the current temperature.

        Parameters
        ----------

        - temperature (temperature): The temperature to use. If None, then
          the current temperature will be used
        - random_seed (int): The random seed to use. If None, then
          a random seed will be generated
        """
        self._d.randomise_velocities(temperature=temperature, random_seed=random_seed)

    def constraint(self):
        """
        Return the constraint used for the dynamics (e.g. constraining
        bonds involving hydrogens etc.)
        """
        return self._d.constraint()

    def perturbable_constraint(self):
        """
        Return the perturbable constraint used for the dynamics (e.g.
        constraining bonds involving hydrogens etc.)
        """
        return self._d.perturbable_constraint()

    def get_constraints(self):
        """
        Return the actual list of constraints that have been applied
        to this system. This is two lists of atoms, plus a list of
        distances. The constraint is atom0[i]::atom1[i] with distance[i]
        """
        return self._d.get_constraints()

    def integrator(self):
        """
        Return the integrator that is used to run dynamics
        """
        return self._d.integrator()

    def context(self):
        """
        Return the underlying OpenMM context that is being driven by this
        dynamics object.
        """
        return self._d._omm_mols

    def info(self):
        """
        Return the information that describes the forcefield that will
        be used for dynamics
        """
        return self._d.info()

    def platform(self):
        """
        Return the name of the OpenMM platform being used for dynamics
        """
        return self._d.platform()

    def timestep(self):
        """
        Return the timestep used for this simulation
        """
        return self._d.timestep()

    def current_step(self):
        """
        Return the current number of completed steps of dynamics
        """
        return self._d.current_step()

    def current_time(self):
        """
        Return the current amount of completed time of dynamics
        """
        return self._d.current_time()

    def elapsed_time(self):
        """
        Return the total amount of elapsed time of dynamics. This
        will be the same as the current_time if this run started
        from time=0. Otherwise this will be the difference between
        the start time and the current time.
        """
        return self._d.elapsed_time()

    def walltime(self):
        """
        Return the walltime (real actual runtime) of dynamics, i.e.
        how much real time has been consumed by the simulation
        """
        return self._d.walltime()

    def step_speed(self, time_unit=None):
        """
        Return the speed of this simulation in number of steps per second
        """
        if time_unit is None:
            from ..units import second

            time_unit = second

        nsteps = self.current_step()
        time = self.walltime().to(time_unit)

        if time < 0.000001:
            return 0
        else:
            return nsteps / time

    def time_speed(self, elapsed_unit=None, time_unit=None):
        """
        Return the speed of this simulation in simulation time per
        real time (e.g. nanoseconds simulated per day)
        """
        if elapsed_unit is None:
            from ..units import nanosecond

            elapsed_unit = nanosecond

        if time_unit is None:
            from ..units import day

            time_unit = day

        elapsed = self.elapsed_time().to(elapsed_unit)
        time = self.walltime().to(time_unit)

        if time < 0.000001:
            return 0
        else:
            return elapsed / time

    def current_space(self):
        """
        Return the current space in which the simulation is being performed
        """
        return self._d.current_space()

    def current_energy(self):
        """
        Return the current total energy
        """
        return self._d.current_energy()

    def current_potential_energy(self, lambda_values=None, rest2_scale_factors=None):
        """
        Return the current potential energy.

        If `lambda_values` is passed (which should be a list of
        lambda values) then this will return the energies
        (as a list) at the requested lambda values

        If `rest2_scale_factors` is passed, then these will be
        used to scale the temperature of the REST2 region at each
        lambda value.
        """
        if lambda_values is None:
            return self._d.current_potential_energy()
        else:
            if not isinstance(lambda_values, list):
                lambda_values = [lambda_values]
            if rest2_scale_factors is None:
                rest2_scale_factors = [1.0] * len(lambda_values)
            else:
                if not isinstance(rest2_scale_factors, list):
                    rest2_scale_factors = [rest2_scale_factors]
                else:
                    if len(rest2_scale_factors) != len(lambda_values):
                        raise ValueError(
                            "len(rest2_scale_factors) must be equal to len(lambda_values)"
                        )

            # save the current value of lambda so we
            # can restore it
            old_lambda = self.get_lambda()
            old_rest2_scale = self.get_rest2_scale()

            nrgs = []

            try:
                for lambda_value, rest2_scale in zip(
                    lambda_values, rest2_scale_factors
                ):
                    self.set_lambda(
                        lambda_value, rest2_scale=rest2_scale, update_constraints=False
                    )
                    nrgs.append(self._d.current_potential_energy())
            except Exception:
                self.set_lambda(old_lambda, old_rest2_scale, update_constraints=False)
                raise

            self.set_lambda(
                old_lambda, rest2_scale=old_rest2_scale, update_constraints=False
            )

            return nrgs

    def current_kinetic_energy(self):
        """
        Return the current kinetic energy
        """
        return self._d.current_kinetic_energy()

    def energy_trajectory(self, to_pandas: bool = False, to_alchemlyb: bool = False):
        """
        Return the energy trajectory. This is the trajectory of
        energy values that have been captured during dynamics.

        If 'to_pandas' is True, (the default) then this will
        be returned as a pandas dataframe, with times and energies
        in the defined default units
        """
        t = self._d.energy_trajectory()

        if to_pandas or to_alchemlyb:
            return t.to_pandas(
                to_alchemlyb=to_alchemlyb,
            )
        else:
            return t

    def to_xml(self, f=None):
        """
        Save the current state of the dynamics to XML.
        This is mostly used for debugging. This will return the
        XML string if 'f' is None. Otherwise it will write the
        XML to 'f' (either a filename, or a FILE object)
        """
        return self._d.to_xml(f=f)

    def commit(self, return_as_system: bool = False):
        """
        Commit the dynamics and return the molecules after the simulation.
        Normally this will return the same view of as was used for
        construction. If `return_as_system` is True, then this will
        return a System object instead.
        """
        if not self._d.is_null():
            return self._d.commit(return_as_system=return_as_system)
        else:
            return None

    def __call__(
        self,
        time,
        save_frequency=None,
        frame_frequency=None,
        energy_frequency=None,
        lambda_windows=None,
        save_velocities: bool = True,
        auto_fix_minimise: bool = True,
    ):
        """
        Perform dynamics on the molecules.

        Parameters

        time: Time
            The amount of simulation time to run, e.g.
            dynamics.run(sr.u("5 ps")) would perform
            5 picoseconds of dynamics. The number of steps
            is determined automatically based on the current
            timestep (e.g. if the timestep was 1 femtosecond,
            then 5 picoseconds would involve running 5000 steps)

        save_frequency: Time
            The amount of simulation time between saving frames
            (coordinates, velocities) and energies from the trajectory. The
            number of timesteps between saves will depend on the
            timestep. For example, if save_frequency was 0.1 picoseconds
            and the timestep was 2 femtoseconds, then the coordinates
            would be saved every 50 steps of dynamics. Note that
            `frame_frequency` or `energy_frequency` can be used
            to override the frequency of saving frames or energies,
            if you want them to be saved with different frequencies.
            Specifying both will mean that the value of
            `save_frequency` will be ignored.

        frame_frequency: Time
            The amount of simulation time between saving
            frames (coordinates, velocities) from the trajectory.
            The number of timesteps between saves will depend on the
            timestep. For example, if save_frequency was 0.1 picoseconds
            and the timestep was 2 femtoseconds, then the coordinates
            would be saved every 50 steps of dynamics. The energies
            will be saved into this object and are accessible via the
            `energy_trajectory` function.

        energy_frequency: Time
            The amount of simulation time between saving
            energies (kinetic and potential) from the trajectory.
            The number of timesteps between saves will depend on the
            timestep. For example, if save_frequency was 0.1 picoseconds
            and the timestep was 2 femtoseconds, then the coordinates
            would be saved every 50 steps of dynamics. The energies
            will be saved into this object and are accessible via the
            `energy_trajectory` function.

        lambda_windows: list[float]
            The values of lambda for which the potential energy will be
            evaluated at every save. If this is None (the default) then
            only the current energy will be saved every `energy_frequency`
            time. If this is not None, then the potential energy for
            each of the lambda values in this list will be saved. Note that
            we always save the potential energy of the simulated lambda
            value, even if it is not in the list of lambda windows.

        save_velocities: bool
            Whether or not to save the velocities when running dynamics.
            By default this is True. Set this to False if you aren't
            interested in saving the velocities.

        auto_fix_minimise: bool
            Whether or not to automatically run minimisation if the
            trajectory exits with an error in the first few steps.
            Such failures often indicate that the system needs
            minimsing. This automatically runs the minimisation
            in these cases, and then runs the requested dynamics.
        """
        return self.run(
            time=time,
            save_frequency=save_frequency,
            frame_frequency=frame_frequency,
            energy_frequency=energy_frequency,
            lambda_windows=lambda_windows,
            save_velocities=save_velocities,
            auto_fix_minimise=auto_fix_minimise,
        ).commit()
