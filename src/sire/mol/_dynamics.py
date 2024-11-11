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
            self._elapsed_time = 0 * nanosecond
            self._walltime = 0 * nanosecond
            self._is_running = False
            self._schedule_changed = False

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

        else:
            self._sire_mols = None
            self._energy_trajectory = None

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
        save_velocities: bool = False,
        delta_lambda: float = None,
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

            # Store the potential energy and accumulated non-equilibrium work.
            if self._is_interpolate:
                nrg = nrgs["potential"]

                if sim_lambda_value != 0.0:
                    self._work += delta_lambda * (nrg - self._nrg_prev)
                self._nrg_prev = nrg
                nrgs["work"] = self._work
            else:
                nrgs[str(sim_lambda_value)] = nrgs["potential"]

                if lambda_windows is not None:
                    for lambda_value in lambda_windows:
                        if lambda_value != sim_lambda_value:
                            self._omm_mols.set_lambda(
                                lambda_value, update_constraints=False
                            )
                            nrgs[str(lambda_value)] = (
                                self._omm_mols.get_potential_energy(
                                    to_sire_units=False
                                ).value_in_unit(openmm.unit.kilocalorie_per_mole)
                                * kcal_per_mol
                            )

                self._omm_mols.set_lambda(sim_lambda_value, update_constraints=False)

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
            self.set_lambda(self._omm_mols.get_lambda())

    def get_lambda(self):
        if self.is_null():
            return None
        else:
            return self._omm_mols.get_lambda()

    def set_lambda(self, lambda_value: float, update_constraints: bool = True):
        if not self.is_null():
            s = self.get_schedule()

            if s is None:
                return

            lambda_value = s.clamp(lambda_value)

            if (not self._schedule_changed) and (
                lambda_value == self._omm_mols.get_lambda()
            ):
                # nothing to do
                return

            self._omm_mols.set_lambda(
                lambda_value=lambda_value, update_constraints=update_constraints
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

        self.run_minimisation()

    def run(
        self,
        time,
        save_frequency=None,
        frame_frequency=None,
        energy_frequency=None,
        lambda_windows=None,
        save_velocities: bool = None,
        auto_fix_minimise: bool = True,
    ):
        if self.is_null():
            return

        orig_args = {
            "time": time,
            "save_frequency": save_frequency,
            "frame_frequency": frame_frequency,
            "energy_frequency": energy_frequency,
            "lambda_windows": lambda_windows,
            "save_velocities": save_velocities,
            "auto_fix_minimise": auto_fix_minimise,
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

        if energy_frequency != 0:
            if energy_frequency is None:
                if self._map.specified("energy_frequency"):
                    energy_frequency = (
                        self._map["energy_frequency"].value().to(picosecond)
                    )
                else:
                    energy_frequency = save_frequency
            else:
                energy_frequency = energy_frequency.to(picosecond)

        if frame_frequency != 0:
            if frame_frequency is None:
                if self._map.specified("frame_frequency"):
                    frame_frequency = (
                        self._map["frame_frequency"].value().to(picosecond)
                    )
                else:
                    frame_frequency = save_frequency
            else:
                frame_frequency = frame_frequency.to(picosecond)

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
        else:
            delta_lambda = None
            if lambda_windows is not None:
                if not isinstance(lambda_windows, list):
                    lambda_windows = [lambda_windows]
            else:
                if self._map.specified("lambda_windows"):
                    lambda_windows = self._map["lambda_windows"].value()

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

        block_size = 50

        state = None
        state_has_cv = (False, False)
        saved_last_frame = False

        # whether the energy or frame were saved after the current block
        have_saved_frame = False
        have_saved_energy = False

        class NeedsMinimiseError(Exception):
            pass

        nsteps_before_run = self._current_step
        if nsteps_before_run == 0:
            self._next_save_frame = frame_frequency_steps
            self._next_save_energy = energy_frequency_steps

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
                        steps_till_frame = self._next_save_frame - (
                            completed + nsteps_before_run
                        )
                        if steps_till_frame <= 0 or steps_till_frame == steps_to_run:
                            save_frame = True
                            self._next_save_frame += frame_frequency_steps
                        else:
                            save_frame = False

                        steps_till_energy = self._next_save_energy - (
                            completed + nsteps_before_run
                        )
                        if steps_till_energy <= 0 or steps_till_energy == steps_to_run:
                            save_energy = True
                            self._next_save_energy += energy_frequency_steps
                        else:
                            save_energy = False

                        nrun_till_save = min(steps_till_frame, steps_till_energy)

                        assert nrun_till_save >= 0

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

                        while nrun_till_save > 0:
                            nrun = block_size

                            if 2 * nrun > nrun_till_save:
                                nrun = nrun_till_save

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
                                nrun_till_save -= nrun
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
                            save_velocities=save_velocities,
                            delta_lambda=delta_lambda,
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
            # try to fix this problem by minimising,
            # then running again
            self._is_running = False
            self._clear_state()
            self._rebuild_and_minimise()
            orig_args["auto_fix_minimise"] = False
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
        save_velocities: bool = None,
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
            By default this is False. Set this to True if you are
            interested in saving the velocities.

        auto_fix_minimise: bool
            Whether or not to automatically run minimisation if the
            trajectory exits with an error in the first few steps.
            Such failures often indicate that the system needs
            minimsing. This automatically runs the minimisation
            in these cases, and then runs the requested dynamics.
        """
        if not self._d.is_null():
            if save_velocities is None:
                if self._d._map.specified("save_velocities"):
                    save_velocities = self._d._map["save_velocities"].value().as_bool()
                else:
                    save_velocities = False

            self._d.run(
                time=time,
                save_frequency=save_frequency,
                frame_frequency=frame_frequency,
                energy_frequency=energy_frequency,
                lambda_windows=lambda_windows,
                save_velocities=save_velocities,
                auto_fix_minimise=auto_fix_minimise,
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

    def set_lambda(self, lambda_value: float, update_constraints: bool = True):
        """
        Set the current value of lambda for this system. This will
        update the forcefield parameters in the context according
        to the data in the LambdaSchedule. This does nothing if
        this isn't a perturbable system.

        If `update_constraints` is True, then this will also update
        the constraint length of any constrained perturbable bonds.
        These will be set to the r0 value for that bond at this
        value of lambda. If `update_constraints` is False, then
        the constraint will not be changed.
        """
        self._d.set_lambda(
            lambda_value=lambda_value, update_constraints=update_constraints
        )

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

    def current_potential_energy(self, lambda_values=None):
        """
        Return the current potential energy.

        If `lambda_values` is passed (which should be a list of
        lambda values) then this will return the energies
        (as a list) at the requested lambda values
        """
        if lambda_values is None:
            return self._d.current_potential_energy()
        else:
            if not isinstance(lambda_values, list):
                lambda_values = [lambda_values]

            # save the current value of lambda so we
            # can restore it
            old_lambda = self.get_lambda()

            nrgs = []

            try:
                for lambda_value in lambda_values:
                    self.set_lambda(lambda_value, update_constraints=False)
                    nrgs.append(self._d.current_potential_energy())
            except Exception:
                self.set_lambda(old_lambda, update_constraints=False)
                raise

            self.set_lambda(old_lambda, update_constraints=False)

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
