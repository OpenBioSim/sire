__all__ = ["Dynamics"]


from typing import Any


class DynamicsData:
    """
    Internal class that is designed to only be used by the Dynamics
    class. This holds the shared state for dynamics on a set
    of molecule(s).
    """

    def __init__(self, mols=None, map=None, **kwargs):
        from ..base import create_map

        map = create_map(map, kwargs)

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
                map.set(
                    "fixed",
                    mols.atoms().find(selection_to_atoms(mols, fixed_atoms)),
                )

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
                self._sire_mols._system.add(
                    mols.molecules().to_molecule_group()
                )
                self._sire_mols._system.set_property(
                    "space", self._ffinfo.space()
                )

            # find the existing energy trajectory - we will build on this
            self._energy_trajectory = self._sire_mols.energy_trajectory(
                to_pandas=False, map=self._map
            )

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

            from ..convert import to

            self._omm_mols = to(self._sire_mols, "openmm", map=self._map)
            self._omm_state = None
            self._omm_state_has_cv = (False, False)

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

        if not state_has_cv[0]:
            # there is no information to update
            return

        from ..legacy.Convert import (
            openmm_extract_coordinates_and_velocities,
            openmm_extract_coordinates,
            openmm_extract_space,
        )

        if state_has_cv[1]:
            # get velocities too
            mols = openmm_extract_coordinates_and_velocities(
                state,
                self._sire_mols.molecules(),
                perturbable_maps=self._omm_mols.get_lambda_lever().get_perturbable_molecule_maps(),
                map=self._map,
            )
        else:
            mols = openmm_extract_coordinates(
                state,
                self._sire_mols.molecules(),
                perturbable_maps=self._omm_mols.get_lambda_lever().get_perturbable_molecule_maps(),
                map=self._map,
            )

        space = openmm_extract_space(state)

        self._current_step = nsteps_completed

        self._sire_mols.update(mols.to_molecules())
        self._sire_mols.set_property("space", space)
        self._sire_mols.set_property("time", self._current_time)
        self._ffinfo.set_space(space)

    def _enter_dynamics_block(self):
        if self._is_running:
            raise SystemError(
                "Cannot start dynamics while it is already running!"
            )

        self._is_running = True
        self._omm_state = None
        self._omm_state_has_cv = (False, False)

    def _exit_dynamics_block(
        self,
        save_frame: bool = False,
        save_energy: bool = False,
        lambda_windows=[],
        save_velocities: bool = False,
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
            self._omm_state.getTime().value_in_unit(openmm.unit.nanosecond)
            * nanosecond
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
            nrgs[str(sim_lambda_value)] = nrgs["potential"]

            if lambda_windows is not None:
                for lambda_value in lambda_windows:
                    if lambda_value != sim_lambda_value:
                        self._omm_mols.set_lambda(lambda_value)
                        nrgs[str(lambda_value)] = (
                            self._omm_mols.get_potential_energy(
                                to_sire_units=False
                            ).value_in_unit(openmm.unit.kilocalorie_per_mole)
                            * kcal_per_mol
                        )

                self._omm_mols.set_lambda(sim_lambda_value)

            self._energy_trajectory.set(
                self._current_time, nrgs, {"lambda": str(sim_lambda_value)}
            )

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

    def set_ensemble(self, ensemble, rescale_velocities: bool = True):
        """
        Set the ensemble for the dynamics. Note that this will
        only let you change the temperature and/or pressure of the
        ensemble. You can't change its fundemental nature.

        If rescalse_velocities is True, then the velocities will
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

    def constraint(self):
        if self.is_null():
            return None
        else:
            if self._map.specified("constraint"):
                return self._map["constraint"].source()
            else:
                return "none"

    def get_schedule(self):
        if self.is_null():
            return None
        else:
            return self._omm_mols.get_lambda_schedule()

    def set_schedule(self, schedule):
        if not self.is_null():
            self._omm_mols.set_lambda_schedule(schedule)

    def get_lambda(self):
        if self.is_null():
            return None
        else:
            return self._omm_mols.get_lambda()

    def set_lambda(self, lambda_value: float):
        if not self.is_null():
            s = self.get_schedule()

            if s is None:
                return

            lambda_value = s.clamp(lambda_value)

            if lambda_value == self._omm_mols.get_lambda():
                # nothing to do
                return

            self._omm_mols.set_lambda(lambda_value)
            self._clear_state()

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

    def run_minimisation(self, max_iterations: int):
        """
        Internal method that runs minimisation on the molecules.

        Parameters:

        - max_iterations (int): The maximum number of iterations to run
        """
        from openmm import LocalEnergyMinimizer
        from concurrent.futures import ThreadPoolExecutor

        if max_iterations <= 0:
            max_iterations = 0

        from ..base import ProgressBar

        def runfunc(max_its):
            try:
                LocalEnergyMinimizer.minimize(
                    self._omm_mols, maxIterations=max_its
                )

                return 0
            except Exception as e:
                return e

        with ProgressBar(text="minimisation") as spinner:
            spinner.set_speed_unit("checks / s")

            with ThreadPoolExecutor() as pool:
                run_promise = pool.submit(runfunc, max_iterations)

                while not run_promise.done():
                    try:
                        result = run_promise.result(timeout=0.2)
                    except Exception:
                        spinner.tick()
                        pass

                if result != 0:
                    raise result

    def _rebuild_and_minimise(self):
        if self.is_null():
            return

        from ..utils import Console

        Console.warning(
            "Something went wrong when running dynamics. Since no steps "
            "were completed, it is likely that the system needs minimising. "
            "The system will be minimised, and then dynamics will be "
            "attempted again. If an error still occurs, then it is likely "
            "that the step size is too large, the molecules are "
            "over-constrained, or there is something more fundemental "
            "going wrong..."
        )

        # rebuild the molecules
        from ..convert import to

        self._omm_mols = to(self._sire_mols, "openmm", map=self._map)

        self.run_minimisation(max_iterations=10000)

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

        if lambda_windows is not None:
            if type(lambda_windows) is not list:
                lambda_windows = [lambda_windows]

        try:
            steps_to_run = int(
                time.to(picosecond) / self.timestep().to(picosecond)
            )
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
                    save_frequency = (
                        self._map["save_frequency"].value().to(picosecond)
                    )
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

        if lambda_windows is None:
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

        completed = 0

        frame_frequency_steps = int(
            frame_frequency / self.timestep().to(picosecond)
        )

        energy_frequency_steps = int(
            energy_frequency / self.timestep().to(picosecond)
        )

        def get_steps_till_save(completed: int, total: int):
            """Internal function to calculate the number of steps
            to run before the next save. This returns a tuple
            of number of steps, and then if a frame should be
            saved and if the energy should be saved
            """
            if completed < 0:
                completed = 0

            if completed >= total:
                return (0, True, True)

            elif frame_frequency_steps <= 0 and energy_frequency_steps <= 0:
                return (total, True, True)

            n_to_end = total - completed

            if frame_frequency_steps > 0:
                n_to_frame = min(
                    frame_frequency_steps
                    - (completed % frame_frequency_steps),
                    n_to_end,
                )
            else:
                n_to_frame = total - completed

            if energy_frequency_steps > 0:
                n_to_energy = min(
                    energy_frequency_steps
                    - (completed % energy_frequency_steps),
                    n_to_end,
                )
            else:
                n_to_energy = total - completed

            if n_to_frame == n_to_energy:
                return (n_to_frame, True, True)
            elif n_to_frame < n_to_energy:
                return (n_to_frame, True, False)
            else:
                return (n_to_energy, False, True)

        block_size = 50

        state = None
        state_has_cv = (False, False)
        saved_last_frame = False

        class NeedsMinimiseError(Exception):
            pass

        nsteps_before_run = self._current_step

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
                        (
                            nrun_till_save,
                            save_frame,
                            save_energy,
                        ) = get_steps_till_save(completed, steps_to_run)

                        assert nrun_till_save > 0

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

                            while not run_promise.done():
                                try:
                                    result = run_promise.result(timeout=1.0)
                                except Exception:
                                    pass

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
                                    except Exception:
                                        pass

                                if completed == 0 and auto_fix_minimise:
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
                        )

                        saved_last_frame = False

                        kinetic_energy = (
                            state.getKineticEnergy().value_in_unit(
                                openmm.unit.kilocalorie_per_mole
                            )
                        )

                        ke_per_atom = kinetic_energy / self._num_atoms

                        if isnan(ke_per_atom) or ke_per_atom > 1000:
                            # The system has blown up!
                            state = None
                            saved_last_frame = True

                            if completed == 0 and auto_fix_minimise:
                                raise NeedsMinimiseError()

                            raise RuntimeError(
                                "The kinetic energy has exceeded 1000 kcal mol-1 "
                                f"per atom (it is {ke_per_atom} kcal mol-1 atom-1,"
                                f" and {kinetic_energy} kcal mol-1 total). This "
                                "suggests that the simulation has become "
                                "unstable. Try reducing the timestep and/or "
                                "minimising the system and run again."
                            )

                self._walltime += (
                    datetime.now() - start_time
                ).total_seconds() * second

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
            self._omm_state = None
            self._omm_state_has_cv = (False, False)
            self._rebuild_and_minimise()
            orig_args["auto_fix_minimise"] = False
            self.run(**orig_args)
            return

    def commit(self):
        if self.is_null():
            return

        self._update_from(
            state=self._get_current_state(
                include_coords=True, include_velocities=True
            ),
            state_has_cv=(True, True),
            nsteps_completed=self._current_step,
        )

        self._sire_mols.set_energy_trajectory(
            self._energy_trajectory, map=self._map
        )


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
        schedule=None,
        lambda_value=None,
        swap_end_states=None,
        ignore_perturbations=None,
        shift_delta=None,
        coulomb_power=None,
        restraints=None,
        fixed=None,
    ):
        from ..base import create_map
        from .. import u

        extras = {}

        _add_extra(extras, "cutoff", cutoff)
        _add_extra(extras, "cutoff_type", cutoff_type)

        if timestep is not None:
            _add_extra(extras, "timestep", u(timestep))

        _add_extra(extras, "constraint", constraint)
        _add_extra(extras, "schedule", schedule)
        _add_extra(extras, "lambda", lambda_value)
        _add_extra(extras, "swap_end_states", swap_end_states)
        _add_extra(extras, "ignore_perturbations", ignore_perturbations)

        if shift_delta is not None:
            _add_extra(extras, "shift_delta", u(shift_delta))

        _add_extra(extras, "coulomb_power", coulomb_power)
        _add_extra(extras, "restraints", restraints)
        _add_extra(extras, "fixed", fixed)

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

    def minimise(self, max_iterations: int = 10000):
        """
        Perform minimisation on the molecules, running a maximum
        of max_iterations iterations.

        Parameters:

        - max_iterations (int): The maximum number of iterations to run
        """
        if not self._d.is_null():
            self._d.run_minimisation(max_iterations=max_iterations)

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
                    save_velocities = (
                        self._d._map["save_velocities"].value().as_bool()
                    )
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

    def set_lambda(self, lambda_value: float):
        """
        Set the current value of lambda for this system. This will
        update the forcefield parameters in the context according
        to the data in the LambdaSchedule. This does nothing if
        this isn't a perturbable system
        """
        self._d.set_lambda(lambda_value)

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

    def constraint(self):
        """
        Return the constraint used for the dynamics (e.g. constraining
        bonds involving hydrogens etc.)
        """
        return self._d.constraint()

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
            if not type(lambda_values) is list:
                lambda_values = [lambda_values]

            # save the current value of lambda so we
            # can restore it
            old_lambda = self.get_lambda()

            nrgs = []

            try:
                for lambda_value in lambda_values:
                    self.set_lambda(lambda_value)
                    nrgs.append(self._d.current_potential_energy())
            except Exception:
                self.set_lambda(old_lambda)
                raise

            self.set_lambda(old_lambda)

            return nrgs

    def current_kinetic_energy(self):
        """
        Return the current kinetic energy
        """
        return self._d.current_kinetic_energy()

    def energy_trajectory(self, to_pandas: bool = True):
        """
        Return the energy trajectory. This is the trajectory of
        energy values that have been captured during dynamics.

        If 'to_pandas' is True, (the default) then this will
        be returned as a pandas dataframe, with times and energies
        in the defined default units
        """
        t = self._d.energy_trajectory()

        if to_pandas:
            return t.to_pandas()
        else:
            return t

    def commit(self):
        if not self._d.is_null():
            self._d.commit()

        return self._d._sire_mols

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
