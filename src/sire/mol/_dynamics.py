__all__ = ["Dynamics"]


class DynamicsData:
    """
    Internal class that is designed to only be used by the Dynamics
    class. This holds the shared state for dynamics on a set
    of molecule(s).
    """

    def __init__(self, mols=None, map=None):
        if mols is not None:
            self._map = map  # this is already a PropertyMap

            # get the forcefield info from the passed parameters
            # and from whatever we can glean from the molecules
            from ..system import ForceFieldInfo

            self._ffinfo = ForceFieldInfo(mols, map=self._map)

            # We want to store the molecules as a System so that
            # we can more easily track the space and time properties
            from ..system import System, ForceFieldInfo

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
            self._is_running = None

            from ..convert import to

            self._omm_mols = to(self._sire_mols, "openmm", map=self._map)
            self._omm_state = None
            self._omm_state_has_cv = False

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

    def is_null(self):
        return self._sire_mols is None

    def _clear_state(self):
        self._omm_state = None
        self._omm_state_has_cv = False

    def _update_from(self, state, nsteps_completed):
        if self.is_null():
            return

        from ..legacy.Convert import (
            openmm_extract_coordinates_and_velocities,
            openmm_extract_space,
        )

        mols = openmm_extract_coordinates_and_velocities(
            state, self._sire_mols.molecules(), self._map
        )

        space = openmm_extract_space(state)

        self._current_step = nsteps_completed

        self._sire_mols.update(mols.to_molecules())
        self._sire_mols.set_property("space", space)
        self._sire_mols.set_property("time", self._current_time)
        self._ffinfo.set_space(space)

    def _start_dynamics_block(self):
        if self._is_running is not None:
            raise SystemError(
                "Cannot start dynamics while it is already running!"
            )

        import datetime

        self._is_running = datetime.datetime.now().timestamp()
        self._omm_state = None
        self._omm_state_has_cv = False

    def _end_dynamics_block(self, coords_and_vels=False):
        if self._is_running is None:
            raise SystemError("Cannot stop dynamics that is not running!")

        import datetime
        import openmm
        from ..units import second, nanosecond

        self._walltime += (
            datetime.datetime.now().timestamp() - self._is_running
        ) * second

        if coords_and_vels:
            self._omm_state = self._omm_mols.getState(
                getEnergy=True,
                getPositions=True,
                getVelocities=True,
                enforcePeriodicBox=self._enforce_periodic_box,
            )
        else:
            self._omm_state = self._omm_mols.getState(getEnergy=True)

        self._omm_state_has_cv = coords_and_vels

        current_time = (
            self._omm_state.getTime().value_in_unit(openmm.unit.nanosecond)
            * nanosecond
        )

        delta = current_time - self._elapsed_time

        self._elapsed_time = current_time
        self._current_time += delta

        self._is_running = None

        return self._omm_state

    def _get_current_state(self, coords_and_vels=False):
        if self._omm_state is not None:
            if self._omm_state_has_cv or (coords_and_vels is False):
                return self._omm_state

        # we need to get this again as either it doesn't exist,
        # or we now want velocities as well
        if coords_and_vels:
            self._omm_state = self._omm_mols.getState(
                getEnergy=True,
                getPositions=True,
                getVelocities=True,
                enforcePeriodicBox=self._enforce_periodic_box,
            )
        else:
            self._omm_state = self._omm_mols.getState(getEnergy=True)

        self._omm_state_has_cv = coords_and_vels

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

    def constraint(self):
        if self.is_null():
            return None
        else:
            if self._map.specified("constraint"):
                return self._map["constraint"].source()
            else:
                return "none"

    def info(self):
        if self.is_null():
            return None
        else:
            return self._ffinfo

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

    def _rebuild_and_minimise(self):
        if self.is_null():
            return

        from concurrent.futures import ThreadPoolExecutor
        from ..utils import Console
        import openmm

        def runfunc(max_its):
            openmm.LocalEnergyMinimizer.minimize(
                self._omm_mols, maxIterations=max_its
            )

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

        with Console.spinner("minimisation") as spinner:
            with ThreadPoolExecutor() as pool:
                run_promise = pool.submit(runfunc, 0)

                while not run_promise.done():
                    try:
                        run_promise.result(timeout=1.0)
                    except Exception:
                        pass

            spinner.success()

    def run(self, time, save_frequency, auto_fix_minimise: bool = True):
        if self.is_null():
            return

        from concurrent.futures import ThreadPoolExecutor
        import openmm

        from ..utils import Console
        from ..units import picosecond

        try:
            steps = int(time.to(picosecond) / self.timestep().to(picosecond))
        except Exception:
            # passed in the number of steps instead
            steps = int(time)

        if steps < 1:
            steps = 1

        if steps <= 0:
            return

        def runfunc(num_steps):
            try:
                integrator = self._omm_mols.getIntegrator()
                integrator.step(num_steps)
                return 0
            except Exception as e:
                return e

        def process_block(state, nsteps_completed):
            self._update_from(state, nsteps_completed)
            self._sire_mols.save_frame(
                map={
                    "space": self.current_space(),
                    "time": self.current_time(),
                }
            )

        completed = 0

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

        if save_frequency <= 0:
            # don't save any intermediate frames
            save_size = steps
        else:
            save_size = int(save_frequency / self.timestep().to(picosecond))

        block_size = 50

        state = None
        saved_last_frame = False

        class NeedsMinimiseError(Exception):
            pass

        nsteps_before_run = self._current_step

        try:
            with Console.progress(transient=True) as progress:
                t = progress.add_task("dynamics", total=steps)

                with ThreadPoolExecutor() as pool:
                    while completed < steps:
                        if steps - completed < save_size:
                            nrun_till_save = steps - completed
                        else:
                            nrun_till_save = save_size

                        self._start_dynamics_block()

                        # process the last block in the foreground
                        if state is not None:
                            process_promise = pool.submit(
                                process_block,
                                state,
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
                                progress.update(t, completed=completed)
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
                        state = self._end_dynamics_block(coords_and_vels=True)
                        saved_last_frame = False

                        kinetic_energy = (
                            state.getKineticEnergy().value_in_unit(
                                openmm.unit.kilocalorie_per_mole
                            )
                        )

                        ke_per_atom = kinetic_energy / self._num_atoms

                        if ke_per_atom > 1000:
                            # The system has blown up!
                            state = None
                            saved_last_frame = True

                            if completed == 0 and auto_fix_minimise:
                                raise NeedsMinimiseError()

                            raise RuntimeError(
                                "The kinetic energy has exceeded 100 kcal mol-1 "
                                f"per atom (it is {ke_per_atom} kcal mol-1 atom-1,"
                                f" and {kinetic_energy} kcal mol-1 total). This "
                                "suggests that the simulation has become "
                                "unstable. Try reducing the timestep and/or "
                                "minimising the system and run again."
                            )

            if state is not None and not saved_last_frame:
                # we can process the last block in the main thread
                process_block(state, nsteps_before_run + completed)

        except NeedsMinimiseError:
            # try to fix this problem by minimising,
            # then running again
            self._is_running = None
            self._omm_state = None
            self._omm_state_has_cv = False
            self._rebuild_and_minimise()
            self.run(
                time=time,
                save_frequency=save_frequency * picosecond,
                auto_fix_minimise=False,
            )
            return

    def commit(self):
        if self.is_null():
            return

        self._update_from(
            self._get_current_state(coords_and_vels=True), self._current_step
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
        save_frequency=None,
        constraint=None,
    ):
        from ..base import create_map

        extras = {}

        _add_extra(extras, "cutoff", cutoff)
        _add_extra(extras, "cutoff_type", cutoff_type)
        _add_extra(extras, "timestep", timestep)
        _add_extra(extras, "save_frequency", save_frequency)
        _add_extra(extras, "constraint", constraint)

        map = create_map(map, extras)

        self._d = DynamicsData(mols=mols, map=map)

        # Save the original view too, as this is the view that
        # will be returned to the user
        self._view = mols

    def __str__(self):
        speed = self.time_speed()

        if speed == 0:
            return f"Dynamics(completed=0)"
        else:
            return (
                f"Dynamics(completed={self.current_time()}, "
                f"energy={self.current_energy()}, speed={speed:.1f} ns day-1)"
            )

    def __repr__(self):
        return self.__str__()

    def run(self, time, save_frequency=None):
        """
        Perform dynamics on the molecules
        """
        if not self._d.is_null():
            self._d.run(time=time, save_frequency=save_frequency)

        return self

    def ensemble(self):
        """
        Return the ensemble in which the simulation is being performed
        """
        return self._d.ensemble()

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

    def current_potential_energy(self):
        """
        Return the current potential energy
        """
        return self._d.current_potential_energy()

    def current_kinetic_energy(self):
        """
        Return the current kinetic energy
        """
        return self._d.current_kinetic_energy()

    def commit(self):
        if not self._d.is_null():
            self._d.commit()

        return self._d._sire_mols
