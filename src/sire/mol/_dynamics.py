__all__ = ["Dynamics"]


class DynamicsData:
    """
    Internal class that is designed to only be used by the Dynamics
    class. This holds the shared state for dynamics on a set
    of molecule(s).
    """

    def __init__(self, mols=None, map=None):
        if mols is not None:
            # eventually want to call 'extract' on this?
            self._sire_mols = mols

            from ..base import create_map

            self._map = create_map(map)

            from ..convert import to

            self._omm_mols = to(self._sire_mols, "openmm", map=self._map)
            self._omm_state = None
            self._omm_state_has_cv = False

        else:
            self._sire_mols = None
            self._map = None
            self._omm_mols = None
            self._omm_state = None
            self._omm_state_has_cv = False

    def is_null(self):
        return self._sire_mols is None

    def _clear_state(self):
        self._omm_state = None
        self._omm_state_has_cv = False

    def _get_current_state(self, coords_and_vels=False):
        if self._omm_state is not None:
            if self._omm_state_has_cv or (coords_and_vels == False):
                return self._omm_state

        if self._omm_mols is None:
            return None

        if coords_and_vels:
            self._omm_state = self._omm_mols.getState(
                getEnergy=True, getPositions=True, getVelocities=True
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

    def current_step(self):
        if self.is_null():
            return 0
        else:
            return self._omm_mols.getStepCount()

    def current_time(self):
        if self.is_null():
            return 0
        else:
            from openmm.unit import picosecond as _omm_ps
            from ..units import picosecond as _sire_ps

            return self._omm_mols.getTime().value_in_unit(_omm_ps) * _sire_ps

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

    def run(self, steps: int):
        if self.is_null():
            return

        from concurrent.futures import ThreadPoolExecutor
        import time

        from ..utils import Console

        if steps <= 0:
            return

        def runfunc(num_steps):
            try:
                integrator = self._omm_mols.getIntegrator()
                integrator.step(num_steps)
                return 0
            except Exception as e:
                return e

        completed = 0

        if self._map.specified("block_size"):
            block_size = self._map["block_size"].value().as_integer()

            if block_size < 1:
                block_size = 1
        else:
            block_size = 100

        state = None
        saved_last_frame = False

        with Console.progress() as progress:
            t = progress.add_task("dynamics", steps)

            with ThreadPoolExecutor() as pool:
                while completed < steps:
                    if steps - completed < block_size:
                        nrun = steps - completed
                    else:
                        nrun = block_size

                    self._clear_state()
                    # run the current block in the background
                    r = pool.submit(runfunc, nrun)

                    # process the last block in the foreground
                    if state is not None:
                        self._update_from(state)
                        self._sire_mols.save_frame()
                        saved_last_frame = True

                    while not r.done():
                        try:
                            result = r.result(timeout=0.5)
                        except Exception:
                            pass

                    if result == 0:
                        completed += nrun
                        progress.update(t, completed=completed)
                    else:
                        # something went wrong - re-raise the exception
                        raise result

                    # get the state, including coordinates and velocities
                    state = self._get_current_state(coords_and_vels=True)

        if state is not None and not saved_last_frame:
            self._update_from(state)
            self._sire_mols.save_frame()

    def _update_from(self, state):
        if self.is_null():
            return

        from ..legacy.Convert import openmm_extract_coordinates_and_velocities

        mols = openmm_extract_coordinates_and_velocities(
            state, self._sire_mols.molecules(), self._map
        )

        self._sire_mols.update(mols.to_molecules())

    def commit(self):
        if self.is_null():
            return

        self._update_from(self._get_current_state(coords_and_vels=True))


class Dynamics:
    """
    Class that runs dynamics on the contained molecule(s). Note that
    this class is not designed to be constructed directly. You should only
    use this class by calling `.dynamics()` on the molecules(s)
    you want to simulate
    """

    def __init__(self, mols=None, map=None):
        self._d = DynamicsData(mols=mols, map=map)

    def run(self, steps: int):
        """
        Perform dynamics on the molecules
        """
        if not self._d.is_null():
            self._d.run(steps=steps)

        return self

    def ensemble(self):
        """
        Return the ensemble in which the simulation is being performed
        """
        return self._d.ensemble()

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
