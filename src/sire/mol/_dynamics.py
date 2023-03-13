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

        else:
            self._sire_mols = None
            self._map = None
            self._omm_mols = None

    def is_null(self):
        return self._sire_mols is None

    def run(self):
        from concurrent.futures import ThreadPoolExecutor
        import time

        from ..utils import Console

        def runfunc(num_steps):
            print("RUN", num_steps)
            self._omm_mols.integrator.step(num_steps)

        with Console.spinner("Dynamics") as spinner:
            with ThreadPoolExecutor() as pool:
                r = pool.submit(runfunc, 100)

                while not r.done():
                    time.sleep(0.2)

            spinner.success()

    def commit(self):
        from ..legacy.Convert import openmm_extract_coordinates_and_velocities

        state = self._omm_mols.getState(getPositions=True, getVelocities=True)

        mols = openmm_extract_coordinates_and_velocities(
            state, self._sire_mols.molecules(), self._map
        )

        self._sire_mols.update(mols.to_molecules())


class Dynamics:
    """
    Class that runs dynamics on the contained molecule(s). Note that
    this class is not designed to be constructed directly. You should only
    use this class by calling `.dynamics()` on the molecules(s)
    you want to simulate
    """

    def __init__(self, mols=None, map=None):
        self._d = DynamicsData(mols=mols, map=map)

    def run(self):
        """
        Perform dynamics on the molecules
        """
        if not self._d.is_null():
            self._d.run()

        return self

    def commit(self):
        if not self._d.is_null():
            self._d.commit()

        return self._d._sire_mols
