__all__ = ["Minimisation"]


class MinimisationData:
    """
    Internal class that is designed to only be used by the Minimisation
    class. This holds the shared state for minimisation on a set
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

    def run(self, max_iterations: int):
        from openmm import LocalEnergyMinimizer
        from concurrent.futures import ThreadPoolExecutor
        import time

        if max_iterations <= 0:
            max_iterations = 0

        from ..utils import Console

        def runfunc(max_its):
            LocalEnergyMinimizer.minimize(
                self._omm_mols, maxIterations=max_its
            )

        with Console.spinner("Minimisation") as spinner:
            with ThreadPoolExecutor() as pool:
                r = pool.submit(runfunc, max_iterations)

                while not r.done():
                    time.sleep(0.2)

            spinner.success()

    def commit(self):
        from ..legacy.Convert import openmm_extract_coordinates

        state = self._omm_mols.getState(getPositions=True)

        mols = openmm_extract_coordinates(
            state, self._sire_mols.molecules(), self._map
        )

        self._sire_mols.update(mols.to_molecules())


class Minimisation:
    """
    Class that runs minimisation on the contained molecule(s). Note that
    this class is not designed to be constructed directly. You should only
    use this class by calling `.minimisation()` on the molecules(s)
    you want to minimise
    """

    def __init__(self, mols=None, map=None):
        self._d = MinimisationData(mols=mols, map=map)

    def run(self, max_iterations: int = 1000):
        """
        Perform minimisation on the molecules, running a maximum
        of max_iterations iterations.
        """
        if not self._d.is_null():
            self._d.run(max_iterations=max_iterations)

        return self

    def commit(self):
        if not self._d.is_null():
            self._d.commit()

        return self._d._sire_mols
