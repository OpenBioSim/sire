__all__ = ["Minimisation"]

from ._dynamics import DynamicsData as _DynamicsData
from ._dynamics import _add_extra


class MinimisationData(_DynamicsData):
    """
    Internal class that is designed to only be used by the Minimisation
    class. This holds the shared state for minimisation on a set
    of molecule(s).
    """

    def __init__(self, mols=None, map=None):
        super().__init__(mols=mols, map=map)

    def run(self, max_iterations: int):
        from openmm import LocalEnergyMinimizer
        from concurrent.futures import ThreadPoolExecutor

        if max_iterations <= 0:
            max_iterations = 0

        from ..base import ProgressBar

        def runfunc(max_its):
            LocalEnergyMinimizer.minimize(
                self._omm_mols, maxIterations=max_its
            )

        with ProgressBar(text="minimisation") as spinner:
            spinner.set_speed_unit("checks / s")

            with ThreadPoolExecutor() as pool:
                run_promise = pool.submit(runfunc, max_iterations)

                while not run_promise.done():
                    try:
                        run_promise.result(timeout=0.2)
                    except Exception:
                        spinner.tick()
                        pass


class Minimisation:
    """
    Class that runs minimisation on the contained molecule(s). Note that
    this class is not designed to be constructed directly. You should only
    use this class by calling `.minimisation()` on the molecules(s)
    you want to minimise
    """

    def __init__(
        self,
        mols=None,
        map=None,
        cutoff=None,
        cutoff_type=None,
        schedule=None,
        lambda_value=None,
    ):
        from ..base import create_map

        extras = {}

        _add_extra(extras, "cutoff", cutoff)
        _add_extra(extras, "cutoff_type", cutoff_type)
        _add_extra(extras, "schedule", schedule)
        _add_extra(extras, "lambda", lambda_value)

        map = create_map(map, extras)

        self._d = MinimisationData(mols=mols, map=map)

    def __str__(self):
        return f"Minimisation()"

    def __repr__(self):
        return self.__str__()

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
