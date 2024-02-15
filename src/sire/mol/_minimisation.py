__all__ = ["Minimisation"]


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
        swap_end_states=None,
        ignore_perturbations=None,
        shift_delta=None,
        shift_coulomb=None,
        coulomb_power=None,
        restraints=None,
        fixed=None,
    ):
        from ..base import create_map
        from ._dynamics import DynamicsData, _add_extra
        from .. import u

        extras = {}

        _add_extra(extras, "cutoff", cutoff)
        _add_extra(extras, "cutoff_type", cutoff_type)
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

        map = create_map(map, extras)

        self._d = DynamicsData(mols=mols, map=map)

    def __str__(self):
        return "Minimisation()"

    def __repr__(self):
        return self.__str__()

    def constraint(self):
        """
        Return the constraint used for the minimisation (e.g. constraining
        bonds involving hydrogens etc.)
        """
        return self._d.constraint()

    def perturbable_constraint(self):
        """
        Return the perturbable constraint used for the minimisation (e.g.
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

    def get_log(self):
        """
        Return the log of the minimisation
        """
        return self._d.get_minimisation_log()

    def run(
        self,
        max_iterations: int = 10000,
        tolerance: float = 10.0,
        max_restarts: int = 10,
        max_ratchets: int = 20,
        ratchet_frequency: int = 500,
        starting_k: float = 100.0,
        ratchet_scale: float = 10.0,
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
            )

        return self

    def commit(self, return_as_system: bool = False):
        """
        Commit the minimisation to the molecules, returning the
        minimised molecules.

        Normally this will return the same view of as was used for
        construction. If `return_as_system` is True, then this will
        return a System object instead.
        """
        if not self._d.is_null():
            return self._d.commit(return_as_system=return_as_system)
        else:
            return None

    def __call__(self, *args, **kwargs):
        """
        Perform minimisation on the molecules, running a maximum
        of max_iterations iterations.

        Parameters:

        max_iterations: int = 10000,
        tolerance: float = 10.0,
        max_restarts: int = 10,
        max_ratchets: int = 20,
        starting_k: float = 100.0,
        """
        return self.run(*args, **kwargs).commit()
