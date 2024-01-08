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
        coulomb_power=None,
        restraints=None,
        fixed=None,
    ):
        from ..base import create_map
        from ._dynamics import DynamicsData, _add_extra

        extras = {}

        _add_extra(extras, "cutoff", cutoff)
        _add_extra(extras, "cutoff_type", cutoff_type)
        _add_extra(extras, "schedule", schedule)
        _add_extra(extras, "lambda", lambda_value)
        _add_extra(extras, "swap_end_states", swap_end_states)
        _add_extra(extras, "ignore_perturbations", ignore_perturbations)
        _add_extra(extras, "shift_delta", shift_delta)
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

    def run(self, max_iterations: int = 10000):
        """
        Perform minimisation on the molecules, running a maximum
        of max_iterations iterations.

        Parameters:

        - max_iterations (int): The maximum number of iterations to run
        """
        if not self._d.is_null():
            self._d.run_minimisation(max_iterations=max_iterations)

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

    def __call__(self, max_iterations: int = 10000):
        """
        Perform minimisation on the molecules, running a maximum
        of max_iterations iterations.

        Parameters:

        - max_iterations (int): The maximum number of iterations to run
        """
        return self.run(max_iterations=max_iterations).commit()
