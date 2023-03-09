__all__ = ["Dynamics"]


class DynamicsData:
    """
    This is a shared data class that holds all of the data for a
    dynamics object
    """

    def __init__(self, molecules=None, map=None):
        if molecules is None:
            self.sire_molecules = None
            self.omm_molecules = None
        else:
            from ..convert import to

            self.sire_molecules = molecules
            self.omm_molecules = to(molecules, "openmm", map)

        self.map = map


class Dynamics:
    """
    This class provides an interface for editing molecule(s) by performing
    dynamics-type moves
    """

    def __init__(self, molecules=None, map=None):
        self._d = DynamicsData(molecules=molecules, map=map)

    def minimise(self, max_steps=1000):
        return self

    def run(self, time):
        return self

    def commit(self):
        return self._d.sire_molecules
