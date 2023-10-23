__all__ = [
    "Integrator",
    "Constraint",
    "PerturbableConstraint",
    "Cutoff",
    "Platform",
]

from ._option import Option as _Option


class Integrator(_Option):
    """
    All of the supported options for the integrator
    """

    VERLET = "verlet"
    LEAPFROG = "leapfrog"
    LANGEVIN = "langevin"
    LANGEVIN_MIDDLE = "langevin_middle"
    NOSE_HOOVER = "nose_hoover"
    BROWNIAN = "brownian"
    ANDERSEN = "andersen"

    @staticmethod
    def create(option: str):
        return _Option._create(Integrator, option)

    @staticmethod
    def options():
        return _Option._options(Integrator)


class Constraint(_Option):
    """
    All of the supported constraint options
    """

    NONE = "none"
    AUTO = "auto"
    HBONDS = "h-bonds"
    BONDS = "bonds"
    HBONDS_HANGLES = "h-bonds-h-angles"
    BOND_HANGLES = "bonds-h-angles"

    @staticmethod
    def create(option: str):
        return _Option._create(Constraint, option)

    @staticmethod
    def options():
        return _Option._options(Constraint)


PerturbableConstraint = Constraint


class Cutoff(_Option):
    """
    All of the support cutoff options
    """

    NONE = "none"
    RF = "rf"
    PME = "pme"
    EWALD = "ewald"

    @staticmethod
    def canonicalise(option: str):
        """
        Convert the passed option string to the canonical form
        """
        option = _Option.canonicalise(option).replace("-", "_")

        if option == "reaction_field":
            return "rf"
        elif option == "particle_mesh_ewald":
            return "pme"
        elif option == "no_cutoff":
            return "none"
        else:
            return option

    @staticmethod
    def create(option: str):
        return _Option._create(Cutoff, option)

    @staticmethod
    def options():
        return _Option._options(Cutoff)


class Platform(_Option):
    """
    All of the supported platforms
    """

    AUTO = "auto"
    CPU = "cpu"
    CUDA = "cuda"
    OPENCL = "opencl"
    METAL = "metal"
    HIP = "hip"

    @staticmethod
    def create(option: str):
        return _Option._create(Platform, option)

    @staticmethod
    def options():
        return _Option._options(Platform)
