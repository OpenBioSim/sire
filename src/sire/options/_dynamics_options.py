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

    AUTO = "auto", "Choose the integrator automatically"
    VERLET = "verlet", "Use the Verlet integrator"
    LEAPFROG = "leapfrog", "Use the Leapfrog integrator"
    LANGEVIN = "langevin", "Use the Langevin integrator"
    LANGEVIN_MIDDLE = (
        "langevin_middle",
        "Use the middle scheme Langevin integrator",
    )
    NOSE_HOOVER = "nose_hoover", "Use the Nose-Hoover integrator"
    BROWNIAN = "brownian", "Use the Brownian integrator"
    ANDERSEN = (
        "andersen",
        "Use the Verlet integrator with an Andersen thermostat",
    )

    @staticmethod
    def create(option: str):
        return _Option._create(Integrator, option)

    @staticmethod
    def options(include_docs: bool = False):
        return _Option._options(Integrator, include_docs=include_docs)


class Constraint(_Option):
    """
    All of the supported constraint options
    """

    NONE = "none", "Do not use constraints"
    AUTO = "auto", "Choose the constraints automatically"
    HBONDS = "h_bonds", "Constrain bonds involving hydrogens"
    HBONDS_NOT_PERTURBED = (
        "h_bonds_not_perturbed",
        "Constrain bonds involving hydrogens, excluding those that are perturbed",
    )
    HBONDS_NOT_HEAVY_PERTURBED = (
        "h_bonds_not_heavy_perturbed",
        "Constrain bonds involving hydrogens, excluding those that are perturbed "
        "but do not involve a hydrogen in any end state.",
    )
    BONDS = "bonds", "Constrain all bonds"
    BONDS_NOT_PERTURBED = (
        "bonds_not_perturbed",
        "Constrain all bonds, excluding those are perturbed",
    )
    BONDS_NOT_HEAVY_PERTURBED = (
        "bonds_not_heavy_perturbed",
        "Constrain all bonds, excluding those that are perturbed but do not "
        "involve a hydrogen in any end state.",
    )
    HBONDS_HANGLES = (
        "h_bonds_h_angles",
        "Constrain bonds and angles involving hydrogens",
    )
    HBONDS_HANGLES_NOT_PERTURBED = (
        "h_bonds_h_angles_not_perturbed",
        "Constrain bonds and angles involving hydrogens, "
        "excluding those that are perturbed.",
    )
    HBONDS_HANGLES_NOT_HEAVY_PERTURBED = (
        "h_bonds_h_angles_not_heavy_perturbed",
        "Constrain bonds and angles involving hydrogens, "
        "excluding those that are perturbed "
        "but do not involve a hydrogen in any end state.",
    )
    BOND_HANGLES = (
        "bonds_h_angles",
        "Constrain all bonds, and angles involving hydrogens",
    )
    BOND_HANGLES_NOT_PERTURBED = (
        "bonds_h_angles_not_perturbed",
        "Constrain all bonds, and angles involving hydrogens, "
        "excluding those that are perturbed",
    )
    BONDS_HANGLES_NOT_HEAVY_PERTURBED = (
        "bonds_h_angles_not_heavy_perturbed",
        "Constrain all bonds, and angles involving hydrogens, "
        "excluding those that are perturbed "
        "but do not involve a hydrogen in any end state.",
    )
    AUTO_BONDS = (
        "auto_bonds",
        "Choose the constraints automatically, constraining bonds based "
        "on whether their predicted vibrational periods are less than a "
        "tenth of the simulaton timestep.",
    )

    @staticmethod
    def create(option: str):
        return _Option._create(Constraint, option)

    @staticmethod
    def options(include_docs: bool = False):
        return _Option._options(Constraint, include_docs=include_docs)


PerturbableConstraint = Constraint


class Cutoff(_Option):
    """
    All of the support cutoff options
    """

    NONE = "none", "Do not use a cutoff"
    AUTO = "auto", "Choose the cutoff automatically"
    RF = "rf", "Use a reaction field cutoff"
    PME = "pme", "Use a Particle Mesh Ewald cutoff"
    EWALD = "ewald", "Use an Ewald cutoff"

    @staticmethod
    def canonicalise(option: str):
        """
        Convert the passed option string to the canonical form
        """
        option = _Option.canonicalise(option)

        if option == "reaction_field" or option == "reaction field":
            return "rf"
        elif option == "particle_mesh_ewald" or option == "particle mesh ewald":
            return "pme"
        elif option == "no_cutoff" or option == "no cutoff":
            return "none"
        else:
            return option

    @staticmethod
    def create(option: str):
        return _Option._create(Cutoff, option)

    @staticmethod
    def options(include_docs: bool = False):
        return _Option._options(Cutoff, include_docs=include_docs)


class Platform(_Option):
    """
    All of the supported platforms
    """

    AUTO = "auto", "Choose the platform automatically"
    CPU = "cpu", "Run on the CPU"
    CUDA = "cuda", "Run on the GPU using CUDA (nVidia)"
    OPENCL = "opencl", "Run on the GPU using OpenCL (all GPUs)"
    METAL = "metal", "Run on the GPU using Metal (Apple)"
    HIP = "hip", "Run on the GPU using HIP (AMD)"
    REFERENCE = "reference", "Run on CPU using the reference implementation"

    @staticmethod
    def create(option: str):
        return _Option._create(Platform, option)

    @staticmethod
    def options(include_docs: bool = False):
        return _Option._options(Platform, include_docs=include_docs)
