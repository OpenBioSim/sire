__all__ = [
    "sire_to_rdkit",
    "rdkit_to_sire",
    "rdkit_to_smiles",
    "rdkit_to_smarts",
    "rdkit_remove_hydrogens",
    "smiles_to_rdkit",
    "smarts_to_rdkit",
    "sire_to_openmm",
    "openmm_to_sire",
    "openmm_extract_coordinates",
    "openmm_extract_coordinates_and_velocities",
    "openmm_extract_space",
    "sire_to_gemmi",
    "gemmi_to_sire",
    "supported_formats",
    "LambdaLever",
    "PerturbableOpenMMMolecule",
    "OpenMMMetaData",
    "SOMMContext",
    "QMEngine",
    "PyQMCallback",
    "PyQMEngine",
]

try:
    from ._SireRDKit import (
        sire_to_rdkit,
        rdkit_to_sire,
        rdkit_to_smiles,
        rdkit_to_smarts,
        rdkit_remove_hydrogens,
        smiles_to_rdkit,
        smarts_to_rdkit,
        _register_smarts_search,
    )

    _has_rdkit = True
    _register_smarts_search()

except Exception as e:
    _rdkit_import_error = e

    # RDKit support is not available
    def _no_rdkit():
        print(_rdkit_import_error)
        raise ModuleNotFoundError(
            "Unable to convert to/from RDKit as it is not installed. "
            "Please install using `conda install -c conda-forge rdkit` "
            "and then re-run this script."
        )

    _has_rdkit = False

    def sire_to_rdkit(*args, **kwargs):
        _no_rdkit()

    def rdkit_to_sire(*args, **kwargs):
        _no_rdkit()

    def rdkit_to_smiles(*args, **kwargs):
        _no_rdkit()

    def rdkit_to_smarts(*args, **kwargs):
        _no_rdkit()

    def rdkit_remove_hydrogens(*args, **kwargs):
        _no_rdkit()

    def smiles_to_rdkit(*args, **kwargs):
        _no_rdkit()

    def smarts_to_rdkit(*args, **kwargs):
        _no_rdkit()


try:
    from ._SireOpenMM import sire_to_openmm_system as _sire_to_openmm_system
    from ._SireOpenMM import openmm_system_to_sire as _openmm_system_to_sire
    from ._SireOpenMM import extract_coordinates as _openmm_extract_coordinates
    from ._SireOpenMM import (
        extract_coordinates_and_velocities as _openmm_extract_coordinates_and_velocities,
    )
    from ._SireOpenMM import extract_space as _openmm_extract_space
    from ._SireOpenMM import minimise_openmm_context as _minimise_openmm_context

    from ._sommcontext import SOMMContext
    from ._perturbablemol import (
        _changed_atoms,
        _changed_bonds,
        _changed_angles,
        _changed_torsions,
        _changed_exceptions,
        _changed_constraints,
        _get_lever_values,
    )

    from ._SireOpenMM import (
        LambdaLever,
        PerturbableOpenMMMolecule,
        OpenMMMetaData,
        QMEngine,
        PyQMCallback,
        PyQMEngine,
    )

    try:
        from ._SireOpenMM import TorchQMEngine
        __all__.append("TorchQMEngine")
    except:
        pass

    from ..._pythonize import _pythonize

    _pythonize(
        [
            LambdaLever,
            PerturbableOpenMMMolecule,
            OpenMMMetaData,
            QMEngine,
            PyQMCallback,
            PyQMEngine,
        ],
        delete_old=True,
    )

    try:
        _pythonize(TorchQMEngine, delete_old=True)
    except:
        pass

    PerturbableOpenMMMolecule.changed_atoms = _changed_atoms
    PerturbableOpenMMMolecule.changed_bonds = _changed_bonds
    PerturbableOpenMMMolecule.changed_angles = _changed_angles
    PerturbableOpenMMMolecule.changed_torsions = _changed_torsions
    PerturbableOpenMMMolecule.changed_exceptions = _changed_exceptions
    PerturbableOpenMMMolecule.changed_constraints = _changed_constraints
    PerturbableOpenMMMolecule.get_lever_values = _get_lever_values

    _has_openmm = True

    def openmm_to_sire(mols, map):
        import openmm

        if type(mols) is not openmm.Context:
            raise TypeError(
                "You can only convert an openmm.Context to sire, not "
                f"a {type(mols)}."
            )

        # Need to be sure that 'mols' is an openmm.System or else
        # we will crash!
        system = mols.getSystem()

        if type(system) is not openmm.System:
            raise TypeError(
                "You can only convert an openmm.System to sire, not "
                f"a {type(system)}"
            )

        sire_mols = _openmm_system_to_sire(system, map)

        # now set the coordinates and velocities...
        # mols.getState()... state.getCoordinates()... state.getVelocities()

        return sire_mols

    def sire_to_openmm(mols, map):
        import openmm

        # OpenMM has system data spread over several objects.
        # The forces / parameters are in an openmm.System.
        # We create this first...
        system = openmm.System()

        # Next, we need to create an openmm.Integrator, so that we can
        # then create an openmm.Context, into which the coordinate
        # and velocity data can be placed

        from ...units import femtosecond, picosecond, kelvin, atm
        from ...move import Ensemble

        if map.specified("timestep"):
            timestep = map["timestep"].value()
        else:
            timestep = 1 * femtosecond

        if not timestep.has_same_units(femtosecond):
            raise TypeError(
                "The timestep should be in units of time. You cannot use "
                f"'{timestep}'"
            )

        timestep_in_fs = timestep.to(femtosecond)
        timestep = timestep.to(picosecond) * openmm.unit.picosecond
        map.set("timestep_in_fs", timestep_in_fs)

        ensemble = Ensemble(map=map)

        if map.specified("cutoff"):
            # we need to make sure that this is a unit
            cutoff = map["cutoff"]

            if cutoff.has_source():
                cutoff = cutoff.source()

                if cutoff.lower() == "none" or cutoff.lower().startswith("infinit"):
                    map.set("cutoff_type", "NONE")
                    map.unset("cutoff")
                elif cutoff.lower() == "auto":
                    map.unset("cutoff")
                elif cutoff != "cutoff":
                    from ... import u

                    map.set("cutoff", u(cutoff))

        if map.specified("integrator"):
            integrator = map["integrator"]

            if integrator.has_value():
                integrator = integrator.value()
            else:
                integrator = integrator.source()
        else:
            integrator = None

        if map.specified("friction"):
            friction = map["friction"].value()
        else:
            friction = 1.0 / picosecond

        friction = friction.to(1.0 / picosecond) / openmm.unit.picosecond

        use_andersen = False
        temperature = None

        if isinstance(integrator, str):
            from ...options import Integrator

            integrator = Integrator.create(integrator)

            if integrator == "verlet" or integrator == "leapfrog":
                if not ensemble.is_nve():
                    raise ValueError(
                        "You cannot use a verlet integrator with the " f"{ensemble}"
                    )

                integrator = openmm.VerletIntegrator(timestep)

            elif integrator != "auto":
                temperature = ensemble.temperature().to(kelvin) * openmm.unit.kelvin

                if ensemble.is_nve():
                    raise ValueError(
                        f"You cannot use a {integrator} integrator "
                        f"with the ensemble {ensemble}"
                    )

                if integrator == "langevin_middle":
                    integrator = openmm.LangevinMiddleIntegrator(
                        temperature, friction, timestep
                    )

                elif integrator == "langevin":
                    integrator = openmm.LangevinIntegrator(
                        temperature, friction, timestep
                    )

                elif integrator == "nose_hoover":
                    integrator = openmm.NoseHooverIntegrator(
                        temperature, friction, timestep
                    )

                elif integrator == "brownian":
                    integrator = openmm.BrownianIntegrator(
                        temperature, friction, timestep
                    )

                elif integrator == "andersen":
                    # use a verlet integrator and switch on the
                    # andersen thermostat with the specified frequency
                    integrator = openmm.VerletIntegrator(timestep)
                    use_andersen = True

                else:
                    raise ValueError(f"Unrecognised integrator {integrator}")

        if integrator is None:
            if ensemble.is_nve():
                integrator = openmm.VerletIntegrator(timestep)
            else:
                integrator = openmm.LangevinMiddleIntegrator(
                    ensemble.temperature().to(kelvin) * openmm.unit.kelvin,
                    friction,
                    timestep,
                )

                temperature = ensemble.temperature().to(kelvin) * openmm.unit.kelvin
        elif openmm.Integrator not in type(integrator).mro():
            raise TypeError(
                f"Cannot cast the integrator {integrator} to the correct "
                "type. It should be a string or an openmm.Integrator object"
            )

        if map.specified("constraint"):
            from ...options import Constraint

            constraint = Constraint.create(map.get_string("constraint"))

            if constraint == "auto":
                # choose the constraint based on the timestep
                if timestep_in_fs > 4:
                    # need constraint on everything
                    constraint = "bonds-not-heavy-perturbed"

                elif timestep_in_fs > 1:
                    # need it just on H bonds and angles
                    constraint = "h-bonds-not-heavy-perturbed"

                else:
                    # can get away with no constraints
                    constraint = "none"

            map.set("constraint", constraint)

        if map.specified("perturbable_constraint"):
            from ...options import PerturbableConstraint

            constraint = PerturbableConstraint.create(
                map.get_string("perturbable_constraint")
            )

            if constraint == "auto":
                # only apply the constraint to non-perturbed hydrogens
                constraint = "h-bonds-not-heavy-perturbed"

            map.set("perturbable_constraint", constraint)

        # Next, convert the sire system to an openmm system

        # system must be an openmm.System() or else we will crash!
        openmm_metadata = _sire_to_openmm_system(system, mols, map)

        # If we want temperature controlled by an Andersen thermostat
        # then add this here
        if use_andersen:
            system.addForce(openmm.AndersenThermostat(temperature, friction))

        # If we want NPT and this is periodic then we have to
        # add the barostat to the system
        if ensemble.is_npt():
            if not system.usesPeriodicBoundaryConditions():
                raise ValueError(
                    "You cannot run a constant pressure simulation "
                    "on a system with a non-periodic space."
                )

            barostat_freq = 25

            if map.specified("barostat_frequency"):
                barostat_freq = map["barostat_frequency"].value().as_integer()

            pressure = ensemble.pressure().to(atm) * openmm.unit.atmosphere

            barostat = openmm.MonteCarloBarostat(pressure, temperature, barostat_freq)

            system.addForce(barostat)

        platform = None

        if map.specified("platform"):
            from ...options import Platform

            desired_platform = Platform.create(map.get_string("platform"))

            # only look for the desired platform if it is not "auto"
            if desired_platform != "auto":
                platforms = []

                for i in range(0, openmm.Platform.getNumPlatforms()):
                    p = openmm.Platform.getPlatform(i)

                    if (p.getName().lower() == desired_platform.lower()) or (
                        p.getName() == "HIP" and desired_platform.lower() == "metal"
                    ):
                        platform = p
                        break
                    else:
                        platforms.append(p.getName().lower())

                if platform is None:
                    platforms = ", ".join(platforms)
                    raise ValueError(
                        f"Cannot create the openmm platform {desired_platform} "
                        "as this is not supported by this installation of "
                        f"openmm. Available platforms are [{platforms}]"
                    )

        if platform is None:
            # just find the fastest platform - this will be "metal" if that
            # is available and we are on Mac, or CUDA if CUDA works,
            # or OpenCL if OpenCL works, or CPU if nothing is left...
            import sys

            platforms = {}

            for i in range(0, openmm.Platform.getNumPlatforms()):
                p = openmm.Platform.getPlatform(i)
                platforms[p.getName().lower()] = p

            platform = None

            if sys.platform == "darwin":
                if "hip" in platforms:
                    platform = platforms["hip"]
                elif "metal" in platforms:
                    platform = platforms["metal"]

            if platform is None:
                if "cuda" in platforms:
                    platform = platforms["cuda"]

                elif "opencl" in platforms:
                    platform = platforms["opencl"]

                elif "cpu" in platforms:
                    platform = platforms["cpu"]

                elif len(platforms) > 0:
                    platform = platforms[list(platforms.keys())[0]]

            if platform is None:
                raise ValueError(
                    "This installation of openmm is broken, as there "
                    "are no available platforms!"
                )

        supported_properties = platform.getPropertyNames()

        if "Precision" in supported_properties and map.specified("precision"):
            precision = map.get_string("precision")
            platform.setPropertyDefaultValue("Precision", precision)

        if "Threads" in supported_properties and map.specified("threads"):
            try:
                threads = map["threads"].value().as_integer()
            except Exception:
                threads = map["threads"].source()

            platform.setPropertyDefaultValue("Threads", str(threads))

        if "DeviceIndex" in supported_properties and map.specified("device"):
            try:
                device_index = map["device"].value().as_integer()
            except Exception:
                device_index = map["device"].source()

            platform.setPropertyDefaultValue("DeviceIndex", str(device_index))

        if map.specified("cpu_pme") and "UseCpuPme" in supported_properties:
            usecpu = int(map["cpu_pme"].value().as_boolean())
            platform.setPropertyDefaultValue("UseCpuPme", str(usecpu).lower())

        try:
            from ._sommcontext import SOMMContext

            context = SOMMContext(
                system=system,
                integrator=integrator,
                platform=platform,
                metadata=openmm_metadata,
                map=map,
            )
        except Exception as e:
            raise ValueError(
                "There was a problem creating the OpenMM context. Perhaps "
                "the platform was not supported for this system, options "
                f"or on this computer? The error message is: {e}"
            )

        return context

    def openmm_extract_coordinates(state, mols, perturbable_maps=None, map=None):
        from ...base import create_map

        map = create_map(map)

        if perturbable_maps is None:
            perturbable_maps = {}

        return _openmm_extract_coordinates(
            state=state, mols=mols, perturbable_maps=perturbable_maps, map=map
        )

    def openmm_extract_coordinates_and_velocities(
        state, mols, perturbable_maps=None, map=None
    ):
        from ...base import create_map

        map = create_map(map)

        if perturbable_maps is None:
            from ..Mol import MolNum

            perturbable_maps = {}

        return _openmm_extract_coordinates_and_velocities(
            state=state, mols=mols, perturbable_maps=perturbable_maps, map=map
        )

    def openmm_extract_space(state):
        return _openmm_extract_space(state)

    def minimise_openmm_context(
        context,
        max_iterations: int = 10000,
        tolerance: float = 10.0,
        max_restarts: int = 10,
        max_ratchets: int = 20,
        ratchet_frequency: int = 500,
        starting_k: float = 100.0,
        ratchet_scale: float = 2.0,
        max_constraint_error: float = 0.01,
        timeout: str = "300s",
    ):
        return _minimise_openmm_context(
            context,
            max_iterations=max_iterations,
            tolerance=tolerance,
            max_restarts=max_restarts,
            max_ratchets=max_ratchets,
            ratchet_frequency=ratchet_frequency,
            starting_k=starting_k,
            ratchet_scale=ratchet_scale,
            max_constraint_error=max_constraint_error,
            timeout=timeout,
        )

except Exception as e:
    _openmm_import_exception = e

    # OpenMM support is not available
    def _no_openmm():
        print(_openmm_import_exception)

        raise ModuleNotFoundError(
            "Unable to convert to/from OpenMM as this code hasn't been "
            "written yet. We hope to support this soon!"
        )

    _has_openmm = False

    class LambdaLever:
        def __init__(self, *args, **kwargs):
            _no_openmm()

    class PerturbableOpenMMMolecule:
        def __init__(self, *args, **kwargs):
            _no_openmm()

    class OpenMMMetaData:
        def __init__(self, *args, **kwargs):
            _no_openmm()

    class SOMMContext:
        def __init__(self, *args, **kwargs):
            _no_openmm()

    def sire_to_openmm(*args, **kwargs):
        _no_openmm()

    def openmm_to_sire(*args, **kwargs):
        _no_openmm()

    def openmm_extract_coordinates(*arg, **kwargs):
        _no_openmm()

    def openmm_extract_coordinates_and_velocities(*args, **kwargs):
        _no_openmm()

    def openmm_extract_space(*args, **kwargs):
        _no_openmm()

    def minimise_openmm_context(*args, **kwargs):
        _no_openmm()


try:
    from ._SireGemmi import sire_to_gemmi, gemmi_to_sire, _register_pdbx_loader

    # make sure we have also import gemmi so that we
    # have the gemmi objects registered with python
    import gemmi as _gemmi  # noqa: F401

    _has_gemmi = True
    _register_pdbx_loader()
except Exception as e:
    _gemmi_import_error = e

    # Gemmi support is not available
    def _no_gemmi():
        print(_gemmi_import_error)
        raise ModuleNotFoundError(
            "Unable to convert to/from Gemmi as it is not installed. "
            "Please install using `conda install -c conda-forge gemmi` "
            "and then re-run this script."
        )

    _has_gemmi = False

    def sire_to_gemmi(*args, **kwargs):
        _no_gemmi()

    def gemmi_to_sire(*args, **kwargs):
        _no_gemmi()


def supported_formats():
    """Return all of the formats supported by this installation"""
    f = ["sire"]

    if _has_openmm:
        f.append("openmm")

    if _has_rdkit:
        f.append("rdkit")

    if _has_gemmi:
        f.append("gemmi")

    import sys

    # BioSimSpace needs to have already been loaded
    # otherwise it can't be imported
    if "BioSimSpace" in sys.modules:
        f.append("biosimspace")

    f.sort()

    return f
