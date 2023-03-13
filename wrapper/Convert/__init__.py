__all__ = [
    "sire_to_rdkit",
    "rdkit_to_sire",
    "rdkit_to_smiles",
    "rdkit_remove_hydrogens",
    "smiles_to_rdkit",
    "sire_to_openmm",
    "openmm_to_sire",
    "openmm_extract_coordinates",
    "supported_formats",
]

try:
    from ._SireRDKit import (
        sire_to_rdkit,
        rdkit_to_sire,
        rdkit_to_smiles,
        rdkit_remove_hydrogens,
        smiles_to_rdkit,
    )

    _has_rdkit = True

except Exception:
    # RDKit support is not available
    def _no_rdkit():
        raise ModuleNotFoundError(
            "Unable to convert to/from RDKit as it is not installed. "
            "Please install using `mamba install -c conda-forge rdkit` "
            "and then re-run this script."
        )

    _has_rdkit = False

    def sire_to_rdkit(*args, **kwargs):
        _no_rdkit()

    def rdkit_to_sire(*args, **kwargs):
        _no_rdkit()

    def rdkit_to_smiles(*args, **kwargs):
        _no_rdkit()

    def rdkit_remove_hydrogens(*args, **kwargs):
        _no_rdkit()

    def smiles_to_rdkit(*args, **kwargs):
        _no_rdkit()


try:
    from ._SireOpenMM import (
        _sire_to_openmm_system,
        _openmm_system_to_sire,
        _set_openmm_coordinates_and_velocities,
        _openmm_extract_coordinates,
        _openmm_extract_coordinates_and_velocities,
    )

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

        from ...units import femtosecond, picosecond, kelvin
        from ...move import Ensemble

        if map.specified("timestep"):
            timestep = map["timestep"]
        else:
            timestep = 2 * femtosecond

        if not timestep.has_same_units(femtosecond):
            raise TypeError(
                "The timestep should be in units of time. You cannot use "
                f"'{timestep}'"
            )

        timestep = timestep.to(picosecond) * openmm.unit.picosecond

        ensemble = Ensemble(map=map)

        if map.specified("integrator"):
            integrator = map["integrator"]
        else:
            integrator = None

        if map.specified("friction"):
            friction = map["friction"]
        else:
            friction = 1.0 / picosecond

        friction = friction.to(1.0 / picosecond) / openmm.unit.picosecond

        if integrator is None:
            if ensemble.is_nve():
                integrator = openmm.VerletIntegrator(timestep)
            else:
                integrator = openmm.LangevinMiddleIntegrator(
                    ensemble.temperature().to(kelvin) * openmm.unit.kelvin,
                    friction,
                    timestep,
                )

        elif type(integrator) is str:
            integrator = integrator.lower()

            if integrator == "verlet":
                if not ensemble.is_nve():
                    raise ValueError(
                        "You cannot use a verlet integrator with the "
                        f"ensemble {ensemble}"
                    )
                integrator = openmm.VerletIntegrator(timestep)

            else:
                temperature = (
                    ensemble.temperature().to(kelvin) * openmm.unit.kelvin
                )

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

                else:
                    raise ValueError(f"Unrecognised integrator {integrator}")

        elif openmm.Integrator not in type(integrator).mro():
            raise TypeError(
                f"Cannot cast the integrator {integrator} to the correct "
                "type. It should be a string or an openmm.Integrator object"
            )

        # Next, convert the sire system to an openmm system

        # system must be an openmm.System() or else we will crash!
        coords_and_vels = _sire_to_openmm_system(system, mols, map)

        context = openmm.Context(system, integrator)

        # place the coordinates and velocities into the context
        _set_openmm_coordinates_and_velocities(context, coords_and_vels)

        # Then there is the OpenMM Topology and Simulation, that
        # hold the above data together with some metadata about the
        # system being simulated. I will need to think about how
        # I should handle this...

        return context

    def openmm_extract_coordinates(state, mols, map):
        return _openmm_extract_coordinates(state, mols, map)

    def openmm_extract_coordinates_and_velocities(state, mols, map):
        return _openmm_extract_coordinates_and_velocities(state, mols, map)

except Exception:
    # OpenMM support is not available
    def _no_openmm():
        raise ModuleNotFoundError(
            "Unable to convert to/from OpenMM as this code hasn't been "
            "written yet. We hope to support this soon!"
        )

    _has_openmm = False

    def sire_to_openmm(*args, **kwargs):
        _no_openmm()

    def openmm_to_sire(*args, **kwargs):
        _no_openmm()

    def openmm_extract_coordinates(state, mols, map):
        _no_openmm()

    def openmm_extract_coordinates_and_velocities(state, mols, map):
        _no_openmm()


def supported_formats():
    """Return all of the formats supported by this installation"""
    f = ["sire"]

    if _has_openmm:
        f.append("openmm")

    if _has_rdkit:
        f.append("rdkit")

    import sys

    # BioSimSpace needs to have already been loaded
    # otherwise it can't be imported
    if "BioSimSpace" in sys.modules:
        f.append("biosimspace")

    f.sort()

    return f
