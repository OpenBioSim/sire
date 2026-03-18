import sire as sr
import pytest


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_openmm_cmap_energy(tmpdir, multichain_cmap, openmm_platform):
    """
    Verify that Sire correctly adds CMAPTorsionForce to the OpenMM context by
    comparing the total potential energy against a context built directly via
    the OpenMM Python API from the same input files.

    The multichain_cmap fixture is a periodic solvated system with three protein
    chains, each carrying CMAP backbone correction terms. Using it exercises the
    multi-molecule CMAP code path in the conversion layer.
    """
    import openmm
    import openmm.app
    import openmm.unit

    mols = sr.system.System()
    mols.add(multichain_cmap[0])
    mols.add(multichain_cmap[1])

    # Sanity-check: at least two molecules must carry CMAP so that the
    # multi-chain code path is exercised.
    cmap_mol_count = sum(1 for mol in mols.molecules() if mol.has_property("cmap"))
    assert (
        cmap_mol_count >= 2
    ), "Expected at least two molecules with CMAP terms in multichain_cmap"

    # Save the Sire system to AMBER files so the direct OpenMM path reads the
    # same topology and coordinates that Sire uses internally.
    dir_path = str(tmpdir.mkdir("cmap_omm"))
    prm7 = str(sr.save(mols, f"{dir_path}/system.prm7")[0])
    rst7 = str(sr.save(mols, f"{dir_path}/system.rst7")[0])

    platform_name = openmm_platform or "CPU"

    # Create and OpenMM context via Sire's conversion layer, then get the
    # potential energy.
    sire_map = {
        "constraint": "none",
        "cutoff": "none",
        "cutoff_type": "none",
        "platform": platform_name,
    }
    omm_sire = sr.convert.to(mols, "openmm", map=sire_map)
    sire_energy = (
        omm_sire.getState(getEnergy=True)
        .getPotentialEnergy()
        .value_in_unit(openmm.unit.kilojoules_per_mole)
    )

    # Create an OpenMM context directly from the AMBER files and get the
    # potential energy.
    prmtop = openmm.app.AmberPrmtopFile(prm7)
    inpcrd = openmm.app.AmberInpcrdFile(rst7)

    omm_system = prmtop.createSystem(
        nonbondedMethod=openmm.app.NoCutoff,
        constraints=None,
    )

    integrator = openmm.VerletIntegrator(0.001)
    platform = openmm.Platform.getPlatformByName(platform_name)
    omm_context = openmm.Context(omm_system, integrator, platform)
    omm_context.setPositions(inpcrd.positions)
    if inpcrd.boxVectors is not None:
        omm_context.setPeriodicBoxVectors(*inpcrd.boxVectors)

    direct_energy = (
        omm_context.getState(getEnergy=True)
        .getPotentialEnergy()
        .value_in_unit(openmm.unit.kilojoules_per_mole)
    )

    # Energies should agree to within 1 kJ/mol.
    assert sire_energy == pytest.approx(direct_energy, abs=1.0)
