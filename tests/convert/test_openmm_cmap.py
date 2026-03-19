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
    mols.add(multichain_cmap[2])

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

    # Create an OpenMM context via Sire's conversion layer, then get the
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


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_openmm_cmap_perturbable(openmm_platform):
    """
    Verify that CMAPTorsionForce grids are correctly updated by the lambda lever
    for a perturbable molecule with a genuine single-residue mutation.

    The pre-merged stream file merged_molecule_cmap.s3 contains a perturbable
    ubiquitin chain with a T9A mutation. cmap0 and cmap1 differ for exactly
    one torsion (the one centred on the mutated residue). The test checks that:

    - CMAP torsions are present at both end states.
    - At least one torsion grid differs between lambda=0 and lambda=1.
    - Most torsion grids are unchanged (only the mutated residue is affected).
    """
    import openmm
    import openmm.unit

    platform_name = openmm_platform or "CPU"

    mols_pert = sr.load_test_files("merged_molecule_cmap.s3")
    mols_pert = sr.morph.link_to_reference(mols_pert)

    omm_map = {
        "constraint": "none",
        "cutoff": "none",
        "cutoff_type": "none",
        "platform": platform_name,
    }

    def get_cmap_torsion_grids(context):
        """Return list of (size, grid) for each CMAP torsion, dereferencing
        the map index. Grid values are plain floats (kJ/mol)."""
        system = context.getSystem()
        for force in system.getForces():
            if isinstance(force, openmm.CMAPTorsionForce):
                maps = []
                for i in range(force.getNumMaps()):
                    size, grid = force.getMapParameters(i)
                    grid_floats = [
                        v.value_in_unit(openmm.unit.kilojoules_per_mole) for v in grid
                    ]
                    maps.append((size, grid_floats))
                result = []
                for t in range(force.getNumTorsions()):
                    map_idx = force.getTorsionParameters(t)[0]
                    result.append(maps[map_idx])
                return result
        return []

    omm_pert = sr.convert.to(mols_pert, "openmm", map=omm_map)

    omm_pert.set_lambda(0.0)
    grids_lam0 = get_cmap_torsion_grids(omm_pert)

    omm_pert.set_lambda(1.0)
    grids_lam1 = get_cmap_torsion_grids(omm_pert)

    assert len(grids_lam0) > 0, "No CMAP torsions at lambda=0"
    assert len(grids_lam1) == len(
        grids_lam0
    ), f"Torsion count differs between end states: {len(grids_lam0)} vs {len(grids_lam1)}"

    differing = sum(
        1
        for (s0, g0), (s1, g1) in zip(grids_lam0, grids_lam1)
        if s0 != s1 or any(round(a, 3) != round(b, 3) for a, b in zip(g0, g1))
    )

    assert differing > 0, (
        "Expected at least one torsion grid to differ between lambda=0 and lambda=1 "
        "for a genuine single-residue mutation"
    )
    assert differing < len(
        grids_lam0
    ), "Expected most torsion grids to be unchanged for a single-residue mutation"

    # Verify that changed_cmaps() correctly identifies the same set of
    # differing torsions as the direct grid comparison above.
    p_omm = mols_pert.molecule(0).perturbation().to_openmm(map=omm_map)
    changed = p_omm.changed_cmaps()
    assert (
        len(changed) == differing
    ), f"changed_cmaps() returned {len(changed)} torsions but expected {differing}"
