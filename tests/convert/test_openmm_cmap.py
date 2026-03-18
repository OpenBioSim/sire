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
def test_openmm_cmap_perturbable(multichain_cmap, openmm_platform):
    """
    Verify that CMAPTorsionForce grids are correctly handled for a perturbable
    molecule. The pre-merged stream file merged_molecule_cmap.s3 contains a
    perturbable molecule whose two end states are identical (an identity
    perturbation of a CHARMM protein chain), so both end states carry the same
    CMAP backbone correction terms. The test checks that the perturbable code
    path correctly applies the same grids at all lambda values and that
    set_lambda does not corrupt the force parameters.
    """
    import openmm

    platform_name = openmm_platform or "CPU"

    mol0 = multichain_cmap[0]

    mols_pert = sr.load_test_files("merged_molecule_cmap.s3")
    mols_pert = sr.morph.link_to_reference(mols_pert)

    omm_map = {
        "constraint": "none",
        "cutoff": "none",
        "cutoff_type": "none",
        "platform": platform_name,
    }

    def get_cmap_torsion_grids(context):
        """
        Return list of (size, grid) for each CMAP torsion, dereferencing the
        map index.  Grid values are returned as plain floats (kJ/mol) so that
        pytest.approx can compare them.  This is map-count-agnostic: the
        non-perturbable path deduplicates maps while the perturbable path
        allocates one map per torsion, but the per-torsion grid values must
        agree.
        """
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

    def unique_grids(torsion_grids, decimals=3):
        """Return the sorted set of unique (size, rounded-grid) tuples.

        Torsion ordering can differ between the perturbable and non-perturbable
        code paths, so we compare the sets of unique grid shapes rather than
        comparing torsion-by-torsion."""
        seen = set()
        result = []
        for size, grid in torsion_grids:
            key = (size, tuple(round(v, decimals) for v in grid))
            if key not in seen:
                seen.add(key)
                result.append(key)
        return sorted(result)

    # Reference: non-perturbable molecule.
    mols_ref = sr.system.System()
    mols_ref.add(mol0)
    omm_ref = sr.convert.to(mols_ref, "openmm", map=omm_map)
    ref_torsion_grids = get_cmap_torsion_grids(omm_ref)

    assert len(ref_torsion_grids) > 0, "Reference context has no CMAP torsions"
    ref_unique = unique_grids(ref_torsion_grids)

    # Perturbable context — one context, lambda changed in place.
    omm_pert = sr.convert.to(mols_pert, "openmm", map=omm_map)

    # At both lambda=0 and lambda=1 the set of unique CMAP grids must match the
    # non-perturbable reference (cmap0 == cmap1 for an identity perturbation).
    # We compare sets of unique grids because the perturbable and non-perturbable
    # code paths may order torsions differently.
    for lam in (0.0, 1.0):
        omm_pert.set_lambda(lam)
        pert_torsion_grids = get_cmap_torsion_grids(omm_pert)

        assert len(pert_torsion_grids) == len(ref_torsion_grids), (
            f"Wrong number of CMAP torsions at lambda={lam}: "
            f"{len(pert_torsion_grids)} != {len(ref_torsion_grids)}"
        )

        pert_unique = unique_grids(pert_torsion_grids)
        assert pert_unique == ref_unique, (
            f"Set of unique CMAP grids differs from reference at lambda={lam}: "
            f"{len(pert_unique)} unique grids vs {len(ref_unique)} in reference"
        )
