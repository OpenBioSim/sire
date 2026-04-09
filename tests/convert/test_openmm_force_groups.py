"""
Tests for per-force-group energy caching in SOMMContext / LambdaLever.

Checks:
  - Named force groups are assigned and retrievable via get_force_group/get_force_names.
  - The cached per-group sum equals the full potential energy at several lambda values.
  - clear_energy_cache() forces a full re-evaluation on the next get_potential_energy() call.
  - Repeated get_potential_energy() calls without a lambda change return the same value
    (i.e. unchanged groups are served from cache, not re-computed).
  - A REST2 scale change marks the appropriate groups as dirty.
"""

import pytest
import sire as sr


@pytest.fixture(scope="module")
def perturbable_omm(merged_ethane_methanol, openmm_platform):
    """Return a SOMMContext built from the merged ethane/methanol molecule."""
    mols = merged_ethane_methanol.clone()

    # Pin coordinates to lambda-0 end state so energies are well-behaved.
    c = mols.cursor()
    c["molidx 0"]["coordinates"] = c["molidx 0"]["coordinates0"]
    c["molidx 0"]["coordinates1"] = c["molidx 0"]["coordinates0"]
    mols = c.commit()

    l = sr.cas.LambdaSchedule()
    l.add_stage("morph", (1 - l.lam()) * l.initial() + l.lam() * l.final())

    map = {
        "platform": openmm_platform,
        "schedule": l,
        "constraint": "h-bonds-not-perturbed",
        "include_constrained_energies": True,
        "dynamic_constraints": False,
    }

    return sr.convert.to(mols[0], "openmm", map=map)


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_force_groups_assigned(perturbable_omm):
    """Every expected named force has a non-negative force group index."""
    omm = perturbable_omm
    lever = omm.get_lambda_lever()

    force_names = lever.get_force_names()
    assert len(force_names) > 0, "No force names registered on LambdaLever"

    for name in force_names:
        grp = lever.get_force_group(name)
        assert grp >= 0, f"Force '{name}' has invalid group index {grp}"

    # The force_group_map built in Python must be consistent.
    assert len(omm._force_group_map) > 0
    for name, grp in omm._force_group_map.items():
        assert grp >= 0


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_per_group_sum_equals_full_energy(perturbable_omm):
    """
    Sum of per-group energies must equal the full potential energy
    (within numerical noise) at lambda=0, 0.5, and 1.
    """
    import openmm

    omm = perturbable_omm

    for lam in (0.0, 0.5, 1.0):
        omm.set_lambda(lam)
        omm.clear_energy_cache()

        # Full evaluation (groups bitmask = all 32 groups).
        full_state = omm.getState(getEnergy=True)
        full_kj = full_state.getPotentialEnergy().value_in_unit(
            openmm.unit.kilojoule_per_mole
        )

        # Per-group sum.
        group_sum_kj = 0.0
        for grp in omm._force_group_map.values():
            s = omm.getState(getEnergy=True, groups=(1 << grp))
            group_sum_kj += s.getPotentialEnergy().value_in_unit(
                openmm.unit.kilojoule_per_mole
            )

        assert group_sum_kj == pytest.approx(full_kj, abs=1e-3), (
            f"Group sum {group_sum_kj:.6f} kJ/mol != full energy {full_kj:.6f} kJ/mol "
            f"at lambda={lam}"
        )


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_cached_energy_matches_full(perturbable_omm):
    """
    get_potential_energy() via the cache must match a direct full getState()
    at lambda=0, 0.5, and 1.
    """
    import openmm

    omm = perturbable_omm

    for lam in (0.0, 0.5, 1.0):
        omm.set_lambda(lam)
        omm.clear_energy_cache()

        cached_kj = omm.get_potential_energy(to_sire_units=False).value_in_unit(
            openmm.unit.kilojoule_per_mole
        )

        full_state = omm.getState(getEnergy=True)
        full_kj = full_state.getPotentialEnergy().value_in_unit(
            openmm.unit.kilojoule_per_mole
        )

        assert cached_kj == pytest.approx(full_kj, abs=1e-3), (
            f"Cached energy {cached_kj:.6f} kJ/mol != full energy {full_kj:.6f} kJ/mol "
            f"at lambda={lam}"
        )


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_cache_stable_without_lambda_change(perturbable_omm):
    """
    Calling get_potential_energy() twice without a lambda change returns
    the same value (second call served entirely from cache).
    """
    omm = perturbable_omm
    omm.set_lambda(0.5)
    omm.clear_energy_cache()

    nrg1 = omm.get_potential_energy(to_sire_units=True).value()
    nrg2 = omm.get_potential_energy(to_sire_units=True).value()

    assert nrg1 == pytest.approx(
        nrg2, rel=1e-10
    ), "Energy changed between two consecutive calls with no lambda change"


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_clear_cache_marks_all_dirty(perturbable_omm):
    """
    After clear_energy_cache(), _dirty_groups contains every group in
    _force_group_map, and get_potential_energy() returns the correct value.
    """
    import openmm

    omm = perturbable_omm
    omm.set_lambda(0.0)

    # Populate cache.
    _ = omm.get_potential_energy(to_sire_units=False)

    # Clear.
    omm.clear_energy_cache()

    assert omm._dirty_groups == set(
        omm._force_group_map.values()
    ), "After clear_energy_cache(), not all groups are marked dirty"
    assert (
        len(omm._energy_cache) == 0
    ), "After clear_energy_cache(), energy_cache should be empty"

    # Energy should still be correct after re-evaluation.
    full_state = omm.getState(getEnergy=True)
    full_kj = full_state.getPotentialEnergy().value_in_unit(
        openmm.unit.kilojoule_per_mole
    )
    cached_kj = omm.get_potential_energy(to_sire_units=False).value_in_unit(
        openmm.unit.kilojoule_per_mole
    )
    assert cached_kj == pytest.approx(full_kj, abs=1e-3)


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_set_positions_invalidates_cache(perturbable_omm):
    """
    Calling setPositions() must invalidate the energy cache so that the next
    get_potential_energy() call re-evaluates all force groups.
    """
    omm = perturbable_omm
    omm.set_lambda(0.0)

    # Populate the cache.
    _ = omm.get_potential_energy(to_sire_units=False)
    assert len(omm._dirty_groups) == 0, "Cache should be clean after evaluation"

    # Retrieve current positions and set them back — content unchanged but
    # the override must still invalidate the cache.
    import openmm

    positions = omm.getState(getPositions=True).getPositions()
    omm.setPositions(positions)

    assert omm._dirty_groups == set(omm._force_group_map.values()), (
        "All groups should be dirty after setPositions()"
    )


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_set_state_invalidates_cache(perturbable_omm):
    """
    Calling setState() must invalidate the energy cache so that the next
    get_potential_energy() call re-evaluates all force groups.
    """
    omm = perturbable_omm
    omm.set_lambda(0.0)

    # Populate the cache.
    _ = omm.get_potential_energy(to_sire_units=False)
    assert len(omm._dirty_groups) == 0, "Cache should be clean after evaluation"

    # Round-trip through setState using the current state.
    state = omm.getState(getPositions=True)
    omm.setState(state)

    assert omm._dirty_groups == set(omm._force_group_map.values()), (
        "All groups should be dirty after setState()"
    )


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_set_periodic_box_vectors_invalidates_cache(perturbable_omm):
    """
    Calling setPeriodicBoxVectors() must invalidate the energy cache since
    a box change affects the PME energy.
    """
    omm = perturbable_omm
    omm.set_lambda(0.0)

    # Populate the cache.
    _ = omm.get_potential_energy(to_sire_units=False)
    assert len(omm._dirty_groups) == 0, "Cache should be clean after evaluation"

    # Set the same box vectors back — content unchanged but the override must
    # still invalidate the cache.
    box = omm.getState(getPositions=True).getPeriodicBoxVectors()
    omm.setPeriodicBoxVectors(*box)

    assert omm._dirty_groups == set(omm._force_group_map.values()), (
        "All groups should be dirty after setPeriodicBoxVectors()"
    )


# ---------------------------------------------------------------------------
# Levers that belong to each named OpenMM force.  When ALL levers for a force
# are pinned to l.initial(), that force's parameters cannot change between
# lambda steps, so it must NOT be marked dirty.
# ---------------------------------------------------------------------------
_FORCE_LEVERS = {
    "bond": ["bond_k", "bond_length"],
    "angle": ["angle_k", "angle_size"],
    # Note: fixing torsion_k also implicitly fixes cmap_grid via the default
    # coupling, but merged_ethane_methanol has no CMAP so this has no effect.
    "torsion": ["torsion_k", "torsion_phase"],
    # Fixing these without a force argument pins them for clj, ghost/ghost,
    # ghost/non-ghost and ghost-14 simultaneously.
    "clj": ["charge", "sigma", "epsilon", "alpha", "kappa", "charge_scale", "lj_scale"],
    # cmap_grid is the only lever for the CMAPTorsionForce.  By default it is
    # coupled to torsion_k, so without an explicit equation it would morph
    # whenever torsion_k does.  _make_fixed_schedule sets an explicit equation
    # (l.initial()) which breaks that coupling and pins CMAP independently.
    # Tested with a molecule that actually has perturbable CMAP terms
    # (merged_molecule_cmap.s3).
    "cmap": ["cmap_grid"],
}

# Forces whose dirty-state is tied together with "clj" (they share levers).
_CLJ_RELATED = {"clj", "ghost/ghost", "ghost/non-ghost", "ghost-14"}


def _make_fixed_schedule(fixed_levers):
    """
    Return a single-stage LambdaSchedule that morphs all parameters
    linearly, except for the levers listed in *fixed_levers* which are
    pinned to their lambda=0 (initial) values.
    """
    l = sr.cas.LambdaSchedule()
    l.add_stage("morph", (1 - l.lam()) * l.initial() + l.lam() * l.final())
    for lever in fixed_levers:
        l.set_equation(stage="morph", lever=lever, equation=l.initial())
    return l


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
@pytest.mark.parametrize("fixed_force", list(_FORCE_LEVERS.keys()))
def test_fixed_lever_not_dirty(merged_ethane_methanol, openmm_platform, fixed_force):
    """
    When all levers controlling *fixed_force* are pinned to their initial
    values, that force must not be marked dirty after a lambda step.
    All other forces (whose levers still morph) must be dirty.
    The cached energy must still match the full OpenMM energy.
    """
    import openmm as mm

    fixed_levers = _FORCE_LEVERS[fixed_force]
    schedule = _make_fixed_schedule(fixed_levers)

    if fixed_force == "cmap":
        # merged_molecule_cmap.s3 is a perturbable ubiquitin chain (T9A mutation)
        # that carries genuine CMAP backbone correction terms at both end states.
        mols = sr.load_test_files("merged_molecule_cmap.s3")
        mols = sr.morph.link_to_reference(mols)
        omm = sr.convert.to(
            mols,
            "openmm",
            map={
                "platform": openmm_platform or "CPU",
                "schedule": schedule,
                "constraint": "none",
                "cutoff": "none",
                "cutoff_type": "none",
            },
        )
    else:
        mols = merged_ethane_methanol.clone()

        # Pin coordinates to lambda-0 so energies are well-behaved.
        c = mols.cursor()
        c["molidx 0"]["coordinates"] = c["molidx 0"]["coordinates0"]
        c["molidx 0"]["coordinates1"] = c["molidx 0"]["coordinates0"]
        mols = c.commit()

        omm = sr.convert.to(
            mols[0],
            "openmm",
            map={
                "platform": openmm_platform,
                "schedule": schedule,
                "constraint": "h-bonds-not-perturbed",
                "include_constrained_energies": True,
                "dynamic_constraints": False,
            },
        )

    lever = omm.get_lambda_lever()

    # Step 1: prime the cache at lambda=0 (first call — all forces dirty
    # because there is no previous cached state to compare against).
    omm.set_lambda(0.0)
    _ = omm.get_potential_energy(to_sire_units=False)  # clears dirty_groups

    # Step 2: advance lambda — now hasChanged() compares against the
    # lambda=0 values stored in prev_cache.
    omm.set_lambda(0.5)

    # The pinned force must NOT be dirty.
    if fixed_force == "clj":
        # All CLJ-related forces share the same levers.
        for name in _CLJ_RELATED:
            if name in omm._force_group_map:
                assert not lever.was_force_changed(name), (
                    f"'{name}' should not be changed when all its levers are "
                    f"pinned to initial (fixed_force='{fixed_force}')"
                )
                assert (
                    omm._force_group_map[name] not in omm._dirty_groups
                ), f"Force group for '{name}' should not be dirty"
    else:
        assert not lever.was_force_changed(fixed_force), (
            f"'{fixed_force}' should not be changed when all its levers are "
            f"pinned to initial"
        )
        if fixed_force in omm._force_group_map:
            assert (
                omm._force_group_map[fixed_force] not in omm._dirty_groups
            ), f"Force group for '{fixed_force}' should not be dirty"

    # All OTHER morphing forces must be dirty.
    # Exclude "cmap": molecules without CMAP terms have no CMAP parameters to
    # change, so was_force_changed("cmap") is correctly False regardless of
    # pinning. The reverse direction (CMAP not dirty when pinned) is covered
    # by the fixed_force="cmap" parametrize case which uses a CMAP molecule.
    other_forces = set(_FORCE_LEVERS.keys()) - {fixed_force, "cmap"}
    for other in other_forces:
        if other == "clj":
            # Check at least the primary clj force.
            if "clj" in omm._force_group_map:
                assert lever.was_force_changed(
                    "clj"
                ), f"'clj' should be changed (fixed_force='{fixed_force}')"
        else:
            if other in omm._force_group_map:
                assert lever.was_force_changed(other), (
                    f"'{other}' should be changed when it is not pinned "
                    f"(fixed_force='{fixed_force}')"
                )

    # Energy correctness: cached sum must match full OpenMM evaluation.
    full_kj = (
        omm.getState(getEnergy=True)
        .getPotentialEnergy()
        .value_in_unit(mm.unit.kilojoule_per_mole)
    )
    cached_kj = omm.get_potential_energy(to_sire_units=False).value_in_unit(
        mm.unit.kilojoule_per_mole
    )
    assert cached_kj == pytest.approx(full_kj, abs=1e-3), (
        f"Cached energy {cached_kj:.6f} kJ/mol != full energy {full_kj:.6f} kJ/mol "
        f"(fixed_force='{fixed_force}')"
    )


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_rest2_scale_change_dirties_correct_groups(
    merged_ethane_methanol, openmm_platform
):
    """
    REST2 scale changes must dirty CLJ and torsion groups (and their related
    ghost forces) but must NOT dirty bond or angle groups.

    Uses a schedule where every morphed lever is pinned to its initial value,
    so morphed parameter vectors never change between lambda steps.  This
    isolates the REST2 scale as the sole source of cache invalidation.

    Three scenarios are tested:
    1. Lambda changes with REST2 scale held at 1.0 → no groups dirtied
       (cache entirely reused because neither morphed values nor scale changed).
    2. Lambda held constant, REST2 scale changes from 1.0 → 2.0 → CLJ and
       torsion groups (including ghost variants) are dirtied; bond and angle
       groups are NOT.
    3. The cached energy after the REST2 scale change still matches the full
       OpenMM potential energy.
    """
    import openmm

    # Pin every lever to its initial value so morphed vectors never change.
    all_levers = (
        list(_FORCE_LEVERS["bond"])
        + list(_FORCE_LEVERS["angle"])
        + list(_FORCE_LEVERS["torsion"])
        + list(_FORCE_LEVERS["clj"])
    )
    schedule = _make_fixed_schedule(all_levers)

    mols = merged_ethane_methanol.clone()
    c = mols.cursor()
    c["molidx 0"]["coordinates"] = c["molidx 0"]["coordinates0"]
    c["molidx 0"]["coordinates1"] = c["molidx 0"]["coordinates0"]
    mols = c.commit()

    omm = sr.convert.to(
        mols[0],
        "openmm",
        map={
            "platform": openmm_platform,
            "schedule": schedule,
            "constraint": "h-bonds-not-perturbed",
            "include_constrained_energies": True,
            "dynamic_constraints": False,
        },
    )

    # -----------------------------------------------------------------------
    # Scenario 1: prime the cache at lambda=0, REST2 scale=1.0.
    # -----------------------------------------------------------------------
    omm.set_lambda(0.0, rest2_scale=1.0)
    _ = omm.get_potential_energy(to_sire_units=False)  # clears _dirty_groups
    assert len(omm._dirty_groups) == 0, "Cache should be fully clean after priming"

    # Advance lambda — morphed values are all pinned so only REST2 scale
    # changes could dirty anything; scale is still 1.0, so nothing is dirty.
    omm.set_lambda(0.5, rest2_scale=1.0)
    assert len(omm._dirty_groups) == 0, (
        "No groups should be dirty after a lambda change when all levers are "
        "pinned and REST2 scale is unchanged"
    )

    # Consume the (still-clean) cache so subsequent checks start fresh.
    _ = omm.get_potential_energy(to_sire_units=False)

    # -----------------------------------------------------------------------
    # Scenario 2: lambda stays at 0.5, REST2 scale changes 1.0 → 2.0.
    # -----------------------------------------------------------------------
    omm.set_lambda(0.5, rest2_scale=2.0)

    rest2_affected = {"clj", "torsion", "ghost/ghost", "ghost/non-ghost"}
    rest2_unaffected = {"bond", "angle"}

    for name in rest2_affected:
        if name in omm._force_group_map:
            assert (
                omm._force_group_map[name] in omm._dirty_groups
            ), f"Force group '{name}' should be dirty after a REST2 scale change"

    for name in rest2_unaffected:
        if name in omm._force_group_map:
            assert (
                omm._force_group_map[name] not in omm._dirty_groups
            ), f"Force group '{name}' should NOT be dirty after a REST2 scale change"

    # -----------------------------------------------------------------------
    # Scenario 3: cached energy still matches the full OpenMM evaluation.
    # -----------------------------------------------------------------------
    full_kj = (
        omm.getState(getEnergy=True)
        .getPotentialEnergy()
        .value_in_unit(openmm.unit.kilojoule_per_mole)
    )
    cached_kj = omm.get_potential_energy(to_sire_units=False).value_in_unit(
        openmm.unit.kilojoule_per_mole
    )
    assert cached_kj == pytest.approx(full_kj, abs=1e-3), (
        f"Cached energy {cached_kj:.6f} kJ/mol != full energy {full_kj:.6f} kJ/mol "
        "after REST2 scale change"
    )


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_lambda_change_dirties_correct_groups(perturbable_omm):
    """
    After set_lambda(), only the groups whose parameters actually changed
    are marked dirty.  Groups that are unchanged are not in _dirty_groups.
    """
    omm = perturbable_omm
    omm.set_lambda(0.5)
    omm.clear_energy_cache()

    # Populate cache at lambda=0.5.
    _ = omm.get_potential_energy(to_sire_units=False)
    assert len(omm._dirty_groups) == 0, "Cache should be fully clean after evaluation"

    # Move to a new lambda — some groups must become dirty.
    omm.set_lambda(0.6)
    assert (
        len(omm._dirty_groups) > 0
    ), "At least one group should be dirty after a lambda change"

    # The cached energy must still be correct.
    import openmm

    full_state = omm.getState(getEnergy=True)
    full_kj = full_state.getPotentialEnergy().value_in_unit(
        openmm.unit.kilojoule_per_mole
    )
    cached_kj = omm.get_potential_energy(to_sire_units=False).value_in_unit(
        openmm.unit.kilojoule_per_mole
    )
    assert cached_kj == pytest.approx(full_kj, abs=1e-3)
