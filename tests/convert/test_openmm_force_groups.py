"""
Tests for force-group assignment and energy caching in SOMMContext / LambdaLever.

Checks:
  - Named forces are assigned unique, valid OpenMM force group indices.
  - The cached energy matches a direct getState() evaluation at several lambda values.
  - Repeated get_potential_energy() calls without a state change return the same value.
  - The energy cache is cleared by setPositions(), setState(), setPeriodicBoxVectors(),
    and set_lambda().
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


pytestmark = pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)


def test_force_groups_assigned(perturbable_omm):
    """Every expected named force has a non-negative force group index."""
    lever = perturbable_omm.get_lambda_lever()

    force_names = lever.get_force_names()
    assert len(force_names) > 0, "No force names registered on LambdaLever"

    for name in force_names:
        grp = lever.get_force_group(name)
        assert grp >= 0, f"Force '{name}' has invalid group index {grp}"


def test_cached_energy_matches_full(perturbable_omm):
    """
    get_potential_energy() must match a direct getState() at lambda=0, 0.5, and 1.
    """
    import openmm

    omm = perturbable_omm

    for lam in (0.0, 0.5, 1.0):
        omm.set_lambda(lam)

        cached_kj = omm.get_potential_energy(to_sire_units=False).value_in_unit(
            openmm.unit.kilojoule_per_mole
        )
        full_kj = (
            omm.getState(getEnergy=True)
            .getPotentialEnergy()
            .value_in_unit(openmm.unit.kilojoule_per_mole)
        )

        assert cached_kj == pytest.approx(full_kj, abs=1e-3), (
            f"Cached energy {cached_kj:.6f} kJ/mol != full energy {full_kj:.6f} kJ/mol "
            f"at lambda={lam}"
        )


def test_cache_stable_without_state_change(perturbable_omm):
    """
    Two consecutive get_potential_energy() calls with no intervening state
    change must return identical values (second call served from cache).
    """
    omm = perturbable_omm
    omm.set_lambda(0.5)

    nrg1 = omm.get_potential_energy(to_sire_units=True).value()
    nrg2 = omm.get_potential_energy(to_sire_units=True).value()

    assert nrg1 == pytest.approx(
        nrg2, rel=1e-10
    ), "Energy changed between two consecutive calls with no state change"


def test_set_positions_invalidates_cache(perturbable_omm):
    """setPositions() must clear the energy cache."""
    omm = perturbable_omm
    omm.set_lambda(0.0)

    # Populate the cache.
    _ = omm.get_potential_energy(to_sire_units=False)
    assert "_total" in omm._energy_cache, "Cache should be populated after evaluation"

    # Set the same positions back — content unchanged but must still invalidate.
    positions = omm.getState(getPositions=True).getPositions()
    omm.setPositions(positions)

    assert len(omm._energy_cache) == 0, "Cache should be empty after setPositions()"


def test_set_state_invalidates_cache(perturbable_omm):
    """setState() must clear the energy cache."""
    omm = perturbable_omm
    omm.set_lambda(0.0)

    # Populate the cache.
    _ = omm.get_potential_energy(to_sire_units=False)
    assert "_total" in omm._energy_cache, "Cache should be populated after evaluation"

    # Round-trip through setState using the current state.
    state = omm.getState(getPositions=True)
    omm.setState(state)

    assert len(omm._energy_cache) == 0, "Cache should be empty after setState()"


def test_set_periodic_box_vectors_invalidates_cache(perturbable_omm):
    """setPeriodicBoxVectors() must clear the energy cache."""
    omm = perturbable_omm
    omm.set_lambda(0.0)

    # Populate the cache.
    _ = omm.get_potential_energy(to_sire_units=False)
    assert "_total" in omm._energy_cache, "Cache should be populated after evaluation"

    # Set the same box vectors back — must still invalidate.
    box = omm.getState(getPositions=True).getPeriodicBoxVectors()
    omm.setPeriodicBoxVectors(*box)

    assert (
        len(omm._energy_cache) == 0
    ), "Cache should be empty after setPeriodicBoxVectors()"


def test_set_lambda_invalidates_cache(perturbable_omm):
    """set_lambda() must clear the energy cache."""
    omm = perturbable_omm
    omm.set_lambda(0.0)

    # Populate the cache.
    _ = omm.get_potential_energy(to_sire_units=False)
    assert "_total" in omm._energy_cache, "Cache should be populated after evaluation"

    omm.set_lambda(0.5)

    assert len(omm._energy_cache) == 0, "Cache should be empty after set_lambda()"


# ---------------------------------------------------------------------------
# Tests for hasChanged() / wasForceChanged() — verifying that the C++ lambda
# lever correctly detects when morphed parameter values actually change vs
# when they are pinned and unchanged.
# ---------------------------------------------------------------------------

# Levers that control each named force.
_FORCE_LEVERS = {
    "bond": ["bond_k", "bond_length"],
    "angle": ["angle_k", "angle_size"],
    "torsion": ["torsion_k", "torsion_phase"],
    "clj": ["charge", "sigma", "epsilon", "alpha", "kappa", "charge_scale", "lj_scale"],
    # cmap_grid is the only lever for CMAPTorsionForce. By default it is coupled
    # to torsion_k, so without an explicit equation it would morph whenever
    # torsion_k does. _make_fixed_schedule sets an explicit equation (l.initial())
    # which breaks that coupling and pins CMAP independently.
    # Tested with a molecule that has perturbable CMAP terms (merged_molecule_cmap.s3).
    "cmap": ["cmap_grid"],
}

# Forces whose changed-state is tied to the "clj" levers.
_CLJ_RELATED = {"clj", "ghost/ghost", "ghost/non-ghost", "ghost-14"}


def _make_fixed_schedule(fixed_levers):
    """
    Return a single-stage LambdaSchedule that morphs all parameters linearly,
    except for the levers listed in *fixed_levers* which are pinned to initial.
    """
    l = sr.cas.LambdaSchedule()
    l.add_stage("morph", (1 - l.lam()) * l.initial() + l.lam() * l.final())
    for lever in fixed_levers:
        l.set_equation(stage="morph", lever=lever, equation=l.initial())
    return l


@pytest.mark.parametrize("fixed_force", list(_FORCE_LEVERS.keys()))
def test_fixed_lever_not_changed(merged_ethane_methanol, openmm_platform, fixed_force):
    """
    When all levers controlling a force are pinned to their initial values,
    wasForceChanged() must return False for that force after a lambda step.
    All other forces (whose levers still morph) must return True.
    """
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

    # Prime at lambda=0 so prev_cache is populated for the next step.
    omm.set_lambda(0.0)

    # Advance lambda — hasChanged() now compares against the lambda=0 values.
    omm.set_lambda(0.5)

    # The pinned force must NOT be marked changed.
    if fixed_force == "clj":
        for name in _CLJ_RELATED:
            if lever.get_force_group(name) >= 0:
                assert not lever.was_force_changed(
                    name
                ), f"'{name}' should not be changed when its levers are pinned"
    else:
        assert not lever.was_force_changed(
            fixed_force
        ), f"'{fixed_force}' should not be changed when its levers are pinned"

    # All other morphing forces must be marked changed.
    # Exclude "cmap": molecules without CMAP terms have no CMAP parameters to
    # change, so was_force_changed("cmap") is correctly False regardless of
    # pinning. The reverse direction (CMAP not changed when pinned) is covered
    # by the fixed_force="cmap" case which uses a CMAP molecule.
    other_forces = set(_FORCE_LEVERS.keys()) - {fixed_force, "cmap"}
    for other in other_forces:
        if other == "clj":
            if lever.get_force_group("clj") >= 0:
                assert lever.was_force_changed(
                    "clj"
                ), f"'clj' should be changed (fixed_force='{fixed_force}')"
        else:
            if lever.get_force_group(other) >= 0:
                assert lever.was_force_changed(
                    other
                ), f"'{other}' should be changed (fixed_force='{fixed_force}')"
