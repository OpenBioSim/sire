"""
Validate that the CustomVolumeForce-based LRC in the perturbable OpenMM system
gives the same total energy as a non-perturbable end-state system that uses the
standard NonbondedForce dispersion correction.

Tests use merged_ethane_methanol (merged_molecule.s3) with the 5 nearest waters
so they run quickly while still exercising a periodic cutoff system.

A separate test exercises the GCMC water LRC (GCMCLRCForce) using the ala_mols
system (alanine-dipeptide in water). GCMC is faked by decrementing n_w by hand
and verifying the energy change against the analytical formula.
"""

import pytest
import sire as sr


def _get_end_state(mol, state, remove_state):
    c = mol.cursor()
    for key in c.keys():
        if key.endswith(state):
            c[key.removesuffix(state)] = c[key]
            del c[key]
        elif key.endswith(remove_state):
            del c[key]
    c["is_perturbable"] = False
    return c.commit()


def _build_systems(mols, platform):
    space = mols.space()

    c = mols.cursor()
    c["molidx 0"]["coordinates"] = c["molidx 0"]["coordinates0"]
    c["molidx 0"]["coordinates1"] = c["molidx 0"]["coordinates0"]
    mols = c.commit()

    merge = mols[0]
    water = mols["closest 5 waters to molidx 0"]

    mols_pert = merge + water
    mols_ref = _get_end_state(merge, "0", "1") + water
    mols_pert_end = _get_end_state(merge, "1", "0") + water

    l = sr.cas.LambdaSchedule()
    l.add_stage("morph", (1 - l.lam()) * l.initial() + l.lam() * l.final())

    map = {
        "platform": platform,
        "schedule": l,
        "constraint": "h-bonds-not-perturbed",
        "include_constrained_energies": True,
        "dynamic_constraints": False,
        "use_dispersion_correction": True,
        "space": space,
    }

    omm_pert = sr.convert.to(mols_pert, "openmm", map=map)
    omm_ref = sr.convert.to(mols_ref, "openmm", map=map)
    omm_pert_end = sr.convert.to(mols_pert_end, "openmm", map=map)

    return omm_pert, omm_ref, omm_pert_end


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_lrc_lambda0_matches_reference_end_state(
    merged_ethane_methanol, openmm_platform
):
    """
    Perturbable system at lambda=0 must give the same energy as a freshly built
    non-perturbable reference (lambda=0) end state with standard NonbondedForce LRC.
    """
    omm_pert, omm_ref, _ = _build_systems(
        merged_ethane_methanol.clone(), openmm_platform
    )

    omm_pert.set_lambda(0.0)
    nrg_pert = omm_pert.get_energy().value()
    nrg_ref = omm_ref.get_energy().value()

    assert nrg_pert == pytest.approx(nrg_ref, rel=1e-3)


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_lrc_lambda1_matches_perturbed_end_state(
    merged_ethane_methanol, openmm_platform
):
    """
    Perturbable system at lambda=1 must give the same energy as a freshly built
    non-perturbable perturbed (lambda=1) end state with standard NonbondedForce LRC.
    """
    omm_pert, _, omm_pert_end = _build_systems(
        merged_ethane_methanol.clone(), openmm_platform
    )

    omm_pert.set_lambda(1.0)
    nrg_pert = omm_pert.get_energy().value()
    nrg_pert_end = omm_pert_end.get_energy().value()

    assert nrg_pert == pytest.approx(nrg_pert_end, rel=1e-3)


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_gcmc_lrc(ala_mols, openmm_platform):
    """
    Test the GCMCLRCForce by faking a GCMC move: decrement n_w by 1 and verify
    the energy change matches the analytical formula.

    E(n_w) = (n_w * lrc_w_solute + n_w * (n_w - 1) * lrc_ww_half) / V

    dE = E(n_w - 1) - E(n_w) = -(lrc_w_solute + 2 * (n_w - 1) * lrc_ww_half) / V
    """
    # Build with GCMC LRC enabled, reserving 1 buffer water so n_w starts at
    # n_waters - 1.
    d = ala_mols.dynamics(
        platform=openmm_platform,
        map={
            "use_dispersion_correction": True,
            "use_gcmc_lrc": True,
            "num_gcmc_waters": 1,
        },
    )
    d.set_lambda(0.0)

    context = d.context()
    omm_system = context.getSystem()

    # Locate the GCMCLRCForce and its force group.
    gcmc_force = None
    for force in omm_system.getForces():
        if force.getName() == "GCMCLRCForce":
            gcmc_force = force
            break
    assert gcmc_force is not None, "GCMCLRCForce not found in system"
    fg = gcmc_force.getForceGroup()

    # Read the pre-computed LRC coefficients and initial n_w from the context.
    state = context.getState(getParameters=True, getPositions=True)
    params = state.getParameters()
    n_w = params["n_w"]
    lrc_w_solute = params["lrc_w_solute"]
    lrc_ww_half = params["lrc_ww_half"]

    # The LJ tail is attractive so both coefficients must be negative.
    assert lrc_w_solute < 0
    assert lrc_ww_half < 0
    assert n_w > 0

    # Get the GCMC LRC energy at the initial n_w.
    nrg1 = (
        context.getState(getEnergy=True, groups=(1 << fg)).getPotentialEnergy()._value
    )

    # "Turn off" one water by decrementing n_w.
    context.setParameter("n_w", n_w - 1)
    nrg2 = (
        context.getState(getEnergy=True, groups=(1 << fg)).getPotentialEnergy()._value
    )

    # Compute the expected energy change analytically.
    # Box vectors are in nm; volume in nm^3 matches the lrc_coeff units (kJ/mol·nm^3).
    box = state.getPeriodicBoxVectors()
    V = box[0][0]._value * box[1][1]._value * box[2][2]._value

    expected_delta = -(lrc_w_solute + 2 * (n_w - 1) * lrc_ww_half) / V
    assert nrg2 - nrg1 == pytest.approx(expected_delta, rel=1e-5)
