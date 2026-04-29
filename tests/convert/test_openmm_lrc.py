"""
Validate that the CustomVolumeForce-based LRC in the perturbable OpenMM system
gives the same total energy as a non-perturbable end-state system that uses the
standard NonbondedForce dispersion correction.

Tests use merged_ethane_methanol (merged_molecule.s3) with the 5 nearest waters
so they run quickly while still exercising a periodic cutoff system.
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
