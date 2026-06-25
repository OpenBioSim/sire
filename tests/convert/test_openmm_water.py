import math

import pytest
import sire as sr

_skip_no_openmm = pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)


@pytest.fixture(scope="module")
def tip3p_mols():
    return sr.load_test_files("tip3p.s3")


@pytest.fixture(scope="module")
def tip4p_mols():
    return sr.load_test_files("tip4p.s3")


@pytest.fixture(scope="module")
def tip5p_mols():
    return sr.load_test_files("tip5p.s3")


@pytest.fixture(scope="module")
def opc_mols():
    return sr.load_test_files("opc.s3")


def _get_vsite_info(omm_context):
    """Return (n_vsites, n_zero_mass, n_bad_constraints) for an OpenMM context."""
    import openmm

    sys = omm_context.getSystem()
    natoms = sys.getNumParticles()

    zero_mass = {
        i
        for i in range(natoms)
        if sys.getParticleMass(i).value_in_unit(openmm.unit.dalton) == 0.0
    }
    n_vsites = sum(1 for i in zero_mass if sys.isVirtualSite(i))

    bad_constraints = sum(
        1
        for k in range(sys.getNumConstraints())
        if sys.getConstraintParameters(k)[0] in zero_mass
        or sys.getConstraintParameters(k)[1] in zero_mass
    )

    return n_vsites, len(zero_mass), bad_constraints


def _potential_energy(omm_context):
    import openmm

    return (
        omm_context.getState(getEnergy=True)
        .getPotentialEnergy()
        .value_in_unit(openmm.unit.kilojoules_per_mole)
    )


@_skip_no_openmm
def test_tip3p_no_virtual_sites(tip3p_mols, openmm_platform):
    """TIP3P has no EP atoms — no virtual sites should be registered."""
    omm_map = {
        "constraint": "h-bonds",
        "cutoff_type": "none",
        "cutoff": "none",
        "platform": openmm_platform or "CPU",
    }

    omm = sr.convert.to(tip3p_mols, "openmm", map=omm_map)
    n_vsites, n_zero_mass, bad_constraints = _get_vsite_info(omm)

    assert (
        n_zero_mass == 0
    ), f"TIP3P should have no zero-mass particles, got {n_zero_mass}"
    assert n_vsites == 0, f"TIP3P should have no virtual sites, got {n_vsites}"
    assert bad_constraints == 0

    e = _potential_energy(omm)
    assert not math.isnan(e), "Potential energy is NaN"


@_skip_no_openmm
def test_tip4p_virtual_sites(tip4p_mols, openmm_platform):
    """TIP4P: one ThreeParticleAverageSite per water, no bad constraints, finite energy."""
    import openmm

    omm_map = {
        "constraint": "h-bonds",
        "cutoff_type": "none",
        "cutoff": "none",
        "platform": openmm_platform or "CPU",
    }

    omm = sr.convert.to(tip4p_mols, "openmm", map=omm_map)
    n_vsites, n_zero_mass, bad_constraints = _get_vsite_info(omm)

    n_waters = tip4p_mols.num_molecules()

    assert (
        n_zero_mass == n_waters
    ), f"Expected {n_waters} zero-mass EP atoms, got {n_zero_mass}"
    assert n_vsites == n_waters, f"Expected {n_waters} virtual sites, got {n_vsites}"
    assert bad_constraints == 0, f"Constraints on virtual sites: {bad_constraints}"

    # All virtual sites should be ThreeParticleAverageSite
    sys = omm.getSystem()
    natoms = sys.getNumParticles()
    for i in range(natoms):
        if sys.isVirtualSite(i):
            vs = sys.getVirtualSite(i)
            assert isinstance(
                vs, openmm.ThreeParticleAverageSite
            ), f"Particle {i}: expected ThreeParticleAverageSite, got {type(vs).__name__}"

    e = _potential_energy(omm)
    assert not math.isnan(e), "Potential energy is NaN"


@_skip_no_openmm
def test_opc_virtual_sites(opc_mols, openmm_platform):
    """OPC: one ThreeParticleAverageSite per water, no bad constraints, finite energy."""
    import openmm

    omm_map = {
        "constraint": "h-bonds",
        "cutoff_type": "none",
        "cutoff": "none",
        "platform": openmm_platform or "CPU",
    }

    omm = sr.convert.to(opc_mols, "openmm", map=omm_map)
    n_vsites, n_zero_mass, bad_constraints = _get_vsite_info(omm)

    n_waters = opc_mols.num_molecules()

    assert (
        n_zero_mass == n_waters
    ), f"Expected {n_waters} zero-mass EP atoms, got {n_zero_mass}"
    assert n_vsites == n_waters, f"Expected {n_waters} virtual sites, got {n_vsites}"
    assert bad_constraints == 0, f"Constraints on virtual sites: {bad_constraints}"

    sys = omm.getSystem()
    natoms = sys.getNumParticles()
    for i in range(natoms):
        if sys.isVirtualSite(i):
            vs = sys.getVirtualSite(i)
            assert isinstance(
                vs, openmm.ThreeParticleAverageSite
            ), f"Particle {i}: expected ThreeParticleAverageSite, got {type(vs).__name__}"

    e = _potential_energy(omm)
    assert not math.isnan(e), "Potential energy is NaN"


@_skip_no_openmm
def test_tip5p_virtual_sites(tip5p_mols, openmm_platform):
    """TIP5P: two OutOfPlaneSite virtual sites per water, no bad constraints, finite energy."""
    import openmm

    omm_map = {
        "constraint": "h-bonds",
        "cutoff_type": "none",
        "cutoff": "none",
        "platform": openmm_platform or "CPU",
    }

    omm = sr.convert.to(tip5p_mols, "openmm", map=omm_map)
    n_vsites, n_zero_mass, bad_constraints = _get_vsite_info(omm)

    n_waters = tip5p_mols.num_molecules()

    assert (
        n_zero_mass == 2 * n_waters
    ), f"Expected {2 * n_waters} zero-mass EP atoms, got {n_zero_mass}"
    assert (
        n_vsites == 2 * n_waters
    ), f"Expected {2 * n_waters} virtual sites, got {n_vsites}"
    assert bad_constraints == 0, f"Constraints on virtual sites: {bad_constraints}"

    sys = omm.getSystem()
    natoms = sys.getNumParticles()
    for i in range(natoms):
        if sys.isVirtualSite(i):
            vs = sys.getVirtualSite(i)
            assert isinstance(
                vs, openmm.OutOfPlaneSite
            ), f"Particle {i}: expected OutOfPlaneSite, got {type(vs).__name__}"

    e = _potential_energy(omm)
    assert not math.isnan(e), "Potential energy is NaN"
