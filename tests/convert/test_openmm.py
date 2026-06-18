import sire as sr
import pytest


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_openmm_single_energy_neopentane(neopentane_methane, openmm_platform):
    mol = neopentane_methane[0]

    # this function will extract the lambda0 or lambda1 end state
    mol0 = sr.morph.extract_reference(mol)
    mol1 = sr.morph.extract_perturbed(mol)

    map = {
        "space": sr.vol.Cartesian(),
        "platform": openmm_platform,
        "constraint": "none",
        "ignore_perturbations": True,
    }

    omm0 = sr.convert.to(mol0, "openmm", map=map)

    state0 = omm0.getState(getEnergy=True)

    energy0 = state0.getPotentialEnergy()

    energy0 = energy0.value_in_unit(energy0.unit)

    assert mol0.energy(map=map).to(sr.units.kJ_per_mol) == pytest.approx(
        energy0, abs=1e-3
    )

    omm1 = sr.convert.to(mol1, "openmm", map=map)

    state1 = omm1.getState(getEnergy=True)

    energy1 = state1.getPotentialEnergy()

    energy1 = energy1.value_in_unit(energy1.unit)

    assert mol1.energy(map=map).to(sr.units.kJ_per_mol) == pytest.approx(
        energy1, abs=1e-3
    )

    assert mol0.dynamics(map=map, cutoff="none").current_potential_energy().to(
        sr.units.kJ_per_mol
    ) == pytest.approx(energy0, abs=1e-3)

    assert mol1.dynamics(map=map, cutoff="none").current_potential_energy().to(
        sr.units.kJ_per_mol
    ) == pytest.approx(energy1, abs=1e-3)


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_openmm_single_energy(kigaki_mols, openmm_platform):
    mols = kigaki_mols

    mol = mols[0]

    map = {
        "space": sr.vol.Cartesian(),
        "platform": openmm_platform,
        "constraint": "bonds-h-angles",
    }

    omm = sr.convert.to(mol, "openmm", map=map)

    state = omm.getState(getPositions=True, getEnergy=True)

    positions = state.getPositions()

    assert len(positions) == mol.num_atoms()

    for i, atom in enumerate(mol.atoms()):
        c = atom.coords()
        p = positions[i]

        assert c.x().to(sr.units.nanometer) == pytest.approx(p.x)
        assert c.y().to(sr.units.nanometer) == pytest.approx(p.y)
        assert c.z().to(sr.units.nanometer) == pytest.approx(p.z)

    energy = state.getPotentialEnergy()

    # get this as a float in kJ mol-1
    energy = energy.value_in_unit(energy.unit)

    # these won't be exactly the same - this is 5227 +/- 0.1 kJ mol-1
    assert mol.energy(map=map).to(sr.units.kJ_per_mol) == pytest.approx(energy, abs=0.1)


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_openmm_multi_energy_small_cart(kigaki_mols, openmm_platform):
    # first, try just 50 molecules in a cartesian space
    mols = kigaki_mols[0:50]

    map = {
        "space": sr.vol.Cartesian(),
        "platform": openmm_platform,
        "constraint": "bonds-h-angles",
    }

    omm = sr.convert.to(mols, "openmm", map=map)

    state = omm.getState(getEnergy=True)

    energy = state.getPotentialEnergy()

    # get this as a float in kJ mol-1
    energy = energy.value_in_unit(energy.unit)

    # these won't be exactly the same - this is 4865.82 +/- 0.04
    assert mols.energy(map=map).to(sr.units.kJ_per_mol) == pytest.approx(
        energy, abs=0.5
    )


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
@pytest.mark.slow
def test_openmm_multi_energy_all_cart(kigaki_mols, openmm_platform):
    # use all of the molecules
    mols = kigaki_mols

    map = {
        "space": sr.vol.Cartesian(),
        "cutoff": 10000 * sr.units.angstrom,
        "cutoff_type": "REACTION_FIELD",
        "dielectric": 1.0,
        "platform": openmm_platform,
        "constraint": "bonds-h-angles",
    }

    omm = sr.convert.to(mols, "openmm", map=map)

    state = omm.getState(getEnergy=True)

    energy = state.getPotentialEnergy()

    # get this as a float in kJ mol-1
    energy = energy.value_in_unit(energy.unit)

    # -127882 +/- 0.5
    assert mols.energy(map=map).to(sr.units.kJ_per_mol) == pytest.approx(
        energy, abs=1.0
    )


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_openmm_multi_energy_all_cart_cutoff(kigaki_mols, openmm_platform):
    # use all of the molecules
    mols = kigaki_mols

    map = {
        "space": sr.vol.Cartesian(),
        "cutoff": 10 * sr.units.angstrom,
        "cutoff_type": "REACTION_FIELD",
        "dielectric": 78.0,
        "platform": openmm_platform,
        "constraint": "bonds-h-angles",
    }

    omm = sr.convert.to(mols, "openmm", map=map)

    state = omm.getState(getEnergy=True)

    energy = state.getPotentialEnergy()

    # get this as a float in kJ mol-1
    energy = energy.value_in_unit(energy.unit)

    # -125622.2 +/- 0.05
    assert mols.energy(map=map).to(sr.units.kJ_per_mol) == pytest.approx(
        energy, abs=0.5
    )


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_openmm_multi_energy_all_periodic_cutoff(kigaki_mols, openmm_platform):
    # use all of the molecules
    mols = kigaki_mols

    map = {
        "cutoff": 10 * sr.units.angstrom,
        "cutoff_type": "REACTION_FIELD",
        "dielectric": 78.0,
        "platform": openmm_platform,
        "constraint": "bonds-h-angles",
    }

    omm = sr.convert.to(mols, "openmm", map=map)

    state = omm.getState(getEnergy=True, enforcePeriodicBox=True)

    energy = state.getPotentialEnergy()

    # get this as a float in kJ mol-1
    energy = energy.value_in_unit(energy.unit)

    # -74975.9 +/- 0.5
    assert mols.energy(map=map).to(sr.units.kJ_per_mol) == pytest.approx(
        energy, abs=0.5
    )


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
@pytest.mark.slow
def test_openmm_dynamics(ala_mols, openmm_platform):
    mols = ala_mols

    map = {
        "cutoff": 10 * sr.units.angstrom,
        "cutoff_type": "REACTION_FIELD",
        "dielectric": 78.0,
        "temperature": 25 * sr.units.celsius,
        "platform": openmm_platform,
        "constraint": "bonds-h-angles",
        # "pressure": 1 * sr.units.atm,   # currently disagree with energies for NPT...
    }

    sire_nrg = mols.energy(map=map)

    # Need to set constraints to `none` so that we
    # get energy agreement with sire - without this
    # we will be missing the bond energies
    # (and some angle energies)

    d = mols.dynamics(
        timestep=1 * sr.units.femtosecond,
        save_frequency=1 * sr.units.picosecond,
        map=map,
        constraint="none",
    )

    omm_nrg = d.current_potential_energy()

    assert sire_nrg.value() == pytest.approx(omm_nrg.value(), abs=0.5)

    assert d.ensemble().is_canonical()
    assert d.ensemble().temperature() == 25 * sr.units.celsius

    assert d.timestep() == 1 * sr.units.femtosecond

    d.run(0.1 * sr.units.picosecond, 0.01 * sr.units.picosecond)

    assert d.current_step() == 100
    # the molecules started from 6 ns
    assert d.current_time().to(sr.units.picosecond) == pytest.approx(6000.1)

    mols = d.commit()

    assert mols.num_frames() == 10

    sire_nrg = mols.energy(map=map)

    omm_nrg = d.current_potential_energy()

    assert sire_nrg.value() == pytest.approx(omm_nrg.value(), abs=0.5)


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_openmm_options(ala_mols, openmm_platform):
    mols = ala_mols

    mol = mols[0]

    m = {
        "integrator": "langevin_middle",
        "temperature": 25 * sr.units.celsius,
        "pressure": 1 * sr.units.atm,
        "friction": 5 / sr.units.picosecond,
        "platform": openmm_platform,
        "constraint": "bonds-h-angles",
    }

    omm = sr.convert.to(mol, "openmm", map=m)

    for platform in ["CPU", "OpenCL", "CUDA"]:
        m["platform"] = platform

        try:
            omm = sr.convert.to(mol, "openmm", map=m)
        except ValueError:
            # maybe OpenCL or CUDA are not supported
            pass


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_openmm_ignore_constrained(ala_mols, openmm_platform):
    mols = ala_mols

    mol = mols[0]

    d = mol.dynamics(
        constraint="bonds-h-angles",
        include_constrained_energies=True,
        platform=openmm_platform,
    )

    nrg1 = d.current_potential_energy()

    d = mol.dynamics(
        constraint="bonds-h-angles",
        include_constrained_energies=False,
        platform=openmm_platform,
    )

    nrg2 = d.current_potential_energy()

    # these two energies should be different, because
    # we should be ignoring the constrained bonds and angles
    assert abs(nrg2.value() - nrg1.value()) > 1.0


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_openmm_no_zero_sigmas(zero_lj_mols, openmm_platform):
    mols = zero_lj_mols

    omm = sr.convert.to(
        mols, "openmm", map={"constraint": "h-bonds", "platform": openmm_platform}
    )

    from openmm import XmlSerializer

    xml = XmlSerializer.serialize(omm.getSystem())

    for line in xml.split("\n"):
        assert 'sig="0"' not in line


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_openmm_skipped_constrained_bonds(zero_lj_mols, openmm_platform):
    mols = zero_lj_mols

    omm1 = sr.convert.to(
        mols,
        "openmm",
        map={
            "constraint": "h-bonds",
            "include_constrained_energies": True,
            "platform": openmm_platform,
        },
    )

    omm2 = sr.convert.to(
        mols,
        "openmm",
        map={
            "constraint": "h-bonds",
            "include_constrained_energies": False,
            "platform": openmm_platform,
        },
    )

    nrg1 = omm1.get_potential_energy().to(sr.units.kcal_per_mol)
    nrg2 = omm2.get_potential_energy().to(sr.units.kcal_per_mol)

    # Check the energies haven't changed
    # (regression check - here are the current values)
    assert nrg1 == pytest.approx(-447.44, 1e-3)
    assert nrg2 == pytest.approx(-3279.87, 1e-3)

    from openmm import XmlSerializer

    xml1 = XmlSerializer.serialize(omm1.getSystem())
    xml2 = XmlSerializer.serialize(omm2.getSystem())

    assert xml1 != xml2

    lines1 = xml1.split("\n")
    lines2 = xml2.split("\n")

    i = 0

    for j in range(0, len(lines2)):
        line1 = lines1[i]
        line2 = lines2[j]
        i += 1

        while line1 != line2:
            line1 = lines1[i]
            assert "Bond" in line1
            i += 1


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_openmm_default_box_vectors(ala_mols, openmm_platform):
    mols = ala_mols.clone()

    # Remove the shared space property.
    mols.remove_shared_property("space")

    # Create a "vacuum" system.
    omm = sr.convert.to(
        mols[0],
        "openmm",
        map={
            "cutoff": "12 A",
            "platform": openmm_platform,
        },
    )

    # Get the box vectors.
    box_vectors = omm.getState().getPeriodicBoxVectors()

    # Get the AABox of the first molecule.
    aabox = sr.legacy.Vol.AABox(mols[0].property("coordinates").to_vector())

    # Work out the box vectors from the AABox, i.e. the box size plus twice the cutoff.
    box = 2.0 * aabox.half_extents() + sr.maths.Vector(24, 24, 24)

    from openmm.unit import angstrom

    # Check that the box vectors match the expected values.
    for i, vec in enumerate(box_vectors):
        assert vec[i].value_in_unit(angstrom) == pytest.approx(box[i].value(), abs=1e-3)


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_openmm_membrane_barostat(ala_mols, openmm_platform):
    from openmm import MonteCarloMembraneBarostat
    from openmm import unit

    mols = ala_mols.clone()

    # Create a dynamics object with a membrane barostat.
    d = mols.dynamics(
        pressure="1atm",
        temperature="298K",
        surface_tension="1 angstrom*atm",
        platform=openmm_platform,
        barostat_frequency=50,
    )

    # Find the barostat.
    barostat = None
    for force in d.context().getSystem().getForces():
        if force.getName() == "MonteCarloMembraneBarostat":
            barostat = force
            break
    assert barostat is not None

    # Check the barostat parameters.
    assert barostat.getDefaultPressure().value_in_unit(
        unit.atmosphere
    ) == pytest.approx(1.0, abs=1e-3)
    assert barostat.getDefaultSurfaceTension().value_in_unit(
        unit.angstrom * unit.atmosphere
    ) == pytest.approx(1.0, abs=1e-3)
    assert barostat.getDefaultTemperature().value_in_unit(unit.kelvin) == pytest.approx(
        298.0, abs=1e-3
    )
    assert barostat.getXYMode() == MonteCarloMembraneBarostat.XYIsotropic
    assert barostat.getZMode() == MonteCarloMembraneBarostat.ZFree
    assert barostat.getFrequency() == 50


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_openmm_ghost_14_bug(ghost_14_bug, openmm_platform):
    """
    Regression test for a bug where a pair involving a ghost atom whose
    intrascale exclusion differs between end states (excluded at lambda=0,
    a real 1-4 scale at lambda=1, because the connectivity path between the
    two atoms only exists at lambda=1) never had a slot allocated in the
    Ghost14BondForce. The construction-time check that decides whether a
    pair needs a ghost-14 slot only looked at the lambda=0 exception scale,
    so a pair that's excluded at lambda=0 but not at lambda=1 was silently
    skipped, then permanently excluded from every nonbonded force for the
    lifetime of the simulation - even at lambda=1, where it should have a
    real interaction.

    The test molecule (cyclopropane -> propane, a minimal "ring-breaking"
    perturbation with no protein or water) reproduces this exactly: atom 2
    (a ring/chain carbon) and atoms 9, 10, 11 (new hydrogens that only
    exist on atom 0 at lambda=1, once the propane chain needs an extra H
    that the cyclopropane ring didn't) are 1-3 excluded at lambda=0 and a
    real 1-4 pair, (0.8333, 0.5), at lambda=1.
    """
    mols = sr.morph.link_to_reference(ghost_14_bug)
    mol = mols[0]

    d = mol.dynamics(
        schedule=sr.cas.LambdaSchedule.standard_morph(),
        platform=openmm_platform,
        cutoff="none",
    )

    def get_force(d, name):
        for force in d.context().getSystem().getForces():
            if force.getName() == name:
                return force
        return None

    # the bug-triggering pairs - atom 2 (C3) with the three new hydrogens
    # on atom 0 (C1) that only exist at lambda=1
    pairs = {(2, 9), (2, 10), (2, 11)}

    def get_ghost14_params(d):
        ghost14 = get_force(d, "Ghost14BondForce")
        assert ghost14 is not None

        found = {}
        for i in range(ghost14.getNumBonds()):
            p1, p2, params = ghost14.getBondParameters(i)
            key = (min(p1, p2), max(p1, p2))
            if key in pairs:
                found[key] = list(params)
        return found

    # at every lambda, all three pairs must have a slot in Ghost14BondForce
    # at all (this is what the bug broke - they were silently missing)
    d.set_lambda(0.0)
    params0 = get_ghost14_params(d)
    assert set(params0.keys()) == pairs

    d.set_lambda(0.5)
    params_mid = get_ghost14_params(d)
    assert set(params_mid.keys()) == pairs

    d.set_lambda(1.0)
    params1 = get_ghost14_params(d)
    assert set(params1.keys()) == pairs

    # at lambda=0 the pair is fully excluded, matching intrascale0 = (0,0)
    # (four_epsilon is parameter index 2)
    for params in params0.values():
        assert params[2] == pytest.approx(0.0, abs=1e-6)

    # the interaction should grow smoothly and monotonically as lambda
    # goes from 0 to 1, ending up clearly nonzero (matching the real
    # intrascale1 = (0.8333, 0.5) 1-4 scale)
    for key in pairs:
        assert params0[key][2] < params_mid[key][2] < params1[key][2]
        assert params1[key][2] > 0.1
