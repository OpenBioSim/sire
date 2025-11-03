import sire as sr
import pytest


def _run_test(mols, is_slow=False, use_taylor=False, precision=1e-3, platform="CPU"):
    try:
        space = mols.space()
    except Exception:
        space = sr.vol.Cartesian()

    c = mols.cursor()

    # can only get the same energies if they have the same coordinates
    c["molidx 0"]["coordinates"] = c["molidx 0"]["coordinates0"]
    c["molidx 0"]["coordinates1"] = c["molidx 0"]["coordinates0"]

    mols = c.commit()

    has_water = len(mols) > 1

    # in this case the perturbable molecule is the first one
    merge = mols[0]

    # and everything else is water
    if has_water:
        if is_slow:
            water = mols["water"]
        else:
            water = mols["closest 5 waters to molidx 0"]

    # this function will extract the lambda0 or lambda1 end state
    def get_end_state(mol, state, remove_state):
        c = mol.cursor()
        for key in c.keys():
            if key.endswith(state):
                c[key.removesuffix(state)] = c[key]
                del c[key]
            elif key.endswith(remove_state):
                del c[key]

        c["is_perturbable"] = False
        return c.commit()

    # this is the perturbable system
    if has_water:
        mols = merge + water

        # these are the non-perturbable systems at lambda 0 and lambda 1
        mols0 = get_end_state(merge, "0", "1") + water
        mols1 = get_end_state(merge, "1", "0") + water
    else:
        mols = merge

        # these are the non-perturbable systems at lambda 0 and lambda 1
        mols0 = get_end_state(merge, "0", "1")
        mols1 = get_end_state(merge, "1", "0")

    # a very basic lambda schedule
    l = sr.cas.LambdaSchedule()
    l.add_stage("morph", (1 - l.lam()) * l.initial() + l.lam() * l.final())

    # need to use the reference platform on GH Actions
    map = {
        "platform": platform,
        "schedule": l,
        "constraint": "h-bonds-not-perturbed",
        "include_constrained_energies": True,
        "dynamic_constraints": False,
        "space": space,
    }

    if use_taylor:
        map["use_taylor_softening"] = True

    # create the perturbable OpenMM system
    omm = sr.convert.to(mols, "openmm", map=map)

    # now the lambda0 and lambda1 non-perturbable end states
    omm0 = sr.convert.to(mols0, "openmm", map=map)
    nrg0 = omm0.get_energy().value()

    omm1 = sr.convert.to(mols1, "openmm", map=map)
    nrg1 = omm1.get_energy().value()

    omm.set_lambda(1.0)
    assert omm.get_energy().value() == pytest.approx(nrg1, precision)

    omm.set_lambda(0.0)
    assert omm.get_energy().value() == pytest.approx(nrg0, precision)

    omm.set_lambda(0.5)
    nrg0_5 = omm.get_energy().value()

    omm.set_lambda(0.5)
    assert omm.get_energy().value() == pytest.approx(nrg0_5, precision)

    omm.set_lambda(0.0)
    assert omm.get_energy().value() == pytest.approx(nrg0, precision)

    omm.set_lambda(0.5)
    assert omm.get_energy().value() == pytest.approx(nrg0_5, precision)

    omm.set_lambda(1.0)
    assert omm.get_energy().value() == pytest.approx(nrg1, precision)

    # now swap the end states - lambda 0 == 1 and lambda 1 == 0
    map["swap_end_states"] = True

    omm = sr.convert.to(mols, "openmm", map=map)
    omm.set_lambda_schedule(l)

    omm.set_lambda(0.0)
    assert omm.get_energy().value() == pytest.approx(nrg1, precision)

    omm.set_lambda(0.5)
    assert omm.get_energy().value() == pytest.approx(nrg0_5, precision)

    omm.set_lambda(1.0)
    assert omm.get_energy().value() == pytest.approx(nrg0, precision)

    omm.set_lambda(0.5)
    assert omm.get_energy().value() == pytest.approx(nrg0_5, precision)

    omm.set_lambda(0.0)
    assert omm.get_energy().value() == pytest.approx(nrg1, precision)

    omm.set_lambda(0.5)
    assert omm.get_energy().value() == pytest.approx(nrg0_5, precision)

    omm.set_lambda(1.0)
    assert omm.get_energy().value() == pytest.approx(nrg0, precision)


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_openmm_scale_lambda_simple(merged_ethane_methanol, openmm_platform):
    _run_test(merged_ethane_methanol.clone(), False, platform=openmm_platform)


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_openmm_scale_lambda_taylor_simple(merged_ethane_methanol, openmm_platform):
    _run_test(merged_ethane_methanol.clone(), False, True, platform=openmm_platform)


@pytest.mark.veryslow
@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_big_openmm_scale_lambda_simple(merged_ethane_methanol, openmm_platform):
    _run_test(merged_ethane_methanol.clone(), True, platform=openmm_platform)


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_openmm_scale_lambda_ligand(merged_zan_ose, openmm_platform):
    _run_test(merged_zan_ose.clone(), False, platform=openmm_platform)


@pytest.mark.veryslow
@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_big_openmm_scale_lambda_ligand(merged_zan_ose, openmm_platform):
    _run_test(merged_zan_ose.clone(), True, platform=openmm_platform)


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_openmm_scale_lambda_dichloroethane(ethane_12dichloroethane, openmm_platform):
    _run_test(ethane_12dichloroethane.clone(), False, platform=openmm_platform)


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_openmm_scale_lambda_cyclopentane(pentane_cyclopentane, openmm_platform):
    mols = pentane_cyclopentane.clone()

    for mol in mols.molecules("property is_perturbable"):
        mol = mol.edit().add_link("connectivity", "connectivity0").commit()
        mols.update(mol)

    _run_test(mols, False, platform=openmm_platform)


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_openmm_scale_lambda_neopentane_methane(neopentane_methane, openmm_platform):
    _run_test(neopentane_methane, False, platform=openmm_platform)


@pytest.mark.veryslow
@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_big_openmm_scale_lambda_neopentane_methane_solv(
    solvated_neopentane_methane, openmm_platform
):
    _run_test(solvated_neopentane_methane, True, platform=openmm_platform)


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_openmm_scale_lambda_neopentane_methane_solv(
    solvated_neopentane_methane, openmm_platform
):
    _run_test(solvated_neopentane_methane, False, platform=openmm_platform)


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_solvated_neopentane_methane_scan(solvated_neopentane_methane, openmm_platform):
    mols = sr.morph.link_to_reference(solvated_neopentane_methane)

    mols = sr.morph.repartition_hydrogen_masses(mols)

    mols = sr.morph.zero_ghost_torsions(mols)

    # these were calculated using somd, no constraints, no cutoff, no space
    expected = {
        0.0: -6845.34,
        1.0: -6759.75,
    }

    d = mols.dynamics(
        constraint="none",
        cutoff="none",
        cutoff_type="RF",
        platform=openmm_platform,
        map={"space": sr.vol.Cartesian()},
    )

    for lam_val, nrg in expected.items():
        d.set_lambda(lam_val)
        calc_nrg = d.current_potential_energy().value()
        assert calc_nrg == pytest.approx(nrg, abs=1e-2)

    # these were calculated using somd, no constraints, 10 A cutoff
    expected = {
        0.0: -8330.43,
        1.0: -8242.37,
    }

    d = mols.dynamics(
        constraint="none",
        cutoff="10 A",
        cutoff_type="RF",
        platform=openmm_platform,
    )

    for lam_val, nrg in expected.items():
        d.set_lambda(lam_val)
        calc_nrg = d.current_potential_energy().value()
        assert calc_nrg == pytest.approx(nrg, abs=1e-2)

    # these were calculated using somd, h-bonds-not-perturbed, 10 A cutoff
    expected = {
        0.0: -8331.14,
        1.0: -8243.08,
    }

    d = mols.dynamics(
        constraint="bonds_not_perturbed",
        cutoff="10 A",
        cutoff_type="RF",
        include_constrained_energies=False,
        platform=openmm_platform,
    )

    for lam_val, nrg in expected.items():
        d.set_lambda(lam_val)
        calc_nrg = d.current_potential_energy().value()
        assert calc_nrg == pytest.approx(nrg, abs=1e-2)


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_neopentane_methane_scan(neopentane_methane, openmm_platform):
    mols = sr.morph.link_to_reference(neopentane_methane)

    mols = sr.morph.repartition_hydrogen_masses(mols)

    # these were calculated using somd, no constraints
    expected_none = {
        0.0: -2.85704,
        0.1: -3.98964,
        0.2: -1.29653,
        0.3: 2.8749,
        0.4: 8.20206,
        0.5: 14.704,
        0.6: 22.4421,
        0.7: 31.4726,
        0.8: 41.8432,
        0.9: 53.5963,
        1.0: 66.773,
    }

    # these were calculated using somd, hbonds
    expected_hbonds = {
        0.0: -3.70711,
        0.1: -5.5007,
        0.2: -4.67283,
        0.3: -3.60548,
        0.4: -2.65583,
        0.5: -1.83948,
        0.6: -1.12953,
        0.7: -0.504399,
        0.8: 0.0490519,
        0.9: 0.538499,
        1.0: 0.970649,
    }

    # these were calculated using somd, hbonds_notperturbed
    expected_hbonds_not_perturbed = {
        0.0: -3.70499,
        0.1: -4.83759,
        0.2: -2.14448,
        0.3: 2.02695,
        0.4: 7.35411,
        0.5: 13.856,
        0.6: 21.5941,
        0.7: 30.6246,
        0.8: 40.9953,
        0.9: 52.7483,
        1.0: 65.9251,
    }

    d = mols.dynamics(constraint="none", cutoff="10 A", platform=openmm_platform)

    calc_none = {}

    for lam_val, nrg in expected_none.items():
        d.set_lambda(lam_val)
        calc_none[lam_val] = d.current_potential_energy().value()

    d = mols.dynamics(
        constraint="bonds_not_perturbed",
        cutoff="10 A",
        include_constrained_energies=True,
        platform=openmm_platform,
    )

    calc_hbonds_not_perturbed = {}

    for lam_val, nrg in expected_hbonds_not_perturbed.items():
        d.set_lambda(lam_val)
        calc_hbonds_not_perturbed[lam_val] = d.current_potential_energy().value()

    d = mols.dynamics(
        constraint="bonds_not_perturbed",
        cutoff="10 A",
        include_constrained_energies=False,
        platform=openmm_platform,
    )

    calc_hbonds_not_perturbed_no_energy = {}

    for lam_val, nrg in expected_hbonds_not_perturbed.items():
        d.set_lambda(lam_val)
        calc_hbonds_not_perturbed_no_energy[lam_val] = (
            d.current_potential_energy().value()
        )

    d = mols.dynamics(
        constraint="bonds-not-perturbed",
        cutoff="10 A",
        include_constrained_energies=False,
        swap_end_states=True,
        platform=openmm_platform,
    )

    calc_hbonds_not_perturbed_no_energy_swap = {}

    for lam_val, nrg in expected_hbonds_not_perturbed.items():
        d.set_lambda(lam_val)
        calc_hbonds_not_perturbed_no_energy_swap[lam_val] = (
            d.current_potential_energy().value()
        )

    d = mols.dynamics(
        constraint="h-bonds",
        cutoff="10 A",
        include_constrained_energies=False,
        platform=openmm_platform,
    )

    calc_hbonds_no_energy = {}

    for lam_val, nrg in expected_hbonds.items():
        d.set_lambda(lam_val)
        calc_hbonds_no_energy[lam_val] = d.current_potential_energy().value()

    # should match the no_constraints somd energy at the end points
    assert calc_none[0.0] == pytest.approx(expected_none[0.0], 1e-3)
    assert calc_none[1.0] == pytest.approx(expected_none[1.0], 1e-3)

    assert calc_hbonds_not_perturbed[0.0] == pytest.approx(expected_none[0.0], 1e-3)
    assert calc_hbonds_not_perturbed[1.0] == pytest.approx(expected_none[1.0], 1e-3)

    assert calc_hbonds_no_energy[0.0] == pytest.approx(expected_hbonds[0.0], 1e-2)
    assert calc_hbonds_no_energy[1.0] == pytest.approx(expected_hbonds[1.0], 1e-2)

    # but the hbonds_not_perturbed energy should be different if constraints
    # are not included - should equal to the somd constraints energy
    # (somd does not calculate energies of constrained bonds)
    assert calc_hbonds_not_perturbed_no_energy[0.0] == pytest.approx(
        expected_hbonds_not_perturbed[0.0], 1e-2
    )
    assert calc_hbonds_not_perturbed_no_energy[1.0] == pytest.approx(
        expected_hbonds_not_perturbed[1.0], 1e-2
    )

    assert calc_hbonds_not_perturbed_no_energy_swap[0.0] == pytest.approx(
        expected_hbonds_not_perturbed[1.0], 1e-2
    )

    assert calc_hbonds_not_perturbed_no_energy_swap[1.0] == pytest.approx(
        expected_hbonds_not_perturbed[0.0], 1e-2
    )

    for lam_val in expected_none.keys():
        # Including the energy should give the same energy regardless
        # of whether constraints are included or not
        assert calc_none[lam_val] == pytest.approx(
            calc_hbonds_not_perturbed[lam_val], 1e-3
        )

        # But not including the constraints should give a different energy
        assert calc_none[lam_val] != calc_hbonds_not_perturbed_no_energy[lam_val]

    # The paths through lambda space for somd and sire will be slightly
    # different - but should be consistently different comparing
    # including and not including constraints
    for lam_val in expected_none.keys():
        assert calc_none[lam_val] - calc_hbonds_not_perturbed_no_energy[
            lam_val
        ] == pytest.approx(
            expected_none[lam_val] - expected_hbonds_not_perturbed[lam_val], 1e-2
        )

    # check backwards is the reverse of forwards
    lamvals_f = list(calc_hbonds_not_perturbed_no_energy.keys())
    lamvals_b = list(calc_hbonds_not_perturbed_no_energy_swap.keys())

    lamvals_b.reverse()

    for lam_f, lam_b in zip(lamvals_f, lamvals_b):
        assert calc_hbonds_not_perturbed_no_energy[lam_f] == pytest.approx(
            calc_hbonds_not_perturbed_no_energy_swap[lam_b], 1e-3
        )


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_neopentane_methane_scan_no_cutoff(neopentane_methane, openmm_platform):
    mols = sr.morph.link_to_reference(neopentane_methane)

    mols = sr.morph.repartition_hydrogen_masses(mols)

    # these were calculated using somd, no constraints
    expected_none = {
        0.0: 0.0158906,
        0.1: -1.55981,
        0.2: 0.72506,
        0.3: 4.523,
        0.4: 9.51135,
        0.5: 15.709,
        0.6: 23.1774,
        0.7: 31.9727,
        0.8: 42.1424,
        0.9: 53.7287,
        1.0: 66.773,
    }

    # these were calculated using somd, hbonds_notperturbed
    expected_hbonds_not_perturbed = {
        0.0: -0.832059,
        0.1: -2.40776,
        0.2: -0.12289,
        0.3: 3.67505,
        0.4: 8.6634,
        0.5: 14.8611,
        0.6: 22.3295,
        0.7: 31.1247,
        0.8: 41.2944,
        0.9: 52.8808,
        1.0: 65.9251,
    }

    d = mols.dynamics(constraint="none", cutoff="infinite", platform=openmm_platform)

    calc_none = {}

    for lam_val, nrg in expected_none.items():
        d.set_lambda(lam_val)
        calc_none[lam_val] = d.current_potential_energy().value()

    d = mols.dynamics(
        constraint="bonds_not_perturbed",
        cutoff="infinite",
        include_constrained_energies=True,
        platform=openmm_platform,
    )

    calc_hbonds_not_perturbed = {}

    for lam_val, nrg in expected_hbonds_not_perturbed.items():
        d.set_lambda(lam_val)
        calc_hbonds_not_perturbed[lam_val] = d.current_potential_energy().value()

    d = mols.dynamics(
        constraint="bonds_not_perturbed",
        cutoff="infinite",
        include_constrained_energies=False,
        platform=openmm_platform,
    )

    calc_hbonds_not_perturbed_no_energy = {}

    for lam_val, nrg in expected_hbonds_not_perturbed.items():
        d.set_lambda(lam_val)
        calc_hbonds_not_perturbed_no_energy[lam_val] = (
            d.current_potential_energy().value()
        )

    d = mols.dynamics(
        constraint="bonds-not-perturbed",
        cutoff="infinite",
        include_constrained_energies=False,
        swap_end_states=True,
        platform=openmm_platform,
    )

    calc_hbonds_not_perturbed_no_energy_swap = {}

    for lam_val, nrg in expected_hbonds_not_perturbed.items():
        d.set_lambda(lam_val)
        calc_hbonds_not_perturbed_no_energy_swap[lam_val] = (
            d.current_potential_energy().value()
        )

    # should match the no_constraints somd energy at the end points
    assert calc_none[0.0] == pytest.approx(expected_none[0.0], abs=1e-2)
    assert calc_none[1.0] == pytest.approx(expected_none[1.0], 1e-2)

    assert calc_hbonds_not_perturbed[0.0] == pytest.approx(expected_none[0.0], abs=1e-3)
    assert calc_hbonds_not_perturbed[1.0] == pytest.approx(expected_none[1.0], 1e-3)

    # but the hbonds_not_perturbed energy should be different if constraints
    # are not included - should equal to the somd constraints energy
    # (somd does not calculate energies of constrained bonds)
    assert calc_hbonds_not_perturbed_no_energy[0.0] == pytest.approx(
        expected_hbonds_not_perturbed[0.0], 1e-2
    )
    assert calc_hbonds_not_perturbed_no_energy[1.0] == pytest.approx(
        expected_hbonds_not_perturbed[1.0], 1e-2
    )

    assert calc_hbonds_not_perturbed_no_energy_swap[0.0] == pytest.approx(
        expected_hbonds_not_perturbed[1.0], 1e-2
    )

    assert calc_hbonds_not_perturbed_no_energy_swap[1.0] == pytest.approx(
        expected_hbonds_not_perturbed[0.0], 1e-2
    )

    for lam_val in expected_none.keys():
        # Including the energy should give the same energy regardless
        # of whether constraints are included or not
        assert calc_none[lam_val] == pytest.approx(
            calc_hbonds_not_perturbed[lam_val], 1e-3
        )

        # But not including the constraints should give a different energy
        assert calc_none[lam_val] != calc_hbonds_not_perturbed_no_energy[lam_val]

    # The paths through lambda space for somd and sire will be slightly
    # different - but should be consistently different comparing
    # including and not including constraints
    for lam_val in expected_none.keys():
        assert calc_none[lam_val] - calc_hbonds_not_perturbed_no_energy[
            lam_val
        ] == pytest.approx(
            expected_none[lam_val] - expected_hbonds_not_perturbed[lam_val], 1e-2
        )

    # check backwards is the reverse of forwards
    lamvals_f = list(calc_hbonds_not_perturbed_no_energy.keys())
    lamvals_b = list(calc_hbonds_not_perturbed_no_energy_swap.keys())

    lamvals_b.reverse()

    for lam_f, lam_b in zip(lamvals_f, lamvals_b):
        assert calc_hbonds_not_perturbed_no_energy[lam_f] == pytest.approx(
            calc_hbonds_not_perturbed_no_energy_swap[lam_b], 1e-3
        )


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_neopentane_methane_no_charge(neopentane_methane, openmm_platform):
    mols = sr.morph.link_to_reference(neopentane_methane)

    mols = sr.morph.repartition_hydrogen_masses(mols)

    # these were calculated using somd, hbonds
    expected_hbonds = {
        0.0: 2.69409,
        0.1: -0.00735767,
        0.2: -0.0258584,
        0.3: 0.2558,
        0.4: 0.479681,
        0.5: 0.629462,
        0.6: 0.731321,
        0.7: 0.806166,
        0.8: 0.866469,
        0.9: 0.919267,
        1.0: 0.970649,
    }

    # these were calculated using somd, hbonds-notperturbed
    expected_hbonds_not_perturbed = {
        0.0: 2.69621,
        0.1: 0.655753,
        0.2: 2.50249,
        0.3: 5.88823,
        0.4: 10.4896,
        0.5: 16.3249,
        0.6: 23.455,
        0.7: 31.9352,
        0.8: 41.8127,
        0.9: 53.1291,
        1.0: 65.9251,
    }

    # these were calculated using somd, no constraints
    expected_none = {
        0.0: 3.54416,
        0.1: 1.5037,
        0.2: 3.35044,
        0.3: 6.73618,
        0.4: 11.3376,
        0.5: 17.1729,
        0.6: 24.3029,
        0.7: 32.7831,
        0.8: 42.6606,
        0.9: 53.977,
        1.0: 66.773,
    }

    d = mols.dynamics(
        constraint="none",
        cutoff="10 A",
        platform=openmm_platform,
    )

    # Use the schedule to set all charges to zero
    s = d.get_schedule()
    s.set_equation(stage="morph", lever="charge", equation=0.0)
    d.set_schedule(s)

    for lam_val, nrg in expected_none.items():
        d.set_lambda(lam_val)
        calc = d.current_potential_energy().value()
        assert calc == pytest.approx(nrg, abs=1e-4)

    d = mols.dynamics(
        constraint="h-bonds",
        cutoff="10 A",
        include_constrained_energies=False,
        platform=openmm_platform,
    )

    # Use the schedule to set all charges to zero
    s = d.get_schedule()
    s.set_equation(stage="morph", lever="charge", equation=0.0)
    d.set_schedule(s)

    for lam_val, nrg in expected_hbonds.items():
        d.set_lambda(lam_val)
        calc = d.current_potential_energy().value()
        assert calc == pytest.approx(nrg, abs=1e-4)

    d = mols.dynamics(
        constraint="h-bonds-not-perturbed",
        cutoff="10 A",
        include_constrained_energies=False,
        platform=openmm_platform,
    )

    # Use the schedule to set all charges to zero
    s = d.get_schedule()
    s.set_equation(stage="morph", lever="charge", equation=0.0)
    d.set_schedule(s)

    for lam_val, nrg in expected_hbonds_not_perturbed.items():
        d.set_lambda(lam_val)
        calc = d.current_potential_energy().value()
        assert calc == pytest.approx(nrg, abs=1e-4)
