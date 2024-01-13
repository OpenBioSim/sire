import sire as sr
import pytest


def _run_test(mols, is_slow=False, use_taylor=False, precision=1e-3, platform="CPU"):
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
        "constraint": "bonds-h-angles",
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

    omm.set_lambda(0.0)
    assert omm.get_energy().value() == pytest.approx(nrg0, precision)

    omm.set_lambda(0.5)
    nrg0_5 = omm.get_energy().value()

    omm.set_lambda(1.0)
    assert omm.get_energy().value() == pytest.approx(nrg1, precision)

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


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_neopentane_methane_scan(neopentane_methane, openmm_platform):
    mols = neopentane_methane.clone()

    for mol in mols.molecules("property is_perturbable"):
        mols.update(mol.perturbation().link_to_reference().commit())

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
        constraint="h_bonds_not_perturbed",
        cutoff="10 A",
        include_constrained_energies=True,
        platform=openmm_platform,
    )

    calc_hbonds_not_perturbed = {}

    for lam_val, nrg in expected_hbonds_not_perturbed.items():
        d.set_lambda(lam_val)
        calc_hbonds_not_perturbed[lam_val] = d.current_potential_energy().value()

    d = mols.dynamics(
        constraint="h_bonds_not_perturbed",
        cutoff="10 A",
        include_constrained_energies=False,
        platform=openmm_platform,
    )

    calc_hbonds_not_perturbed_no_energy = {}

    for lam_val, nrg in expected_hbonds_not_perturbed.items():
        d.set_lambda(lam_val)
        calc_hbonds_not_perturbed_no_energy[
            lam_val
        ] = d.current_potential_energy().value()

    d = mols.dynamics(
        constraint="h_bonds_not_perturbed",
        cutoff="10 A",
        include_constrained_energies=False,
        swap_end_states=True,
        platform=openmm_platform,
    )

    calc_hbonds_not_perturbed_no_energy_swap = {}

    for lam_val, nrg in expected_hbonds_not_perturbed.items():
        d.set_lambda(lam_val)
        calc_hbonds_not_perturbed_no_energy_swap[
            lam_val
        ] = d.current_potential_energy().value()

    # should match the no_constraints somd energy at the end points
    assert calc_none[0.0] == pytest.approx(expected_none[0.0], 1e-3)
    assert calc_none[1.0] == pytest.approx(expected_none[1.0], 1e-3)

    assert calc_hbonds_not_perturbed[0.0] == pytest.approx(expected_none[0.0], 1e-3)
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
