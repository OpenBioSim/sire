import sire as sr
import pytest


def _run_test(mols, is_slow=False, use_taylor=False):
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
        "platform": "cpu",
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
    assert omm.get_energy().value() == pytest.approx(nrg0)

    omm.set_lambda(0.5)
    nrg0_5 = omm.get_energy().value()

    omm.set_lambda(1.0)
    assert omm.get_energy().value() == pytest.approx(nrg1)

    omm.set_lambda(0.5)
    assert omm.get_energy().value() == pytest.approx(nrg0_5)

    omm.set_lambda(0.0)
    assert omm.get_energy().value() == pytest.approx(nrg0)

    # now swap the end states - lambda 0 == 1 and lambda 1 == 0
    map["swap_end_states"] = True

    omm = sr.convert.to(mols, "openmm", map=map)
    omm.set_lambda_schedule(l)

    omm.set_lambda(0.0)
    assert omm.get_energy().value() == pytest.approx(nrg1)

    omm.set_lambda(0.5)
    assert omm.get_energy().value() == pytest.approx(nrg0_5)

    omm.set_lambda(1.0)
    assert omm.get_energy().value() == pytest.approx(nrg0)

    omm.set_lambda(0.5)
    assert omm.get_energy().value() == pytest.approx(nrg0_5)

    omm.set_lambda(0.0)
    assert omm.get_energy().value() == pytest.approx(nrg1)


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_openmm_scale_lambda_simple(merged_ethane_methanol):
    _run_test(merged_ethane_methanol.clone(), False)


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_openmm_scale_lambda_taylor_simple(merged_ethane_methanol):
    _run_test(merged_ethane_methanol.clone(), False, True)


@pytest.mark.veryslow
@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_big_openmm_scale_lambda_simple(merged_ethane_methanol):
    _run_test(merged_ethane_methanol.clone(), True)


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_openmm_scale_lambda_ligand(merged_zan_ose):
    _run_test(merged_zan_ose.clone(), False)


@pytest.mark.veryslow
@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_big_openmm_scale_lambda_ligand(merged_zan_ose):
    _run_test(merged_zan_ose.clone(), True)


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_openmm_scale_lambda_dichloroethane(ethane_12dichloroethane):
    _run_test(ethane_12dichloroethane.clone(), False)


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_openmm_scale_lambda_cyclopentane(pentane_cyclopentane):
    mols = pentane_cyclopentane.clone()

    for mol in mols.molecules("property is_perturbable"):
        mol = mol.edit().add_link("connectivity", "connectivity0").commit()
        mols.update(mol)

    _run_test(mols, False)
