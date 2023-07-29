import sire as sr
import pytest


def _run_test(mols, is_slow=False):
    c = mols.cursor()

    # can only get the same energies if they have the same coordinates
    c["molidx 0"]["coordinates"] = c["molidx 0"]["coordinates0"]
    c["molidx 0"]["coordinates1"] = c["molidx 0"]["coordinates0"]

    mols = c.commit()

    # in this case the perturbable molecule is the first one
    merge = mols[0]

    # and everything else is water
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
    mols = merge + water

    # these are the non-perturbable systems at lambda 0 and lambda 1
    mols0 = get_end_state(merge, "0", "1") + water
    mols1 = get_end_state(merge, "1", "0") + water

    # create the perturbable OpenMM system
    omm = sr.convert.to(mols, "openmm")

    # a very basic lambda schedule
    l = omm.get_lambda_schedule()
    l.add_stage("morph", (1 - l.lam()) * l.initial() + l.lam() * l.final())
    omm.set_lambda_schedule(l)

    def e(nrg):
        return nrg.value_in_unit(nrg.unit)

    # now the lambda0 and lambda1 non-perturbable end states
    omm0 = sr.convert.to(mols0, "openmm")
    nrg0 = e(omm0.get_energy())

    omm1 = sr.convert.to(mols1, "openmm")
    nrg1 = e(omm1.get_energy())

    omm.set_lambda(0.0)
    assert e(omm.get_energy()) == pytest.approx(nrg0)

    omm.set_lambda(0.5)
    nrg0_5 = e(omm.get_energy())

    omm.set_lambda(1.0)
    assert e(omm.get_energy()) == pytest.approx(nrg1)

    omm.set_lambda(0.5)
    assert e(omm.get_energy()) == pytest.approx(nrg0_5)

    omm.set_lambda(0.0)
    assert e(omm.get_energy()) == pytest.approx(nrg0)

    # now swap the end states - lambda 0 == 1 and lambda 1 == 0
    omm = sr.convert.to(mols, "openmm", map={"swap_end_states": True})
    omm.set_lambda_schedule(l)

    omm.set_lambda(0.0)
    assert e(omm.get_energy()) == pytest.approx(nrg1)

    omm.set_lambda(0.5)
    assert e(omm.get_energy()) == pytest.approx(nrg0_5)

    omm.set_lambda(1.0)
    assert e(omm.get_energy()) == pytest.approx(nrg0)

    omm.set_lambda(0.5)
    assert e(omm.get_energy()) == pytest.approx(nrg0_5)

    omm.set_lambda(0.0)
    assert e(omm.get_energy()) == pytest.approx(nrg1)


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_openmm_scale_lambda(merged_molecule):
    _run_test(merged_molecule.clone(), False)


@pytest.mark.veryslow
@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_big_openmm_scale_lambda(merged_molecule):
    _run_test(merged_molecule.clone(), True)
