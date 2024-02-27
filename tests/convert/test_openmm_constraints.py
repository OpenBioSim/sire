import sire as sr
import pytest


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_h_bond_constraints(merged_ethane_methanol, openmm_platform):
    mols = sr.morph.link_to_reference(merged_ethane_methanol)

    d = mols[0].dynamics(constraint="none", platform=openmm_platform)

    constraints = d.get_constraints()

    # no constraints
    assert len(constraints) == 0

    d = mols[0].dynamics(
        constraint="h-bonds", dynamic_constraints=False, platform=openmm_platform
    )

    constraints = d.get_constraints()

    # there are 6 bonds involving hydrogen
    assert len(constraints) == 6

    for constraint in constraints:
        # check for a single H atom
        assert constraint[0].atom("element H").element().symbol() == "H"

    d = mols[0].dynamics(
        constraint="bonds", dynamic_constraints=False, platform=openmm_platform
    )

    constraints = d.get_constraints()

    assert len(constraints) == len(mols[0].bonds())

    d = mols[0].dynamics(
        constraint="bonds", perturbable_constraint="none", platform=openmm_platform
    )

    # should be no constraints
    assert len(d.get_constraints()) == 0

    d = mols[0].dynamics(
        constraint="bonds",
        perturbable_constraint="h-bonds",
        ignore_perturbations=True,
        platform=openmm_platform,
    )

    # should be full bonds constraints
    assert len(d.get_constraints()) == len(mols[0].bonds())

    d = mols[0].dynamics(constraint="h-bonds-not-perturbed", platform=openmm_platform)

    constraints = d.get_constraints()

    # there should only be 2 bonds that are not constrained (C-C/O and C/O-H/G)
    assert len(constraints) == len(mols[0].bonds()) - 2

    for constraint in constraints:
        # check for a single H atom
        assert constraint[0].atom("element H").element().symbol() == "H"

    d = mols[0].dynamics(constraint="bonds-not-perturbed", platform=openmm_platform)

    # there should only be 2 bonds that are constrained (C-C/O and C/O-H/G)
    constraints = d.get_constraints()

    assert len(d.get_constraints()) == len(mols[0].bonds()) - 2

    d = mols[0].dynamics(
        constraint="none",
        perturbable_constraint="h-bonds-not-perturbed",
        platform=openmm_platform,
    )

    assert len(d.get_constraints()) == len(mols[0].bonds()) - 2

    d = mols[0].dynamics(
        constraint="none",
        perturbable_constraint="bonds-not-perturbed",
        platform=openmm_platform,
    )

    assert len(d.get_constraints()) == len(mols[0].bonds()) - 2


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_neo_constraints(neopentane_methane, openmm_platform):
    mols_fwds = sr.morph.link_to_reference(neopentane_methane)
    mols_bwds = sr.morph.link_to_perturbed(neopentane_methane)

    mols_fwds = sr.morph.repartition_hydrogen_masses(mols_fwds)
    mols_bwds = sr.morph.repartition_hydrogen_masses(mols_bwds)

    d_fwds = mols_fwds.dynamics(constraint="none", platform=openmm_platform)
    d_bwds = mols_bwds.dynamics(
        constraint="none", swap_end_states=True, platform=openmm_platform
    )

    c_fwds = d_fwds.get_constraints()
    c_bwds = d_bwds.get_constraints()

    assert len(c_fwds) == len(c_bwds) == 0

    d_fwds = mols_fwds.dynamics(constraint="bonds", platform=openmm_platform)
    d_bwds = mols_bwds.dynamics(
        constraint="bonds", swap_end_states=True, platform=openmm_platform
    )

    c_fwds = d_fwds.get_constraints()
    c_bwds = d_bwds.get_constraints()

    assert len(c_fwds) == len(c_bwds) == len(mols_fwds[0].bonds())

    for f, b in zip(c_fwds, c_bwds):
        assert f[0][0].name() == b[0][0].name()
        assert f[0][1].name() == b[0][1].name()

    d_fwds = mols_fwds.dynamics(constraint="h-bonds", platform=openmm_platform)
    d_bwds = mols_bwds.dynamics(
        constraint="h-bonds", swap_end_states=True, platform=openmm_platform
    )

    c_fwds = d_fwds.get_constraints()
    c_bwds = d_bwds.get_constraints()

    # this ends up constraining all bonds, as in either end state
    # they will always have a hydrogen
    assert len(c_fwds) == len(c_bwds) == len(mols_fwds[0].bonds())

    for f, b in zip(c_fwds, c_bwds):
        assert f[0][0].name() == b[0][0].name()
        assert f[0][1].name() == b[0][1].name()

    # there should be one dynamic constraint, of the C-C/H bond
    p = mols_fwds[0].perturbation().to_openmm(constraint="h-bonds")
    assert len(p.changed_constraints(to_pandas=False)) == 1

    d_fwds = mols_fwds.dynamics(
        constraint="bonds-not-perturbed", platform=openmm_platform
    )

    d_bwds = mols_bwds.dynamics(
        constraint="bonds-not-perturbed", swap_end_states=True, platform=openmm_platform
    )

    c_fwds = d_fwds.get_constraints()
    c_bwds = d_bwds.get_constraints()

    # there should be no dynamic constraints
    p = mols_fwds[0].perturbation().to_openmm(constraint="bonds-not-perturbed")
    assert len(p.changed_constraints(to_pandas=False)) == 0

    # the perturbing C-C/H bond should not be constrained
    assert len(c_fwds) == len(c_bwds) == len(mols_fwds[0].bonds()) - 1

    for f, b in zip(c_fwds, c_bwds):
        assert f[0][0].name() == b[0][0].name()
        assert f[0][1].name() == b[0][1].name()

    d_fwds = mols_fwds.dynamics(
        constraint="h-bonds-not-perturbed",
        platform=openmm_platform,
    )

    d_bwds = mols_bwds.dynamics(
        constraint="h-bonds-not-perturbed",
        platform=openmm_platform,
        swap_end_states=True,
    )

    c_fwds = d_fwds.get_constraints()
    c_bwds = d_bwds.get_constraints()

    assert len(c_fwds) == len(c_bwds)

    for f, b in zip(c_fwds, c_bwds):
        assert f[0][0].name() == b[0][0].name()
        assert f[0][1].name() == b[0][1].name()

    d_fwds = mols_fwds.dynamics(
        constraint="h-bonds-not-heavy-perturbed",
        platform=openmm_platform,
    )

    d_bwds = mols_bwds.dynamics(
        constraint="h-bonds-not-heavy-perturbed",
        platform=openmm_platform,
        swap_end_states=True,
    )

    c_fwds = d_fwds.get_constraints()
    c_bwds = d_bwds.get_constraints()

    # this should have constrained everything, as everything
    # is either a H bond or a perturbable bond containing H
    # at either end state
    assert len(c_fwds) == len(c_bwds) == len(mols_fwds[0].bonds())

    for f, b in zip(c_fwds, c_bwds):
        assert f[0][0].name() == b[0][0].name()
        assert f[0][1].name() == b[0][1].name()

    # check that there is only one perturbable constraint
    p = mols_fwds[0].perturbation().to_openmm(constraint="h-bonds-not-heavy-perturbed")
    assert (len(p.changed_constraints(to_pandas=False))) == 1


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_dynamic_constraints(merged_ethane_methanol, openmm_platform):
    mols = sr.morph.link_to_reference(merged_ethane_methanol)

    d = mols[0].dynamics(
        constraint="h-bonds",
        dynamic_constraints=True,
        platform=openmm_platform,
        lambda_value=0.0,
    )

    constraints = d.get_constraints()

    nrg = d.current_potential_energy().value()

    assert nrg == pytest.approx(1.65783, abs=0.001)

    # there are 6 bonds involving hydrogen - they should all be the same
    assert len(constraints) == 6

    for constraint in constraints:
        assert constraint[1].value() == pytest.approx(1.0969)

    d = mols[0].dynamics(
        constraint="h-bonds",
        dynamic_constraints=True,
        platform=openmm_platform,
        lambda_value=1.0,
    )

    constraints = d.get_constraints()

    # there are 6 bonds involving hydrogen - they should all be the same
    assert len(constraints) == 6

    # one of the constraints should be 0.973 A
    cons = {}

    for constraint in constraints:
        dist = constraint[1].value()

        if dist not in cons:
            cons[dist] = 0

        cons[dist] += 1

    assert len(cons) == 2

    for key, value in cons.items():
        if value == 1:
            assert key == pytest.approx(0.973)
        elif value == 5:
            assert key == pytest.approx(1.0969)
        else:
            assert False

    # make sure that the bond parameters are correct - can only do this
    # via potential energy
    nrg = d.current_potential_energy()

    assert nrg.value() == pytest.approx(13.8969, abs=0.001)
