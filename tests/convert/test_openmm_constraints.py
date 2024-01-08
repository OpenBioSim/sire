import sire as sr
import pytest


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_h_bond_constraints(merged_ethane_methanol, openmm_platform):
    mols = merged_ethane_methanol.clone()

    for mol in mols.molecules("property is_perturbable"):
        mols.update(mol.perturbation().link_to_reference().commit())

    d = mols[0].dynamics(constraint="none", platform=openmm_platform)

    constraints = d.get_constraints()

    # no constraints
    assert len(constraints) == 0

    d = mols[0].dynamics(constraint="h-bonds", platform=openmm_platform)

    constraints = d.get_constraints()

    # there are 6 bonds involving hydrogen
    assert len(constraints) == 6

    for constraint in constraints:
        # check for a single H atom
        assert constraint[0].atom("element H").element().symbol() == "H"

        # check that the constraint distance is close to the
        # current distance
        assert sr.measure(constraint[0][0], constraint[0][1]).value() == pytest.approx(
            constraint[1].value()
        )

    d = mols[0].dynamics(constraint="bonds", platform=openmm_platform)

    constraints = d.get_constraints()

    assert len(constraints) == len(mols[0].bonds())

    for constraint in constraints:
        # check that the constraint distance is close to the
        # current distance
        assert sr.measure(constraint[0][0], constraint[0][1]).value() == pytest.approx(
            constraint[1].value()
        )

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

    # there should only be 3 bonds that are constrained
    assert len(constraints) == 3

    for constraint in constraints:
        # check for a single H atom
        assert constraint[0].atom("element H").element().symbol() == "H"

        # check that the constraint distance is close to the
        # current distance
        assert sr.measure(constraint[0][0], constraint[0][1]).value() == pytest.approx(
            constraint[1].value()
        )

    d = mols[0].dynamics(constraint="bonds-not-perturbed", platform=openmm_platform)

    # there should only be 3 bonds that are constrained
    constraints = d.get_constraints()

    assert len(d.get_constraints()) == 3
