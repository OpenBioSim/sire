import sire as sr
import pytest


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_dynamic_constraints(merged_ethane_methanol, openmm_platform):
    mols = merged_ethane_methanol.clone()

    for mol in mols.molecules("property is_perturbable"):
        mols.update(mol.perturbation().link_to_reference().commit())

    d = mols[0].dynamics(constraint="bonds", platform=openmm_platform)

    d.set_lambda(0.0)

    constraints_0 = d.get_constraints()

    def get_r0(bond, state="0"):
        return sr.legacy.MM.AmberBond(
            bond.potential(map={"bond": f"bond{state}"}), sr.legacy.CAS.r
        ).r0()

    for constraint in constraints_0:
        bond = sr.mm.Bond(constraint[0][0], constraint[0][1])
        r0_0 = get_r0(bond, "0")
        assert constraint[1].value() == pytest.approx(r0_0)

    d.set_lambda(1.0)

    constraints_1 = d.get_constraints()

    assert len(constraints_0) == len(constraints_1)

    for constraint in constraints_1:
        bond = sr.mm.Bond(constraint[0][0], constraint[0][1])
        r0_1 = get_r0(bond, "1")
        assert constraint[1].value() == pytest.approx(r0_1)

    d.set_lambda(0.5)

    constraints_05 = d.get_constraints()

    assert len(constraints_0) == len(constraints_05)

    for constraint in constraints_05:
        bond = sr.mm.Bond(constraint[0][0], constraint[0][1])
        r0_0 = get_r0(bond, "0")
        r0_1 = get_r0(bond, "1")
        assert constraint[1].value() == pytest.approx(0.5 * (r0_0 + r0_1))

    d.set_lambda(0.0, update_constraints=False)

    assert constraints_05 == d.get_constraints()

    d = mols[0].dynamics(
        constraint="bonds", dynamic_constraints=False, platform=openmm_platform
    )

    d.set_lambda(0.0)

    constraints_0 = d.get_constraints()

    for constraint in constraints_0:
        bond = sr.mm.Bond(constraint[0][0], constraint[0][1])
        r0_0 = get_r0(bond, "0")
        assert constraint[1].value() == pytest.approx(r0_0)

    d.set_lambda(1.0)

    constraints_1 = d.get_constraints()

    for constraint in constraints_1:
        bond = sr.mm.Bond(constraint[0][0], constraint[0][1])
        r0_0 = get_r0(bond, "0")
        assert constraint[1].value() == pytest.approx(r0_0)
        assert constraint in constraints_0

    d.set_lambda(0.5)

    constraints_05 = d.get_constraints()

    for constraint in constraints_05:
        bond = sr.mm.Bond(constraint[0][0], constraint[0][1])
        r0_0 = get_r0(bond, "0")
        assert constraint[1].value() == pytest.approx(r0_0)
        assert constraint in constraints_0
