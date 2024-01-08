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
        print(constraint[0])
        h = constraint[0].atom("element H")
        print(h)

        assert len(h) == 1
        assert h[0].element().symbol() == "H"
