import sire as sr
import pytest


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_perturbable_constraints(merged_ethane_methane, openmm_platform):
    mols = sr.morph.link_to_reference(merged_ethane_methane)

    d = mols[0].dynamics(constraint="none", platform=openmm_platform)

    constraints = d.get_constraints()

    # no constraints
    assert len(constraints) == 0

    d = mols[0].dynamics(constraint="bonds", platform=openmm_platform)

    constraints = d.get_constraints()

    # all bonds should be constrained
    assert len(constraints) == len(mols[0].bonds())

    d = mols[0].dynamics(
        constraint="bonds", perturbable_constraint="none", platform=openmm_platform
    )

    # should be no constraints
    assert len(d.get_constraints()) == 0

    d = mols[0].dynamics(
        constraint="bonds",
        perturbable_constraint="h-bonds",
        platform=openmm_platform,
    )

    constraints = d.get_constraints()

    # there are 7 bonds involving hydrogen (includes the C-C to C-H bond)
    assert len(constraints) == 7

    d = mols[0].dynamics(constraint="h-bonds", platform=openmm_platform)

    constraints = d.get_constraints()

    # there are 7 bonds involving hydrogen (includes the C-C to C-H bond)
    assert len(constraints) == 7

    d = mols[0].dynamics(constraint="bonds-not-perturbed", platform=openmm_platform)

    constraints = d.get_constraints()

    # this should be the 6 C-H bonds, as the central C-C to C-H bond is not constrained
    assert len(constraints) == 6

    d = mols[0].dynamics(constraint="h-bonds-not-perturbed", platform=openmm_platform)

    constraints = d.get_constraints()

    # this should be the 6 C-H bonds, as the central C-C to C-H bond is not constrained
    assert len(constraints) == 6

    d = mols[0].dynamics(
        constraint="h-bonds-not-heavy-perturbed", platform=openmm_platform
    )

    constraints = d.get_constraints()

    # this should be the 6 C-H bonds, plus the C-C to C-H bond, because this is
    # a non-heavy atom in one of the end states
    assert len(constraints) == 7

    d = mols[0].dynamics(
        constraint="bonds-not-heavy-perturbed", platform=openmm_platform
    )

    constraints = d.get_constraints()

    # this should be the 6 C-H bonds, plus the C-C to C-H bond, because this is
    # a non-heavy atom in one of the end states
    assert len(constraints) == 7

    # this was at lambda=0, so check that the constraint lengths are correct
    for constraint in constraints:
        dist = constraint[1].value()

        if dist == pytest.approx(1.0969):
            # this is a C-H bond
            assert (
                constraint[0][0].element().num_protons()
                + constraint[0][1].element().num_protons()
                == 7
            )
        elif dist == pytest.approx(1.5375):
            # this is a C-C bond
            assert (
                constraint[0][0].element().num_protons()
                + constraint[0][1].element().num_protons()
                == 12
            )
        else:
            assert False

    d = mols[0].dynamics(
        constraint="bonds-not-heavy-perturbed",
        lambda_value=1.0,
        platform=openmm_platform,
    )

    constraints = d.get_constraints()

    # this should be the 6 C-H bonds, plus the C-C to C-H bond, because this is
    # a non-heavy atom in one of the end states
    assert len(constraints) == 7

    # this was at lambda=1, so check that the constraint lengths are correct
    # (they should all be C-H bond lengths)
    for constraint in constraints:
        assert constraint[1].value() == pytest.approx(1.0969)


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

    d = mols[0].dynamics(constraint="h-bonds", platform=openmm_platform)

    constraints = d.get_constraints()

    # there are 6 bonds involving hydrogen
    assert len(constraints) == 6

    for constraint in constraints:
        # check for a single H atom
        assert constraint[0].atom("element H").element().symbol() == "H"

    d = mols[0].dynamics(
        constraint="bonds",
        platform=openmm_platform,
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


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_auto_constraints(ala_mols, openmm_platform):
    mols = ala_mols
    mol = mols[0]

    NA = 6.02214076e23
    CONV = 0.695039

    periods = {}

    for bond in mol.bonds():
        mass0 = bond.atom0().mass().value() / (1000.0 * NA)
        mass1 = bond.atom1().mass().value() / (1000.0 * NA)
        k = sr.mm.AmberBond(bond.potential(), sr.cas.Symbol("r")).k() * CONV

        mu = (mass0 * mass1) / (mass0 + mass1)

        # period in fs
        period = 1e15 * 2.0 * 3.14159 * (mu / k) ** 0.5
        periods[bond.id()] = period

    for factor in [5.0, 10.0, 15.0]:
        for timestep in [sr.u("1fs"), sr.u("2fs"), sr.u("4fs")]:
            fs = timestep.to("fs")

            constrained = []

            for bond in mol.bonds():
                period = periods[bond.id()]

                if period < factor * fs:
                    constrained.append(bond.id())

            d = mol.dynamics(
                constraint="auto-bonds",
                platform=openmm_platform,
                timestep=timestep,
                temperature="25oC",
                map={"auto_constraints_factor": factor},
            )

            constraints = d.get_constraints()

            assert len(constraints) == len(constrained)

            for constraint in constraints:
                bond = sr.bondid(
                    constraint[0].atom(0).index(), constraint[0].atom(1).index()
                )
                assert bond in constrained


@pytest.mark.xfail(reason="Unresolved bug.")
def test_asymmetric_constraints():
    # This test is for debugging a peculiar issue with one of the perturbations
    # from the MCL1 test suite. Here there are no ghost atoms and a single atom
    # changes type during the perturbation, from H to Cl. The constraints are
    # different for the two end states. Currently, the minimised energy at
    # lambda=1 does not match the minimised energy at lambda=0 when the end
    # states are swapped. From debugging, it seems that this is the caused by
    # calling context.applyConstraints() for the final constraint projection
    # following succesful minimisation. It's not clear if the bug lies in Sire,
    # or OpenMM.

    from math import isclose

    # Load the MCL1 perturbation. (Perturbable ligand is the last molecule.)
    mol = sr.load_test_files("mcl1_60_61.s3")[-1]

    # Create dynamics objects for the forward and backward perturbations.
    d_forwards = mol.dynamics(
        perturbable_constraint="h_bonds_not_heavy_perturbed",
        dynamic_constraints=True,
        include_constrained_energies=False,
    )
    d_backwards = mol.dynamics(
        perturbable_constraint="h_bonds_not_heavy_perturbed",
        include_constrained_energies=False,
        dynamic_constraints=True,
        swap_end_states=True,
    )

    # Set lambda so the dynamics states are equivalent.
    d_forwards.set_lambda(1.0, update_constraints=True)
    d_backwards.set_lambda(0.0, update_constraints=True)

    # Get the initial potential energies.
    nrg_forwards = d_forwards.current_potential_energy().value()
    nrg_backwards = d_backwards.current_potential_energy().value()

    # Check the potential energies are the same.
    assert isclose(nrg_forwards, nrg_backwards, rel_tol=1e-5)

    # Minimise both dynamics objects.
    d_forwards.minimise()
    d_backwards.minimise()

    # Get the minimisation logs.
    log_forwards = d_forwards._d.get_minimisation_log()
    log_backwards = d_backwards._d.get_minimisation_log()

    lines_forward = log_forwards.split("\n")
    for line in lines_forward:
        if "Final energy" in line:
            nrg_forwards = float(line.split()[2])

    lines_backward = log_backwards.split("\n")
    for line in lines_backward:
        if "Final energy" in line:
            nrg_backwards = float(line.split()[2])

    # Check the final energies from the logs are the same.
    assert isclose(nrg_forwards, nrg_backwards, rel_tol=1e-3)

    # Now get the final potential energies. (Post constraint projection.)
    nrg_forwards = d_forwards.current_potential_energy().value()
    nrg_backwards = d_backwards.current_potential_energy().value()

    # Check the minimised potential energies are the same.
    assert isclose(nrg_forwards, nrg_backwards, rel_tol=1e-3)
