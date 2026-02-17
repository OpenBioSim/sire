import sire as sr
import pytest


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
@pytest.mark.parametrize("molecules", ["kigaki_mols", "merged_ethane_methanol"])
def test_openmm_positional_restraints(molecules, openmm_platform, request):
    mols = request.getfixturevalue(molecules)

    if mols[0].is_perturbable():
        mols = sr.morph.link_to_reference(mols)

    mol = mols[0]

    map = {"space": sr.vol.Cartesian(), "platform": openmm_platform}

    # test restraining all C atoms
    restraints = sr.restraints.positional(mol, atoms="element C")

    assert len(restraints) == len(mol.atoms("element C"))

    coords = []

    for restraint, atom in zip(restraints, mol.atoms("element C")):
        assert restraint.is_atom_restraint()
        assert restraint.position() == atom.coordinates()
        assert restraint.r0().is_zero()
        coords.append(atom.coordinates())

    d = mol.dynamics(restraints=restraints, timestep="4fs", map=map)

    d.run("1ps")

    mol = d.commit()

    for atom, coords in zip(mol.atoms("element C"), coords):
        assert (atom.coordinates() - coords).length() < 0.1 * sr.units.angstrom


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_openmm_distance_restraints(ala_mols, openmm_platform):
    mols = ala_mols

    mols = mols[0:2]

    map = {"space": sr.vol.Cartesian(), "platform": openmm_platform}

    # test restraining the distance between the first and last molecule
    restraints = sr.restraints.distance(
        mols, atoms0=mols[0][0], atoms1=mols[-1][0], r0="5A"
    )

    assert len(restraints) == 1

    # check that we get the same result using bond restraints
    restraints2 = sr.restraints.bond(
        mols, atoms0=mols[0][0], atoms1=mols[-1][0], r0="5A", use_pbc=True
    )

    assert len(restraints2) == 1

    assert restraints == restraints2

    assert restraints[0].is_atom_restraint()
    assert restraints[0].r0() == sr.u("5A")

    # need to minimise the energy to get the restraint to work
    mols = mols.minimisation(restraints=restraints, map=map)()

    d = mols.dynamics(restraints=restraints, timestep="4fs", map=map)

    d.run("10ps")

    mols = d.commit()

    new_coords = [mols[0][0].coordinates(), mols[-1][0].coordinates()]

    assert (new_coords[0] - new_coords[1]).length().value() == pytest.approx(5.0, 1e-2)


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_openmm_fixed_atoms(kigaki_mols, openmm_platform):
    mols = kigaki_mols

    mol = mols[0]

    map = {"space": sr.vol.Cartesian(), "platform": openmm_platform}

    # test fixing all C atoms
    coords = []

    for atom in mol.atoms("element C"):
        coords.append(atom.coordinates())

    d = mol.dynamics(fixed="element C", timestep="1fs", map=map)

    d.run("1ps")

    mol = d.commit()

    for atom, coords in zip(mol.atoms("element C"), coords):
        assert (atom.coordinates() - coords).length() < 0.001 * sr.units.angstrom


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_openmm_alchemical_restraints(ala_mols, openmm_platform):
    mols = ala_mols

    mol = mols[0]

    map = {"space": sr.vol.Cartesian(), "platform": openmm_platform}

    # test scaling a positional restraint
    restraints = sr.restraints.positional(mol, atoms="element C")

    # move the molecule so that the restraints have some energy
    mol = mol.move().translate(sr.maths.Vector(1, 1, 1)).commit()

    d = mol.dynamics(timestep="1fs", restraints=None, map=map)

    nrg_0 = d.current_potential_energy()

    d = mol.dynamics(timestep="1fs", restraints=restraints, map=map)

    nrg_1 = d.current_potential_energy()

    assert nrg_1 != nrg_0

    l = sr.cas.LambdaSchedule()

    l.add_stage("restraints", l.lam() * l.initial())
    l.set_equation(
        stage="restraints", lever="restraint", equation=l.lam() * l.initial()
    )

    d = mol.dynamics(timestep="1fs", restraints=restraints, schedule=l, map=map)

    d.set_lambda(0)

    assert d.current_potential_energy().value() == pytest.approx(nrg_0.value(), 1e-6)

    d.set_lambda(1)

    assert d.current_potential_energy().value() == pytest.approx(nrg_1.value(), 1e-6)

    d.set_lambda(0.3)

    assert d.current_potential_energy().value() == pytest.approx(
        nrg_0.value() + 0.3 * (nrg_1.value() - nrg_0.value()), 1e-6
    )


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_openmm_named_restraints(ala_mols, openmm_platform):
    mols = ala_mols

    mol = mols[0]

    map = {"space": sr.vol.Cartesian(), "platform": openmm_platform}

    # test using named restraints, that we can scale these independently
    posrests = sr.restraints.positional(mol, atoms="element C", name="positional")

    dstrests = sr.restraints.distance(
        mol, atoms0=mol[0], atoms1=mol[-1], name="distance", r0="5A"
    )

    restraints = [posrests, dstrests]

    # move the molecule so that the restraints have some energy
    mol = mol.move().translate(sr.maths.Vector(1, 1, 1)).commit()

    d = mol.dynamics(timestep="1fs", restraints=None, map=map)

    nrg_0 = d.current_potential_energy()

    d = mol.dynamics(timestep="1fs", restraints=posrests, map=map)

    nrg_1_posrests = d.current_potential_energy()

    d = mol.dynamics(timestep="1fs", restraints=dstrests, map=map)

    nrg_1_dstrests = d.current_potential_energy()

    d = mol.dynamics(timestep="1fs", restraints=restraints, map=map)

    nrg_1_1 = d.current_potential_energy()

    assert nrg_0 != nrg_1_1
    assert nrg_0 != nrg_1_dstrests
    assert nrg_0 != nrg_1_posrests
    assert nrg_1_dstrests != nrg_1_posrests

    l = sr.cas.LambdaSchedule()

    l.add_stage("1", 0)
    l.set_equation(stage="1", lever="positional", equation=l.lam() * l.initial())

    l.add_stage("2", 0)
    l.set_equation(stage="2", lever="distance", equation=l.lam() * l.initial())

    l.add_stage("3", 0)
    l.set_equation(stage="3", lever="positional", equation=l.lam() * l.initial())
    l.set_equation(stage="3", lever="distance", equation=l.lam() * l.initial())

    d = mol.dynamics(timestep="1fs", restraints=restraints, schedule=l, map=map)

    d.set_lambda(0)

    assert d.current_potential_energy().value() == pytest.approx(nrg_0.value(), 1e-6)

    d.set_lambda(1)

    assert d.current_potential_energy().value() == pytest.approx(nrg_1_1.value(), 1e-6)

    d.set_lambda(0.99999999999 / 3.0)

    assert d.current_potential_energy().value() == pytest.approx(
        nrg_1_posrests.value(), 1e-6
    )

    d.set_lambda(1.0 / 3.0)

    assert d.current_potential_energy().value() == pytest.approx(nrg_0.value(), 1e-6)

    d.set_lambda(1.99999999999 / 3.0)

    assert d.current_potential_energy().value() == pytest.approx(
        nrg_1_dstrests.value(), 1e-6
    )

    d.set_lambda(2.0 / 3.0)

    assert d.current_potential_energy().value() == pytest.approx(nrg_0.value(), 1e-6)


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
@pytest.mark.parametrize(
    "restraint,default_pbc",
    [
        ("distance", True),
        ("bond", False),
    ],
)
def test_openmm_restraints_pbc(ala_mols, restraint, default_pbc, openmm_platform):
    mols = ala_mols

    # get the restraint class
    restraint_class = getattr(sr.restraints, restraint)

    # create a distance restraint between the first and last atom
    restraints = restraint_class(
        mols,
        atoms0=mols[0][0],
        atoms1=mols[-1][0],
        r0="5A",
    )

    # make sure the _use_pbc flag is set to the default value
    assert restraints.uses_pbc() == default_pbc

    # create a dynamics object
    d = mols.dynamics(restraints=restraints, platform=openmm_platform)

    # find the restraint force
    for force in d.context().getSystem().getForces():
        if force.getName() == "BondRestraintForce":
            restraint_force = force
            break

    # check that the force is using periodic boundary conditions
    assert restraint_force.usesPeriodicBoundaryConditions() == default_pbc

    # now create a distance restraint with use_pbc set to the non-default value
    restraints = restraint_class(
        mols, atoms0=mols[0][0], atoms1=mols[-1][0], r0="5A", use_pbc=not default_pbc
    )

    # make sure the _use_pbc flag is set to False
    assert restraints.uses_pbc() == (not default_pbc)

    # create a dynamics object
    d = mols.dynamics(restraints=restraints, platform=openmm_platform)

    # find the restraint force
    for force in d.context().getSystem().getForces():
        if force.getName() == "BondRestraintForce":
            restraint_force = force
            break

    # check that the force is not using periodic boundary conditions
    assert restraint_force.usesPeriodicBoundaryConditions() == (not default_pbc)
