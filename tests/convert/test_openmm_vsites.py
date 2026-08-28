from math import isclose

import sire as sr
import pytest


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_vsite_params(ethane_12dichloroethane, openmm_platform):
    # Can we create an openmm system with the correct vsite parameters and charges
    mols = ethane_12dichloroethane

    # Just dichloroethane
    mol = mols[0].property("molecule1")

    # Set vsite properties
    vsite_dict = {
        "0": {
            "vs_indices": [0, 1, 2],
            "vs_ows": [1, 0, 0],
            "vs_xs": [1, -1, 0],
            "vs_ys": [0, 1, -1],
            "vs_local": [0.03, 0, 0],
        },
        "1": {
            "vs_indices": [3, 2, 1],
            "vs_ows": [1, 0, 0],
            "vs_xs": [1, -1, 0],
            "vs_ys": [0, 1, -1],
            "vs_local": [0.03, 0, 0],
        },
    }

    parents_dict = {str(atom_i): [] for atom_i in range(mol.num_atoms())}
    for v, vs in enumerate(vsite_dict):
        parent = vsite_dict[vs]["vs_indices"][0]
        parents_dict[str(parent)].append(v)

    n_virtual_sites = len(vsite_dict)
    vs_charges = [0.2, 0.2]

    cursor = mol.cursor()
    cursor.set("n_virtual_sites", n_virtual_sites)
    cursor.set("vs_charges", vs_charges)
    cursor.set("virtual_sites", vsite_dict)
    cursor.set("parents", parents_dict)
    mol = cursor.commit()

    map = {"platform": openmm_platform}

    # Create openmm system
    omm = sr.convert.to(mol, "openmm", map=map)

    omm_system = omm.getSystem()

    # Check the right number of virtual sites are added
    assert omm_system.getNumParticles() == mol.num_atoms() + n_virtual_sites

    # Check the VS parameters and charges are correct

    from openmm import LocalCoordinatesSite, Vec3, unit

    nb_force = next(
        force for force in omm_system.getForces() if force.getName() == "NonbondedForce"
    )

    for vs_index in range(n_virtual_sites):
        omm_vs = omm_system.getVirtualSite(mol.num_atoms() + vs_index)
        assert isinstance(omm_vs, LocalCoordinatesSite)

        origin_weights = omm_vs.getOriginWeights()
        assert list(origin_weights) == vsite_dict[str(vs_index)]["vs_ows"]

        x_weights = omm_vs.getXWeights()
        assert list(x_weights) == vsite_dict[str(vs_index)]["vs_xs"]

        y_weights = omm_vs.getYWeights()
        assert list(y_weights) == vsite_dict[str(vs_index)]["vs_ys"]

        local_position = omm_vs.getLocalPosition()
        assert local_position.value_in_unit(unit.nanometer) == Vec3(
            *[x for x in vsite_dict[str(vs_index)]["vs_local"]]
        )

        nb_params = nb_force.getParticleParameters(mol.num_atoms() + vs_index)
        assert vs_charges[vs_index] == nb_params[0].value_in_unit(
            unit.elementary_charge
        )


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_vsite_pertubation(ethane_12dichloroethane, openmm_platform):
    # Are vsite parameters scaled correctly by lambda
    mols = ethane_12dichloroethane

    # Just dichloroethane
    mol = mols[0]

    # Set vsite properties
    vsite_dict = {
        "0": {
            "vs_indices": [0, 1, 2],
            "vs_ows": [1, 0, 0],
            "vs_xs": [1, -1, 0],
            "vs_ys": [0, 1, -1],
            "vs_local": [0.03, 0, 0],
        },
        "1": {
            "vs_indices": [3, 2, 1],
            "vs_ows": [1, 0, 0],
            "vs_xs": [1, -1, 0],
            "vs_ys": [0, 1, -1],
            "vs_local": [0.03, 0, 0],
        },
    }

    parents_dict = {str(atom_i): [] for atom_i in range(mol.num_atoms())}
    for v, vs in enumerate(vsite_dict):
        parent = vsite_dict[vs]["vs_indices"][0]
        parents_dict[str(parent)].append(v)

    n_virtual_sites = len(vsite_dict)
    vs_charges0 = [0.1, 0.1]
    vs_charges1 = [0.2, 0.2]

    cursor = mol.cursor()
    cursor.set("n_virtual_sites", n_virtual_sites)
    cursor.set("vs_charges0", vs_charges0)
    cursor.set("vs_charges1", vs_charges1)
    cursor.set("virtual_sites", vsite_dict)
    cursor.set("parents", parents_dict)
    mol = cursor.commit()

    mol = sr.morph.link_to_reference(mol)

    from openmm import unit

    for lam in [0.0, 0.5, 1.0]:
        d = mol.dynamics(lambda_value=lam, platform=openmm_platform)
        system = d.context().getSystem()
        nb_force = next(
            force for force in system.getForces() if force.getName() == "NonbondedForce"
        )
        for vs_index in range(n_virtual_sites):
            expected_charge = (1 - lam) * vs_charges0[vs_index] + lam * vs_charges1[
                vs_index
            ]
            nb_charge = nb_force.getParticleParameters(mol.num_atoms() + vs_index)[0]
            assert expected_charge == nb_charge.value_in_unit(unit.elementary_charge)


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
@pytest.mark.parametrize(
    ["rest2_selection", "rest2_vsites"],
    [
        # No selection, so the whole perturbable molecule is the REST2 region
        # and both virtual sites are scaled.
        (None, [0, 1]),
        # Only the parent atom of the first virtual site is in the REST2
        # region, so only that virtual site is scaled.
        ("atomidx 0,1,2", [0]),
    ],
)
def test_vsite_rest2(
    ethane_12dichloroethane, openmm_platform, rest2_selection, rest2_vsites
):
    # Do virtual sites inherit the REST2 flag of their parent atom?
    mols = ethane_12dichloroethane

    # Just dichloroethane
    mol = mols[0]

    # Set vsite properties. The parent atoms are 0 and 3, i.e. the first
    # index of each vs_indices list.
    vsite_dict = {
        "0": {
            "vs_indices": [0, 1, 2],
            "vs_ows": [1, 0, 0],
            "vs_xs": [1, -1, 0],
            "vs_ys": [0, 1, -1],
            "vs_local": [0.03, 0, 0],
        },
        "1": {
            "vs_indices": [3, 2, 1],
            "vs_ows": [1, 0, 0],
            "vs_xs": [1, -1, 0],
            "vs_ys": [0, 1, -1],
            "vs_local": [0.03, 0, 0],
        },
    }

    parents_dict = {str(atom_i): [] for atom_i in range(mol.num_atoms())}
    for v, vs in enumerate(vsite_dict):
        parent = vsite_dict[vs]["vs_indices"][0]
        parents_dict[str(parent)].append(v)

    n_virtual_sites = len(vsite_dict)
    vs_charges0 = [0.1, 0.1]
    vs_charges1 = [0.2, 0.2]

    cursor = mol.cursor()
    cursor.set("n_virtual_sites", n_virtual_sites)
    cursor.set("vs_charges0", vs_charges0)
    cursor.set("vs_charges1", vs_charges1)
    cursor.set("virtual_sites", vsite_dict)
    cursor.set("parents", parents_dict)
    mol = cursor.commit()

    mol = sr.morph.link_to_reference(mol)

    d = mol.dynamics(
        lambda_value=0.0,
        platform=openmm_platform,
        rest2_selection=rest2_selection,
    )

    nb_force = next(
        force
        for force in d.context().getSystem().getForces()
        if force.getName() == "NonbondedForce"
    )

    # Store the unscaled charges at the same lambda value, so that the
    # comparison isolates the REST2 scaling from the lambda lever.
    d.set_lambda(0.0, rest2_scale=1.0)
    charges_initial = [
        nb_force.getParticleParameters(i)[0]._value
        for i in range(nb_force.getNumParticles())
    ]

    # Update the REST2 scaling factor.
    d.set_lambda(0.0, rest2_scale=2.0)

    # Store the scaling factor.
    scale = 0.5

    for vs_index in range(n_virtual_sites):
        parent = vsite_dict[str(vs_index)]["vs_indices"][0]

        # The virtual sites are appended after the atoms of the molecule.
        for i in [parent, mol.num_atoms() + vs_index]:
            charge = nb_force.getParticleParameters(i)[0]._value

            if vs_index in rest2_vsites:
                assert isclose(charge, charges_initial[i] * scale**0.5)
            else:
                assert isclose(charge, charges_initial[i])


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_vsite_restraints(solvated_neopentane_methane, openmm_platform):
    # Are restraints added correctly to vsite systems
    mols = solvated_neopentane_methane

    mol0 = mols[0]

    # Arbitrary vsite definition, we just want to check that the restraints
    # are correctly mapped to the new system with vsites
    vsite_dict = {
        "0": {
            "vs_indices": [0, 1, 2],
            "vs_ows": [1, 0, 0],
            "vs_xs": [1, -1, 0],
            "vs_ys": [0, 1, -1],
            "vs_local": [0.03, 0, 0],
        }
    }

    parents_dict = {str(atom_i): [] for atom_i in range(mol0.num_atoms())}
    for v, vs in enumerate(vsite_dict):
        parent = vsite_dict[vs]["vs_indices"][0]
        parents_dict[str(parent)].append(v)

    n_virtual_sites = len(vsite_dict)
    vs_charges0 = [0.1, 0.1]
    vs_charges1 = [0.2, 0.2]

    cursor = mol0.cursor()
    cursor.set("n_virtual_sites", n_virtual_sites)
    cursor.set("vs_charges0", vs_charges0)
    cursor.set("vs_charges1", vs_charges1)
    cursor.set("virtual_sites", vsite_dict)
    cursor.set("parents", parents_dict)
    mol0 = cursor.commit()
    mols.update(mol0)

    mols = sr.morph.link_to_reference(mols)

    start_mol1 = mol0.num_atoms() + n_virtual_sites
    base_particles = sum(m.num_atoms() for m in mols) + n_virtual_sites

    # Applying restraints to the first N water atoms (not realistic but we
    # just want to check the mapping is correct)
    restr_atoms = mols[1:].atoms()

    restraint_cases = [
        (
            "bond",
            lambda system: sr.restraints.bond(
                system, atoms0=restr_atoms[0], atoms1=restr_atoms[1]
            ),
            "BondRestraintForce",
            lambda force: (
                force.getBondParameters(0)[:2] == [start_mol1 + 0, start_mol1 + 1]
            ),
        ),
        (
            "inverse_bond",
            lambda system: sr.restraints.inverse_bond(
                system, atoms0=restr_atoms[1], atoms1=restr_atoms[2]
            ),
            "InverseBondRestraintForce",
            lambda force: (
                force.getBondParameters(0)[:2] == [start_mol1 + 1, start_mol1 + 2]
            ),
        ),
        (
            "morse_potential",
            lambda system: sr.restraints.morse_potential(
                system,
                atoms0=restr_atoms[0],
                atoms1=restr_atoms[1],
                r0="1.5A",
                k="100 kcal mol-1 A-2",
                de="25 kcal mol-1",
                auto_parametrise=False,
            )[0],
            "MorsePotentialRestraintForce",
            lambda force: (
                force.getBondParameters(0)[:2] == [start_mol1 + 0, start_mol1 + 1]
            ),
        ),
        (
            "positional",
            lambda system: sr.restraints.positional(system, atoms=restr_atoms[0]),
            "PositionalRestraintForce",
            lambda force: (
                force.getBondParameters(0)[0] == start_mol1 + 0
                and force.getBondParameters(0)[1] >= base_particles
                or force.getBondParameters(0)[1] == start_mol1 + 0
                and force.getBondParameters(0)[0] >= base_particles
            ),
        ),
        (
            "angle",
            lambda system: sr.restraints.angle(system, atoms=restr_atoms[:3]),
            "AngleRestraintForce",
            lambda force: (
                force.getAngleParameters(0)[:3]
                == [start_mol1 + 0, start_mol1 + 1, start_mol1 + 2]
            ),
        ),
        (
            "dihedral",
            lambda system: sr.restraints.dihedral(system, atoms=restr_atoms[:4]),
            "TorsionRestraintForce",
            lambda force: (
                force.getTorsionParameters(0)[:4]
                == [start_mol1 + 0, start_mol1 + 1, start_mol1 + 2, start_mol1 + 3]
            ),
        ),
        (
            "boresch",
            lambda system: sr.restraints.boresch(
                system,
                receptor=restr_atoms[0:3],
                ligand=restr_atoms[3:6],
            ),
            "BoreschRestraintForce",
            lambda force: (
                list(force.getBondParameters(0)[0])
                == [
                    start_mol1 + 0,
                    start_mol1 + 1,
                    start_mol1 + 2,
                    start_mol1 + 3,
                    start_mol1 + 4,
                    start_mol1 + 5,
                ]
            ),
        ),
    ]

    for restraint_name, build_restraints, force_name, validate in restraint_cases:
        restraints = build_restraints(mols)
        d = mols.dynamics(restraints=restraints, platform=openmm_platform)
        system = d.context().getSystem()
        force = next(
            force for force in system.getForces() if force.getName() == force_name
        )

        assert validate(force), f"Incorrect atom mapping for {restraint_name}"
