from math import isclose

import sire as sr


def test_rest2():
    """
    Test that REST2 modifications are correctly applied to the system.
    """

    # Load the test perturbation.
    mol = sr.load_test_files("toluene_methane.s3")

    # Link to the reference state.
    mol = sr.morph.link_to_reference(mol)

    # Work out the number of dihedrals in the system.
    num_dihedrals = len(mol.dihedrals())

    # Create a dynamics object.
    d = mol.dynamics()

    # Extract the OpenMM system.
    omm_system = d.context().getSystem()

    # Find the PeriodicTorsionForce.
    for force in omm_system.getForces():
        if force.getName() == "PeriodicTorsionForce":
            break

    # Store the initial parameters.
    torsion_params_initial = []
    for i in range(force.getNumTorsions()):
        torsion_params_initial.append(force.getTorsionParameters(i)[-1])

    # Find the NonbondedForce.
    for force in omm_system.getForces():
        if force.getName() == "NonbondedForce":
            break

    # Store the initial parameters.
    nonbonded_params_initial = []
    for i in range(force.getNumParticles()):
        charge, sigma, epsilon = force.getParticleParameters(i)
        nonbonded_params_initial.append((charge, epsilon))
    exception_params_initial = []
    for i in range(force.getNumExceptions()):
        exception_params_initial.append(force.getExceptionParameters(i)[-3::2])

    # Find the ghost/ghost nonbonded interaction.
    for force in omm_system.getForces():
        if force.getName() == "GhostGhostNonbondedForce":
            break

    # Store the initial parameters.
    ghost_ghost_params_initial = []
    for i in range(force.getNumParticles()):
        charge, half_sigma, two_sqrt_epsilon, alpha, kappa = (
            force.getParticleParameters(i)
        )
        ghost_ghost_params_initial.append((charge, two_sqrt_epsilon))

    # Find the ghost/non-ghost nonbonded interaction.
    for force in omm_system.getForces():
        if force.getName() == "GhostNonGhostNonbondedForce":
            break

    # Store the initial parameters.
    ghost_non_ghost_params_initial = []
    for i in range(force.getNumParticles()):
        charge, half_sigma, two_sqrt_epsilon, alpha, kappa = (
            force.getParticleParameters(i)
        )
        ghost_non_ghost_params_initial.append((charge, two_sqrt_epsilon))

    # Update the REST2 scaling factor.
    d.set_lambda(0.0, rest2_scale=2.0)

    # Extract the OpenMM system.
    omm_system = d.context().getSystem()

    # Find the PeriodicTorsionForce.
    for force in omm_system.getForces():
        if force.getName() == "PeriodicTorsionForce":
            break

    # Store the modified parameters.
    torsion_params_modified = []
    for i in range(force.getNumTorsions()):
        torsion_params_modified.append(force.getTorsionParameters(i)[-1])

    # Find the NonbondedForce.
    for force in omm_system.getForces():
        if force.getName() == "NonbondedForce":
            break

    # Store the modified parameters.
    nonbonded_params_modified = []
    for i in range(force.getNumParticles()):
        charge, sigma, epsilon = force.getParticleParameters(i)
        nonbonded_params_modified.append((charge, epsilon))
    exception_params_modified = []
    for i in range(force.getNumExceptions()):
        exception_params_modified.append(force.getExceptionParameters(i)[-3::2])

    # Find the ghost/ghost nonbonded interaction.
    for force in omm_system.getForces():
        if force.getName() == "GhostGhostNonbondedForce":
            break

    # Store the modified parameters.
    ghost_ghost_params_modified = []
    for i in range(force.getNumParticles()):
        charge, half_sigma, two_sqrt_epsilon, alpha, kappa = (
            force.getParticleParameters(i)
        )
        ghost_ghost_params_modified.append((charge, two_sqrt_epsilon))

    # Find the ghost/non-ghost nonbonded interaction.
    for force in omm_system.getForces():
        if force.getName() == "GhostNonGhostNonbondedForce":
            break

    # Store the modified parameters.
    ghost_non_ghost_params_modified = []
    for i in range(force.getNumParticles()):
        charge, half_sigma, two_sqrt_epsilon, alpha, kappa = (
            force.getParticleParameters(i)
        )
        ghost_non_ghost_params_modified.append((charge, two_sqrt_epsilon))

    # Store the scaling factor.
    scale = 0.5

    # All dihedral force constants should be halved.
    for i in range(num_dihedrals):
        assert isclose(
            torsion_params_modified[i]._value, torsion_params_initial[i]._value * scale
        )

    # All impropers should be unchanged.
    for i in range(num_dihedrals, len(torsion_params_initial)):
        assert isclose(
            torsion_params_modified[i]._value, torsion_params_initial[i]._value
        )

    # Nonbonded charages should be scaled by the square root of the scaling
    # factor and epsilon should be scaled by the scaling factor.
    for i in range(len(nonbonded_params_initial)):
        charge, epsilon = nonbonded_params_initial[i]
        charge_modified, epsilon_modified = nonbonded_params_modified[i]
        assert isclose(charge_modified._value, charge._value * scale**0.5)
        assert isclose(epsilon_modified._value, epsilon._value * scale)

    # For exceptions, both the charge and epsilon should be scaled by the
    # scaling factor.
    for i in range(len(exception_params_initial)):
        charge, epsilon = exception_params_initial[i]
        charge_modified, epsilon_modified = exception_params_modified[i]
        assert isclose(charge_modified._value, charge._value * scale)
        if epsilon._value > 1e-6:
            assert isclose(epsilon_modified._value, epsilon._value * scale)

    # Ghost/ghost nonbonded charges should be scaled by the square root of the
    # scaling factor and epsilon should be scaled by the scaling factor.
    # (Note that epsilon is stored as sqrt(epsilon) so the scaling factor is
    # also square rooted.)
    for i in range(len(ghost_ghost_params_initial)):
        charge, two_sqrt_epsilon = ghost_ghost_params_initial[i]
        charge_modified, two_sqrt_epsilon_modified = ghost_ghost_params_modified[i]
        assert isclose(charge_modified, charge * scale**0.5)
        if two_sqrt_epsilon > 1e-6:
            assert isclose(two_sqrt_epsilon_modified, two_sqrt_epsilon * scale**0.5)

    # Ghost/non-ghost nonbonded charges should be scaled by the square root of
    # the scaling factor and epsilon should be scaled by the scaling factor.
    # (Note that epsilon is stored as sqrt(epsilon) so the scaling factor is
    # also square rooted.)
    for i in range(len(ghost_non_ghost_params_initial)):
        charge, two_sqrt_epsilon = ghost_non_ghost_params_initial[i]
        charge_modified, two_sqrt_epsilon_modified = ghost_non_ghost_params_modified[i]
        assert isclose(charge_modified, charge * scale**0.5)
        if two_sqrt_epsilon > 1e-6:
            assert isclose(two_sqrt_epsilon_modified, two_sqrt_epsilon * scale**0.5)
