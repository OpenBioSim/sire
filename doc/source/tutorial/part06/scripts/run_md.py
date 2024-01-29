import sire as sr

timestep = "4fs"
energy_frequency = "0.4ps"

constraint = "h-bonds-not-perturbed"
perturbable_constraint = "h-bonds-not-perturbed"

cutoff_type = "RF"

equil_time = "2ps"
run_time = "25ps"

lambda_values = [x / 100.0 for x in range(0, 101, 5)]

mols = sr.load(sr.expand(sr.tutorial_url, "merged_molecule.s3"))

for mol in mols.molecules("molecule property is_perturbable"):
    mol = sr.morph.link_to_reference(mol)
    mol = sr.morph.repartition_hydrogen_masses(mol, mass_factor=1.5)
    mols.update(mol)

mols = mols.minimisation(cutoff_type=cutoff_type).run().commit()

print("\nWater leg")

for i, lambda_value in enumerate(lambda_values):
    print(f"Simulating lambda={lambda_value:.2f}")
    # minimise the system at this lambda value
    min_mols = (
        mols.minimisation(
            cutoff_type=cutoff_type,
            lambda_value=lambda_value,
            constraint=constraint,
            perturbable_constraint="none",
        )
        .run()
        .commit()
    )

    # create a dynamics object for the system
    d = min_mols.dynamics(
        timestep=timestep,
        temperature="25oC",
        cutoff_type=cutoff_type,
        lambda_value=lambda_value,
        constraint=constraint,
        perturbable_constraint=perturbable_constraint,
    )

    # generate random velocities
    d.randomise_velocities()

    # equilibrate, not saving anything
    d.run(equil_time, save_frequency=0)
    print("Equilibration complete")
    print(d)

    # get the values of lambda for neighbouring windows
    lambda_windows = lambda_values[max(i - 1, 0) : min(len(lambda_values), i + 2)]

    # run the dynamics, saving the energy every 0.1 ps
    d.run(
        run_time,
        energy_frequency=energy_frequency,
        frame_frequency=0,
        lambda_windows=lambda_windows,
    )
    print("Dynamics complete")
    print(d)

    # stream the EnergyTrajectory to a sire save stream object
    sr.stream.save(
        d.commit().energy_trajectory(), f"energy_water_{lambda_value:.2f}.s3"
    )

print("\nVacuum leg")

for i, lambda_value in enumerate(lambda_values):
    print(f"Simulating lambda={lambda_value:.2f}")
    # minimise the system at this lambda value using the
    # same constraints as the dynamics (but switching
    # off the perturbable constraint)
    min_mols = (
        mols[0]
        .minimisation(
            lambda_value=lambda_value,
            vacuum=True,
            constraint=constraint,
            perturbable_constraint="none",
        )
        .run()
        .commit(return_as_system=True)
    )

    # create a dynamics object for the system
    d = min_mols.dynamics(
        timestep=timestep,
        temperature="25oC",
        lambda_value=lambda_value,
        constraint=constraint,
        perturbable_constraint=perturbable_constraint,
        vacuum=True,
    )

    # generate random velocities
    d.randomise_velocities()

    # equilibrate, not saving anything
    d.run(equil_time, save_frequency=0)
    print("Equilibration complete")
    print(d)

    # get the values of lambda for neighbouring windows
    lambda_windows = lambda_values[max(i - 1, 0) : min(len(lambda_values), i + 2)]

    # run the dynamics, saving the energy every 0.1 ps
    d.run(
        run_time,
        energy_frequency=energy_frequency,
        frame_frequency=0,
        lambda_windows=lambda_windows,
    )
    print("Dynamics complete")
    print(d)

    # stream the EnergyTrajectory to a sire save stream object
    sr.stream.save(
        d.commit().energy_trajectory(),
        f"energy_vacuum_{lambda_value:.2f}.s3",
    )

from alchemlyb.estimators import BAR

df = sr.morph.to_alchemlyb("energy_water_*.s3")

b = BAR()
b.fit(df)

dG_solv = b.delta_f_.loc[0.00, 1.00]

df = sr.morph.to_alchemlyb("energy_vacuum_*.s3")

b = BAR()
b.fit(df)

dG_vac = b.delta_f_.loc[0.00, 1.00]

print(f"Water phase:  {dG_solv} kcal mol-1")
print(f"Vacuum phase: {dG_vac} kcal mol-1")
print(f"Hydration free energy: {dG_solv - dG_vac} kcal mol-1")
