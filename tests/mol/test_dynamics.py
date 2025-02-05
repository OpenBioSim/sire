import sire as sr
import pytest


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_dynamics_return_type(ala_mols):
    mols = ala_mols

    m = mols.dynamics(platform="Reference").run("1 fs").commit()

    assert isinstance(m, type(mols))

    mol = mols[0]

    m = mol.dynamics(platform="Reference").run("1 fs").commit()

    assert isinstance(m, type(mol))

    atom = mol[0]

    a = atom.dynamics(platform="Reference").run("1 fs").commit()

    assert isinstance(a, type(atom))


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_minimisation_return_type(ala_mols):
    mols = ala_mols

    m = mols.minimisation(platform="Reference").run(1).commit()

    assert isinstance(m, type(mols))

    mol = mols[0]

    m = mol.minimisation(platform="Reference").run(1).commit()

    assert isinstance(m, type(mol))

    atom = mol[0]

    a = atom.minimisation(platform="Reference").run(1).commit()

    assert isinstance(a, type(atom))


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_cutoff_options(ala_mols):
    mols = ala_mols

    d = mols.dynamics(platform="Reference", cutoff="10 A")

    assert d.info().cutoff() == sr.u("10 A")

    d = mols[0].dynamics(platform="Reference", cutoff="infinite", vacuum=True)

    # OpenMM cannot have no cutoff, so sets it to a very large number
    assert d.info().cutoff() > sr.u("1000 A")

    d = mols[0].dynamics(platform="Reference", cutoff="none", vacuum=True)

    # OpenMM cannot have no cutoff, so sets it to a very large number
    assert d.info().cutoff() > sr.u("1000 A")

    d = mols.dynamics(platform="Reference", cutoff=sr.u("7.5A"), cutoff_type="PME")

    assert d.info().cutoff() == sr.u("7.5 A")
    assert d.info().cutoff_type() == "PME"


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_sample_frequency(ala_mols):
    """
    Test that energies and frames are saved at the correct frequency.
    """

    from sire.base import ProgressBar

    ProgressBar.set_silent()

    mols = ala_mols

    d = mols.dynamics(platform="Reference", timestep="1 fs")

    # Create a list of lambda windows.
    lambdas = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

    # Run 10 cycles of dynamics, saving energies every 2 fs and frames every 5 fs.
    for i in range(10):
        d.run(
            "1 fs",
            energy_frequency="2 fs",
            frame_frequency="5 fs",
            lambda_windows=lambdas,
        )

    # Get the energy trajectory.
    nrg_traj = d.energy_trajectory()

    # Make sure the trajectory has 5 frames.
    assert len(nrg_traj) == 5

    # Get the updated system.
    new_mols = d.commit()

    # Check that the trajectory has 2 frames.
    assert new_mols.num_frames() == 2

    # Now check when we request that a trajectory frame is saved when run exits.

    # Recreate the dynamics object.
    d = mols.dynamics(platform="Reference", timestep="1 fs")

    # Run 10 cycles of dynamics, saving energies frames on exit.
    for i in range(10):
        d.run(
            "1 fs",
            energy_frequency="2 fs",
            frame_frequency="5 fs",
            lambda_windows=lambdas,
            save_frame_on_exit=True,
            save_energy_on_exit=True,
        )

    # Get the energy trajectory.
    nrg_traj = d.energy_trajectory()

    # Make sure the trajectory has 10 frames.
    assert len(nrg_traj) == 10

    # Get the updated system.
    new_mols = d.commit()

    # Check that the trajectory has 10 frames.
    assert new_mols.num_frames() == 10
