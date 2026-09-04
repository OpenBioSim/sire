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
def test_sample_frequency(ala_mols, openmm_platform):
    """
    Test that energies and frames are saved at the correct frequency.
    """

    from sire.base import ProgressBar

    ProgressBar.set_silent()

    mols = ala_mols

    d = mols.dynamics(platform=openmm_platform, timestep="1 fs")

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


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_crash_report(merged_ethane_methanol, openmm_platform):
    """
    Test that energies and frames are saved at the correct frequency.
    """

    import os
    import glob
    import tempfile
    from sire.base import ProgressBar

    ProgressBar.set_silent()

    mols = merged_ethane_methanol.clone()
    mols = sr.morph.link_to_reference(mols)

    d = mols.dynamics(platform=openmm_platform)

    # Run a short simulation within a temporary directory.
    tmpdir = tempfile.TemporaryDirectory()

    # Save the current directory.
    old_dir = os.getcwd()

    try:
        # Change to the temporary directory.
        os.chdir(tmpdir.name)

        # Run a short simulation, forcing a crash.
        d.run("1ps", save_crash_report=True)

        # Glob for the crash report files.
        crash_log = glob.glob("crash_*.log")
        crash_system = glob.glob("system_*.xml")
        crash_positions = glob.glob("positions_*.txt")

        # Make sure we have one of each file.
        assert len(crash_log) == 1
        assert len(crash_system) == 1
        assert len(crash_positions) == 1
    except:
        # Ingore exceptions raised during the dynamics run.
        pass
    finally:
        # Change back to the old directory.
        os.chdir(old_dir)


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_clock_and_energy_trajectory_swap(ala_mols):
    """
    Test that a single dynamics object can propagate several independent
    trajectories by swapping the clock and energy trajectory between blocks.

    This underpins replica exchange runs that re-use a bounded number of
    OpenMM contexts across a larger number of replicas.

    The NVE ensemble is used on the Reference platform so that the integrator
    is deterministic and has no random number stream. A thermostat would make
    the comparison meaningless, since the cached run interleaves the replicas
    through a single integrator and so consumes its RNG stream in a different
    order to separate dynamics objects.
    """

    from sire.base import ProgressBar

    ProgressBar.set_silent()

    mols = ala_mols.clone()
    mols.delete_all_frames()

    num_cycles = 5
    num_replicas = 2

    # No temperature, so this is NVE and the integrator has no RNG.
    kwargs = dict(platform="Reference", timestep="1 fs")
    run_kwargs = dict(
        energy_frequency="2 fs",
        frame_frequency="2 fs",
        lambda_windows=[0.0, 0.5, 1.0],
    )

    def potentials(traj):
        return [round(float(v), 8) for v in traj.to_pandas()["potential"]]

    # Build distinct starting states, so that the replicas follow genuinely
    # different trajectories and the test is not vacuous.
    seed = mols.dynamics(**kwargs)
    start_states = []
    for r in range(num_replicas):
        if r > 0:
            seed.run("10 fs", energy_frequency=0, frame_frequency=0)
        start_states.append(
            seed.context().getState(getPositions=True, getVelocities=True)
        )

    # Reference: one dynamics object per replica.
    ref = []
    for r in range(num_replicas):
        d = mols.dynamics(**kwargs)
        d.context().setState(start_states[r])
        d._d._clear_state()
        ref.append(d)

    for i in range(num_cycles):
        for d in ref:
            d.run("2 fs", **run_kwargs)

    ref_nrgs = [potentials(d.energy_trajectory()) for d in ref]
    ref_steps = [d.current_step() for d in ref]

    # Cached: a single dynamics object, with a clock and energy trajectory
    # per replica.
    slot = mols.dynamics(**kwargs)

    # Seed the per-replica trajectories from the slot's own, so that the
    # "ensemble" property is carried over.
    trajs = [slot._d.energy_trajectory() for _ in range(num_replicas)]
    assert all(len(t) == 0 for t in trajs)

    clocks = [slot._get_clock() for _ in range(num_replicas)]
    states = list(start_states)

    for i in range(num_cycles):
        for r in range(num_replicas):
            slot.context().setState(states[r])
            slot._set_clock(clocks[r])
            slot.set_energy_trajectory(trajs[r])
            slot._d._clear_state()

            slot.run("2 fs", **run_kwargs)

            clocks[r] = slot._get_clock()
            states[r] = slot.context().getState(getPositions=True, getVelocities=True)
            slot._d._sire_mols.delete_all_frames()

    cache_nrgs = [potentials(t) for t in trajs]
    cache_steps = [c["_current_step"] for c in clocks]

    # The replicas must be distinct, otherwise nothing is being tested.
    assert ref_nrgs[0] != ref_nrgs[1]

    # Each replica must have accumulated its own energies, and the clock must
    # have advanced as if it had a dynamics object to itself.
    assert cache_steps == ref_steps
    for r in range(num_replicas):
        assert len(cache_nrgs[r]) == num_cycles
        assert cache_nrgs[r] == ref_nrgs[r]


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_energy_cache_cleared_after_dynamics(ala_mols):
    """
    The context's energy cache must be invalidated after every dynamics block,
    not just one that saved energies. The integrator advances the positions
    without going through setPositions(), so nothing else clears it.
    """
    import openmm

    mols = ala_mols

    d = mols.dynamics(timestep="1fs", temperature="300K", platform="Reference")

    def direct():
        return (
            d.context()
            .getState(getEnergy=True)
            .getPotentialEnergy()
            .value_in_unit(openmm.unit.kilocalorie_per_mole)
        )

    assert d.current_potential_energy().value() == pytest.approx(direct())

    # A block that doesn't record an energy.
    d.run("50fs")
    assert d.current_potential_energy().value() == pytest.approx(direct())

    # A block that does.
    d.run("50fs", energy_frequency="10fs")
    assert d.current_potential_energy().value() == pytest.approx(direct())

    # A block that doesn't, again, now that a trajectory exists.
    d.run("50fs")
    assert d.current_potential_energy().value() == pytest.approx(direct())
