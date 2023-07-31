"""
Test the standardstatecorrection script, which is used to 
compute the standard state correction when using multiple 
distance restraints.
"""

import shlex
import subprocess
import os
import sys
import pytest

try:
    import scipy

    have_scipy = True
except ImportError:
    have_scipy = False


try:
    import mdtraj

    have_mdtraj = True
except ImportError:
    have_mdtraj = False


import sire as sr


try:
    _wget = sr.legacy.Base.findExe("wget")
    have_wget = True
except Exception:
    have_wget = False


@pytest.fixture(scope="session")
def input_tmpdir(tmpdir_factory):
    """
    Create a temporary directory with the required
    configuration, topology, coordinate and trajectory
    files.
    """
    input_dir = tmpdir_factory.mktemp("input")

    # Write the SOMD configuration file.
    with open(input_dir.join("sim.cfg"), "w") as f:
        f.write(
            "perturbed residue number = 1\n"
            "nmoves = 25000\n"
            "ncycles = 80\n"
            "buffered coordinates frequency = 5000\n"
            "save coordinates = True\n"
            "timestep = 4 * femtosecond\n"
            "constraint = allbonds\n"
            "hydrogen mass repartitioning factor = 4.0\n"
            "cutoff type = cutoffperiodic\n"
            "cutoff distance = 12*angstrom\n"
            "barostat = True\n"
            "andersen = True\n"
            "energy frequency = 100\n"
            "precision = mixed\n"
            "minimise = True\n"
            "equilibrate = False\n"
            "equilibration iterations = 5000\n"
            "center solute = True\n"
            "reaction field dielectric = 78.3\n"
            "minimal coordinate saving = False\n"
            "lambda array = 0.000, 0.025, 0.050, 0.075, 0.100, 0.125, 0.150, 0.175, 0.200, 0.225, 0.250, 0.275, 0.300, 0.325, 0.350, 0.375, 0.400, 0.425, 0.450, 0.475, 0.500, 0.525, 0.550, 0.575, 0.600, 0.625, 0.650, 0.675, 0.700, 0.725, 0.750, 0.800, 0.850, 0.900, 0.950, 1.000\n"
            "use distance restraints = True\n"
            "distance restraints dictionary = {(21, 4950): (2.72, 20, 0.5), (18, 961): (3.09, 20, 0.5), (17, 512): (3.25, 20, 0.5), (19, 512): (3.69, 20, 0.5)}\n"
        )

    # Download the input files.
    DCD = "https://github.com/OpenBioSim/sire_bigtests/raw/main/io/lj_cor/traj.dcd"
    TOP = "https://github.com/OpenBioSim/sire_bigtests/raw/main/io/receptor_ligand_restraints/SYSTEM.top"

    for input_file in [DCD, TOP]:
        subprocess.check_call(
            shlex.split(
                f"wget -O {input_dir.join(input_file.split('/')[-1])} {input_file}"
            )
        )

    return input_dir


@pytest.mark.skipif(sys.platform == "win32", reason="Not supported on Windows")
@pytest.mark.skipif(not have_scipy, reason="Cannot run as scipy not available")
@pytest.mark.skipif(not have_wget, reason="Cannot run as wget is not available")
@pytest.mark.skipif(not have_mdtraj, reason="Cannot run as mdtraj not available")
@pytest.mark.slow
def test_standard_state_correction(input_tmpdir):
    """
    Check that the standardstatecorrection yields the
    the expected value.
    """
    # Note that far too few integration points (in particular far too few orientations)
    # have been used in order to run the test quickly. This results in a ridiculously large
    # standard state correction, but the purpose of this test is to check that the script
    # runs correctly, not to check that the correction is correct.
    cmd = f"standardstatecorrection -C {input_tmpdir.join('sim.cfg')} -r {input_tmpdir.join('traj.dcd')} -t {input_tmpdir.join('SYSTEM.top')} -s 1 -b 1 -d 0.5 -o 4"
    print(cmd)

    # Get the output of the command
    output = subprocess.check_output(shlex.split(cmd)).decode("utf-8")
    print(output)

    # Check that the correction is correct
    data = output.split(
        "Free energy change upon removing the restraint and applying standard state conditions ="
    )

    print(data)
    assert len(data) > 1

    dg = float(data[1].strip().split(" ")[0])

    # Note that this correction is extremely large and incorrect - far too few integration points
    # have been used in order to run the test correctly.
    assert dg == pytest.approx(-119.84, 0.01)
