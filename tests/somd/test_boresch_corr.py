"""Test the Boresch analytical and numerical correction methods."""

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


def _write_config(d):
    # Write the SOMD configuration file.
    filename = d.join("example_sim.cfg")
    with open(filename, "w") as f:
        f.write(
            "morphfile = ../../input/MORPH.vanish.pert\n"
            "topfile= ../../input/SYSTEM.top\n"
            "crdfile= ../../input/SYSTEM.crd\n"
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
            "use boresch restraints = True\n"
            'boresch restraints dictionary = {"anchor_points":{"r1":4946, "r2":4944, "r3":4949, "l1":11, "l2":2, "l3":3},"equilibrium_values":{"r0":5.92, "thetaA0":1.85, "thetaB0":1.59,"phiA0":-0.30, "phiB0":-1.55, "phiC0":2.90},"force_constants":{"kr":25.49, "kthetaA":66.74, "kthetaB":38.39, "kphiA":215.36, "kphiB":49.23, "kphiC":49.79}}\n'
        )

    return filename


@pytest.mark.skipif(sys.platform == "win32", reason="Not supported on Windows")
@pytest.mark.skipif(not have_scipy, reason="Not supported if scipy not installed")
def test_boresch_analytical_correction(tmpdir):
    """
    Check that the Boresch analytical correction yields
    the expected value.
    """
    d = tmpdir.mkdir("test_boresch_analytical_correction")

    cfg_path = _write_config(d)

    cmd = f"boresch_analytical_correction -C {cfg_path}"
    print(cmd)

    # Get the output of the command.
    proc = subprocess.run(shlex.split(cmd), capture_output=True)
    output = proc.stdout.decode("utf-8")

    print(output)

    # Check that the correction is correct
    data = output.split("Analytical correction for releasing Boresch restraints =")

    assert len(data) > 1

    dg = float(data[1].strip().split(" ")[0])

    assert dg == pytest.approx(-10.98, 0.01)


@pytest.mark.skipif(sys.platform == "win32", reason="Not supported on Windows")
@pytest.mark.skipif(not have_scipy, reason="Not supported if scipy not installed")
def test_boresch_numerical_correction(tmpdir):
    """
    Check that the Boresch numerical correction yields
    the expected value.
    """
    d = tmpdir.mkdir("test_boresch_analytical_correction")

    cfg_path = _write_config(d)

    cmd = f"boresch_numerical_correction -C {cfg_path}"

    print(cmd)

    # Get the output of the command.
    proc = subprocess.run(shlex.split(cmd), capture_output=True)
    output = proc.stdout.decode("utf-8")

    print(output)

    data = output.split("Numerical correction for releasing Boresch restraints =")

    assert len(data) > 1

    # Check that the correction is correct
    dg = float(data[1].strip().split(" ")[0])

    assert dg == pytest.approx(-10.98, 0.01)
