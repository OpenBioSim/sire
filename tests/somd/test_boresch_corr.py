"""Test the Boresch analytical and numerical correction methods."""

import shlex
import subprocess
import os
import pytest

def test_boresch_analytical_correction():
    """Check that the Boresch analytical correction yields
    the expected value.
    """
    cfg_path = os.path.join(os.path.abspath(os.path.dirname(__file__)), "example_sim.cfg")
    cmd = f"boresch_analytical_correction -C {cfg_path}" 
    # Get the output of the command.
    proc = subprocess.run(shlex.split(cmd), capture_output=True)
    output = proc.stdout.decode("utf-8")
    # Check that the correction is correct
    dg = float(output.split("Analytical correction for releasing Boresch restraints =")[1].strip().split(" ")[0])
    assert dg == pytest.approx(-10.98, 0.01)

def test_boresch_numerical_correction():
    """Check that the Boresch numerical correction yields
    the expected value.
    """
    cfg_path = os.path.join(os.path.abspath(os.path.dirname(__file__)), "example_sim.cfg")
    cmd = f"boresch_numerical_correction -C {cfg_path}" 
    # Get the output of the command.
    proc = subprocess.run(shlex.split(cmd), capture_output=True)
    output = proc.stdout.decode("utf-8")
    # Check that the correction is correct
    dg = float(output.split("Numerical correction for releasing Boresch restraints =")[1].strip().split(" ")[0])
    assert dg == pytest.approx(-10.98, 0.01)
