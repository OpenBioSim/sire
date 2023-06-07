import sire as sr

import pytest
import tempfile

from openmm.app import AmberInpcrdFile, Topology


@pytest.fixture
def box():
    return sr.vol.TriclinicBox(
        1,
        1,
        1,
        120.0000229 * sr.units.degrees,
        120.0000229 * sr.units.degrees,
        90.0000000 * sr.units.degrees,
    )


def test_reduction_vectors(box):
    """
    Test that repeat lattice reductions give consistent box dimensions
    and angles when the triclinic space is instantiated from box dimension
    vectors.
    """

    # Store the initial values.
    alpha = box.alpha()
    beta = box.beta()
    gamma = box.gamma()

    for i in range(0, 20):
        # Re-create object using the previous box vectors as input.
        box = sr.vol.TriclinicBox(box.vector0(), box.vector1(), box.vector2())

        # Make sure that the angles haven't changed.
        assert alpha == box.alpha()
        assert beta == box.beta()
        assert gamma == box.gamma()


def test_reduction_angles(box):
    """
    Test that repeat lattice reductions give consistent box dimensions
    and angles when the triclinic space is instantiated from box angles.
    """

    # Store the initial values.
    alpha = box.alpha()
    beta = box.beta()
    gamma = box.gamma()

    for i in range(0, 20):
        # Re-create object using the previous box angles as input.
        box = sr.vol.TriclinicBox(
            1,
            1,
            1,
            box.alpha() * sr.units.degrees,
            box.beta() * sr.units.degrees,
            box.gamma() * sr.units.degrees,
        )

        # Make sure that the angles haven't changed.
        assert alpha == box.alpha()
        assert beta == box.beta()
        assert gamma == box.gamma()


def test_cresset_box():
    """
    Test that OpenMM is happy with Cresset's default triclinic box vectors after
    application of our triclinic lattice reduction.
    """

    # Create a temporary working directory.
    tmp_dir = tempfile.TemporaryDirectory()
    tmp_path = tmp_dir.name
    rst7_file = f"{tmp_path}/test.rst7"

    # Load the test system.
    system = sr.load(
        sr.expand(
            sr.tutorial_url, "cresset_triclinic_box.prm7", "cresset_triclinic_box.rst7"
        ),
        directory=tmp_path,
    )

    # Write back to the temporary path.
    sr.save(system, rst7_file)

    # Load a minimal coordinate file containing the default Cresset triclinic box.
    inpcrd = AmberInpcrdFile(rst7_file)

    # Create an empty topology.
    topology = Topology()

    # Try to set the periodic box vectors using those from the inpcrd file.
    topology.setPeriodicBoxVectors(inpcrd.boxVectors)
