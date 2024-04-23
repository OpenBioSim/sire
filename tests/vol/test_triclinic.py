import sire as sr

import pytest
import tempfile

from openmm.app import AmberInpcrdFile, Topology


@pytest.fixture
def box():
    tri_box = sr.vol.TriclinicBox(
        1,
        1,
        1,
        120.0000229 * sr.units.degrees,
        120.0000229 * sr.units.degrees,
        90.0000000 * sr.units.degrees,
    )
    # Rotate the box using the minimimum precision from supported molecular
    # input files.
    tri_box.rotate(1e-5)
    # Reduce the lattice vectors, biasing towards left-tilting boxes.
    tri_box.reduce(-1e-8)

    return tri_box


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
        # Perform a stable rotation and reduction.
        box.rotate(1e-5)
        box.reduce(-1e-8)

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
        # Perform a stable rotation and reduction.
        box.rotate(1e-5)
        box.reduce(-1e-8)

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
            sr.tutorial_url,
            "cresset_triclinic_box.prm7",
            "cresset_triclinic_box.rst7",
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


def test_stream():
    """
    Test that a rotated and reduced TriclinicBox is correctly recovered when
    streaming to and from a file.
    """

    # Create a temporary working directory.
    tmp_dir = tempfile.TemporaryDirectory()
    tmp_path = tmp_dir.name
    s3_file = f"{tmp_path}/box.s3"

    # Create a triclinic box.
    box = sr.vol.TriclinicBox.truncated_octahedron(1.0, True, True)

    # Make sure the box has been rotated and reduced.
    assert box.is_rotated()
    assert box.is_reduced()

    # Stream to file.
    sr.stream.save(box, s3_file)

    # Load the file.
    recovered_box = sr.stream.load(s3_file)

    # Make sure the boxes are the same.
    assert recovered_box == box


def test_max_cutoff(ala_mols):
    """
    Test that the maximum cutoff is set correctly.
    """

    # Create a local copy of the molecules.
    mols = ala_mols.clone()

    # Create a cubic triclinic space.

    # Set the vectors.
    v0 = sr.maths.Vector(50, 0, 0)
    v1 = sr.maths.Vector(0, 50, 0)
    v2 = sr.maths.Vector(0, 0, 50)

    # Create the space.
    space = sr.vol.TriclinicBox(v0, v1, v2)

    # Check the maximum cutoff.
    assert space.maximum_cutoff() == 25 * sr.units.angstroms

    # Now set the space property on the molecules.
    mols.set_property("space", space)

    # Create a ForceFieldInfo object.
    ffinfo = sr.system.ForceFieldInfo(mols)

    # Check the cutoff. This is the maximum cutoff minus 1 angstrom.
    assert ffinfo.cutoff() == 24 * sr.units.angstroms
