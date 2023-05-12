import sire as sr

import pytest


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
    and angles when the triclinic space is instantiated from box dimension
    vectors.
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
