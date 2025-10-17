import sire as sr
import numpy as np


def test_set_coordinates(ala_mols):

    # Clone the input molecules.
    mols = ala_mols.clone()

    # Store the existing coordinates as a NumPy array.
    coords = sr.io.get_coords_array(mols)

    # Modify the system to multiply all coordinates by 2.
    new_mols = sr.legacy.IO.setCoordinates(mols._system, (coords * 2.0).tolist())

    # Get the new coordinates as a NumPy array.
    new_coords = sr.io.get_coords_array(new_mols)

    # Make sure the new coordinates are as expected.
    assert (
        np.sum(np.round(new_coords / coords)) == 6.0 * mols.num_atoms()
    ), "Coordinates were not set correctly."
