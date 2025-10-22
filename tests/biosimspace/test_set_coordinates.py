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


def test_set_coordinates_perturbable(merged_ethane_methanol):

    # Clone the input molecules.
    mols = merged_ethane_methanol.clone()

    # First link to the reference state.
    mols = sr.morph.link_to_reference(mols)

    # Store the existing coordinates as a NumPy array.
    coords = sr.io.get_coords_array(mols)

    # Modify the system to multiply all coordinates by 2.
    new_mols = sr.legacy.IO.setCoordinates(mols._system, (coords * 2.0).tolist())

    # Link to the reference state.
    new_mols = sr.system.System(new_mols)
    new_mols = sr.morph.link_to_reference(new_mols)

    # Get the new coordinates as a NumPy array.
    new_coords = sr.io.get_coords_array(new_mols)

    # Divide the new coordinates by the old coordinates.
    ratio = new_coords / coords

    # Set any NaN values to 2.0 (in case of zero coordinates).
    ratio = np.where(np.isnan(ratio), 2.0, ratio)

    # Make sure the new coordinates are as expected.
    assert (
        np.sum(np.round(ratio)) == 6.0 * mols.num_atoms()
    ), "Coordinates were not set correctly."

    # Now link to the perturbable state.
    mols = sr.morph.link_to_perturbed(mols)

    # Store the existing coordinates as a NumPy array.
    coords = sr.io.get_coords_array(mols)

    # Modify the system to multiply all coordinates by 2.
    new_mols = sr.legacy.IO.setCoordinates(
        mols._system, (coords * 2.0).tolist(), is_lambda1=True
    )

    # Link to the perturbable state.
    new_mols = sr.system.System(new_mols)
    new_mols = sr.morph.link_to_perturbed(new_mols)

    # Get the new coordinates as a NumPy array.
    new_coords = sr.io.get_coords_array(new_mols)

    # Divide the new coordinates by the old coordinates.
    ratio = new_coords / coords

    # Set any NaN values to 2.0 (in case of zero coordinates).
    ratio = np.where(np.isnan(ratio), 2.0, ratio)

    # Make sure the new coordinates are as expected.
    assert (
        np.sum(np.round(ratio)) == 6.0 * mols.num_atoms()
    ), "Coordinates were not set correctly."
