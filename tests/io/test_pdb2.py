import pytest
import sire as sr
import tempfile


def test_pdb_space():
    """
    Test that we can parse CRYST1 records in a PDB file and that they
    are preserved when writing the file back out.
    """

    # Create a temporary working directory.
    tmp_dir = tempfile.TemporaryDirectory()
    tmp_path = tmp_dir.name
    pdb_file = f"{tmp_path}/test.pdb"

    # Load the test system.
    mols = sr.load(
        sr.expand(
            sr.tutorial_url,
            "cresset_triclinic_box.prm7",
            "cresset_triclinic_box.rst7",
        ),
        directory=tmp_path,
    )

    # Get the space from the system.
    space = mols.space()

    # Write back to the temporary path.
    sr.save(mols, pdb_file)

    # Load the PDB file.
    pdb_mols = sr.load(pdb_file)

    # Get the space from the PDB molecules.
    pdb_space = pdb_mols.space()

    # Make sure the spaces are the same to the precision of the PDB format.
    # Note that although the PDB uses 3dp for box dimensions, angles are only
    # written to 2dp.

    vec0 = space.vector0()
    vec1 = space.vector1()
    vec2 = space.vector2()
    pdb_vec0 = pdb_space.vector0()
    pdb_vec1 = pdb_space.vector1()
    pdb_vec2 = pdb_space.vector2()

    # X, Y, Z vector components should be the same.
    for i in range(3):
        assert vec0[i].value() == pytest.approx(pdb_vec0[i].value(), abs=1e-2)
        assert vec1[i].value() == pytest.approx(pdb_vec1[i].value(), abs=1e-2)
        assert vec2[i].value() == pytest.approx(pdb_vec2[i].value(), abs=1e-2)
