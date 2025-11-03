import pytest
import sire as sr


@pytest.mark.parametrize("save_velocities", [True, False])
def test_velocities(tmpdir, triclinic_protein, save_velocities):
    mols = triclinic_protein

    # write out a AMBER files file

    d = tmpdir.mkdir("test_velocities")
    f = sr.save(
        mols, d.join("test"), format=["RST7", "PRM7"], save_velocities=save_velocities
    )

    # load the files back in
    check = sr.load(f[0], f[1])

    # check that the velocities are present if requested
    assert check[0].has_property("velocity") == save_velocities
