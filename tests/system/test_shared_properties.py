import pytest
import sire as sr


def test_shared_properties(ala_traj, ala_mols):
    mols = ala_traj.clone()
    extra_mol = ala_mols[0]

    # This test will fail if another test has accidentally changed
    # 'ala_traj', e.g. by changing the frame number. You can
    # test for this by uncommenting the below line. If this test
    # passes, then some test using `ala_traj` has forgotten to
    # use 'clone' to create their own copy.
    # mols.load_frame(0)

    first_space = mols.property("space")
    first_time = mols.property("time")

    assert mols.shared_properties()["space"] == mols.space()
    assert mols.shared_properties()["time"] == mols.time()

    assert mols.space() == mols[10].property("space")
    assert mols.time() == mols[100].property("time")

    mols.load_frame(10)

    assert mols.shared_properties()["space"] == mols.space()
    assert mols.shared_properties()["time"] == mols.time()

    assert mols.space() == mols[10].property("space")
    assert mols.time() == mols[100].property("time")

    mols.set_property("space", sr.vol.Cartesian())

    assert mols.space() == sr.vol.Cartesian()

    assert mols.shared_properties()["space"] == mols.space()
    assert mols.space() == mols[10].property("space")
    assert mols.space() == mols[100].property("space")

    mols.load_frame(5)

    assert mols.space() != sr.vol.Cartesian()

    assert mols.shared_properties()["space"] == mols.space()
    assert mols.shared_properties()["time"] == mols.time()

    assert mols.space() == mols[10].property("space")
    assert mols.time() == mols[100].property("time")

    assert extra_mol.property("space") != mols.space()
    assert extra_mol.property("time") != mols.time()

    assert mols[0].property("space") != first_space
    assert mols[0].property("time") != first_time

    mols.load_frame(0)

    assert mols.space() == first_space
    assert mols.time() == first_time

    assert mols[0].property("space") == first_space
    assert mols[0].property("time") == first_time
    assert mols[-1].property("space") == first_space
    assert mols[-1].property("time") == first_time

    mols.load_frame(10)

    mols.add(extra_mol)

    assert mols[-1].number() == extra_mol.number()

    assert mols[-1].property("space") == mols.space()
    assert mols[-1].property("time") == mols.time()

    mols.add_shared_property("cat", "meow")

    assert mols.shared_properties()["cat"] == "meow"

    assert mols.property("cat") == mols[0].property("cat")
    assert mols.property("cat") == mols[-1].property("cat")

    mols.remove_shared_property("cat")

    assert not mols.shared_properties().has_property("cat")
    assert not mols.has_property("cat")
    assert not mols[0].has_property("cat")
    assert not mols[-1].has_property("cat")

    mols.set_property("cat", "meow")

    assert mols.property("cat") == "meow"
    assert not mols.shared_properties().has_property("cat")
    assert not mols[0].has_property("cat")
    assert not mols[-1].has_property("cat")

    mols.remove_all_shared_properties()

    assert mols.shared_properties().is_empty()

    assert not mols.has_property("space")
    assert not mols.has_property("time")
    assert not mols[0].has_property("space")
    assert not mols[0].has_property("time")
    assert not mols[-1].has_property("space")
    assert not mols[-1].has_property("time")

    mols.load_frame(0)

    assert mols[0].property("space") == first_space
    assert mols[0].property("time") == first_time

    # these won't be set as they aren't shared any more
    assert not mols.has_property("space")
    assert not mols.has_property("time")
