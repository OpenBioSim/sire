import pytest
import sire as sr


def _assert_correct_com(com):
    correct_com = sr.maths.Vector(31.07343482, 31.89108144, 0.87194654)

    assert com.x().value() == pytest.approx(correct_com.x().value(), 1e-3)

    assert com.y().value() == pytest.approx(correct_com.y().value(), 1e-3)

    assert com.z().value() == pytest.approx(correct_com.z().value(), 1e-3)


def test_make_whole(wrapped_mols):
    mols = wrapped_mols.clone()

    wrong_center = mols[0].evaluate().center_of_mass()

    mols.make_whole()

    right_center = mols[0].evaluate().center_of_mass()

    assert wrong_center != right_center

    _assert_correct_com(right_center)

    mols = wrapped_mols

    mol = mols[0]

    assert mol.evaluate().center_of_mass() == wrong_center

    mol = mol.move().make_whole().commit()

    assert mol.evaluate().center_of_mass() == right_center

    c = wrapped_mols.cursor()

    assert c.commit()[0].evaluate().center_of_mass() == wrong_center

    c.make_whole()

    mol = c.commit()[0].evaluate().center_of_mass() == right_center


def test_auto_make_whole_on_load_frame(wrapped_mols):
    mols = wrapped_mols.clone()

    wrong_center = mols[0].evaluate().center_of_mass()

    mols.load_frame(0, map={"make_whole": True})

    right_center = mols[0].evaluate().center_of_mass()

    assert right_center != wrong_center

    _assert_correct_com(right_center)


def test_auto_make_whole_on_load():
    mols = sr.load_test_files("wrapped.rst7", "wrapped.prm7", map={"make_whole": True})

    _assert_correct_com(mols[0].evaluate().center_of_mass())


def test_auto_make_whole_on_load_no_breakage(kigaki_mols):
    mols = sr.load_test_files(
        "kigaki.gro",
        "kigaki.top",
        map={"make_whole": True},
    )

    assert (
        kigaki_mols[0].evaluate().center_of_mass()
        == mols[0].evaluate().center_of_mass()
    )

    mols.load_frame(0, map={"make_whole": True})

    assert (
        kigaki_mols[0].evaluate().center_of_mass()
        == mols[0].evaluate().center_of_mass()
    )


def test_make_whole_center_args(ala_mols):
    mols = ala_mols

    c = mols[0].cursor()
    c.make_whole(center="origin")

    c = mols.cursor()
    c.make_whole(center=0)

    c = mols[0].atoms().cursor()
    c.make_whole(center=(1, 2, 3))

    mols = mols.clone()
    mols.make_whole(center=("1A", "2A", "3A"))
