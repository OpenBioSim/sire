import pytest
import sire as sr


def test_make_whole(wrapped_mols):
    mols = wrapped_mols.clone()

    wrong_center = mols[0].evaluate().center_of_mass()

    mols.make_whole()

    right_center = mols[0].evaluate().center_of_mass()

    assert wrong_center != right_center

    correct_com = sr.maths.Vector(31.07343482, 31.89108144, 0.87194654)

    assert right_center.x().value() == pytest.approx(
        correct_com.x().value(), 1e-3
    )

    assert right_center.y().value() == pytest.approx(
        correct_com.y().value(), 1e-3
    )

    assert right_center.z().value() == pytest.approx(
        correct_com.z().value(), 1e-3
    )

    mols = wrapped_mols

    mol = mols[0]

    assert mol.evaluate().center_of_mass() == wrong_center

    mol = mol.move().make_whole().commit()

    assert mol.evaluate().center_of_mass() == right_center

    c = wrapped_mols.cursor()

    assert c.commit()[0].evaluate().center_of_mass() == wrong_center

    c.make_whole()

    mol = c.commit()[0].evaluate().center_of_mass() == right_center
