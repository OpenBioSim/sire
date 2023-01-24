import sire as sr
import pytest


def test_nullangles(tmpdir, ala_mols):
    mols = ala_mols
    mol = mols[0]

    angles = mol.property("angle")
    angle = mol.angles()[0]
    angles.set(angle.id(), 0)

    mol = mol.edit().set_property("angle", angles).commit()

    params = sr.legacy.MM.AmberParams(mol)

    # check default is to keep the null angles
    assert params.angle_functions().num_functions() == angles.num_functions()

    p = params.angle_functions().potential(angle.id())

    assert p.is_zero()

    func = params.get_parameter(angle.id())

    # this will be a null angle with r0 equal to the angle size
    assert func.k() == 0
    assert func.theta0() == angle.size().to(sr.units.radians)

    params = sr.legacy.MM.AmberParams(
        mol, map=sr.base.create_map({"keep_null_angles": True})
    )

    assert params.angle_functions().num_functions() == angles.num_functions()
    assert params.get_parameter(angle.id()).k() == 0

    # now test removing the null angles
    params = sr.legacy.MM.AmberParams(
        mol, map=sr.base.create_map({"keep_null_angles": False})
    )

    assert (
        params.angle_functions().num_functions() == angles.num_functions() - 1
    )

    # this will be a default-constructed empty angle
    func = params.get_parameter(angle.id())

    assert func.k() == 0
    assert func.theta0() == 0

    # now test writing to a top file
    d = tmpdir.mkdir("test_nullangles")

    f = sr.save(mol, d.join("default"), format="PRM7")

    mol2 = sr.load(f[0])[0]

    params = mol2.property("parameters")
    assert params.angle_functions().num_functions() == angles.num_functions()
    func = params.get_parameter(angle.id())
    assert func.k() == 0
    assert func.theta0() == pytest.approx(angle.size().to(sr.units.radians))

    f = sr.save(
        mol, d.join("skip"), format="PRM7", map={"keep_null_angles": False}
    )

    mol2 = sr.load(f[0])[0]

    params = mol2.property("parameters")
    assert (
        params.angle_functions().num_functions() == angles.num_functions() - 1
    )
    func = params.get_parameter(angle.id())
    assert func.k() == 0
    assert func.theta0() == 0
