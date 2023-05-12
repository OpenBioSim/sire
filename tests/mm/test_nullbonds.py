import sire as sr
import pytest


def test_nullbonds(tmpdir, ala_mols):
    mols = ala_mols
    mol = mols[0]

    bonds = mol.property("bond")
    bond = mol.bonds()[0]
    bonds.set(bond.id(), 0)

    mol = mol.edit().set_property("bond", bonds).commit()

    params = sr.legacy.MM.AmberParams(mol)

    # check default is to keep the null bonds
    assert params.bond_functions().num_functions() == bonds.num_functions()

    p = params.bond_functions().potential(bond.id())

    assert p.is_zero()

    func = params.get_parameter(bond.id())

    # this will be a null bond with r0 equal to the bond length
    assert func.k() == 0
    assert func.r0() == bond.length().to(sr.units.angstrom)

    params = sr.legacy.MM.AmberParams(
        mol, map=sr.base.create_map({"keep_null_bonds": True})
    )

    assert params.bond_functions().num_functions() == bonds.num_functions()
    assert params.get_parameter(bond.id()).k() == 0

    # now test removing the null bonds
    params = sr.legacy.MM.AmberParams(
        mol, map=sr.base.create_map({"keep_null_bonds": False})
    )

    assert params.bond_functions().num_functions() == bonds.num_functions() - 1

    # this will be a default-constructed empty bond
    func = params.get_parameter(bond.id())

    assert func.k() == 0
    assert func.r0() == 0

    # now test writing to a top file
    d = tmpdir.mkdir("test_nullbonds")

    f = sr.save(mol, d.join("default"), format="PRM7")

    mol2 = sr.load(f[0])[0]

    params = mol2.property("parameters")
    assert params.bond_functions().num_functions() == bonds.num_functions()
    func = params.get_parameter(bond.id())
    assert func.k() == 0
    assert func.r0() == pytest.approx(bond.length().to(sr.units.angstrom))

    f = sr.save(mol, d.join("skip"), format="PRM7", map={"keep_null_bonds": False})

    mol2 = sr.load(f[0])[0]

    params = mol2.property("parameters")
    assert params.bond_functions().num_functions() == bonds.num_functions() - 1
    func = params.get_parameter(bond.id())
    assert func.k() == 0
    assert func.r0() == 0
