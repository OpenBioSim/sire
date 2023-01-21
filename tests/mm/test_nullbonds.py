import sire as sr
import pytest


def test_nullbonds(ala_mols):
    mols = ala_mols
    mol = mols[0]

    bonds = mol.property("bond")
    bonds.set(mol.bonds()[0].id(), 0)

    mol = mol.edit().set_property("bond", bonds).commit()

    params = sr.legacy.MM.AmberParams(mol)

    assert params.bond_functions().num_functions() == bonds.num_functions() - 1

    params = sr.legacy.MM.AmberParams(
        mol, map=sr.base.create_map({"include_null_bonds": True})
    )

    assert params.bond_functions().num_functions() == bonds.num_functions()

    p = params.bond_functions().potential(mol.bonds()[0].id())

    assert p.is_zero()
