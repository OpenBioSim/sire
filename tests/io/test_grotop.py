import sire as sr

import pytest


def test_grospace(tmpdir, kigaki_mols):
    mols = kigaki_mols.clone()

    energy = mols.energy()

    mol = mols[0]

    mol = mol.edit().rename(" \t\t  NAME \t\t WITH     SPACE    ").commit()

    mols.update(mol)

    d = tmpdir.mkdir("test_grospace")

    f = sr.save(mols, d.join("test"), format=["Gro87", "GroTop"])

    lines = open(f[1], "r").readlines()

    found_correct_line = False

    for line in lines:
        # look for the name with underscores replacing spaces
        if line.find("NAME_WITH_SPACE") != -1:
            found_correct_line = True
            break

    assert found_correct_line

    mols = sr.load(f, show_warnings=False)

    assert energy.value() == pytest.approx(mols.energy().value())
