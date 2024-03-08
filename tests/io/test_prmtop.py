import sire as sr

import pytest


def test_single_molecule_prmtop(tmpdir, ose_mols):
    mols = ose_mols

    # this makes sure it is now "OSE "
    assert mols[0].residues()[0].name() == sr.resid(name="OSE")

    dir = tmpdir.mkdir("test_single_molecule_prmtop")

    f = sr.save(mols, dir.join("test"), format=["prm"])

    mols2 = sr.load(f)

    assert mols[0].residues()[0].name() == sr.resid(name="OSE")


def test_reverse_dihedral(tmpdir, openmm_interchange_mols):
    mols = openmm_interchange_mols

    dir = tmpdir.mkdir("test_reverse_dihedral")

    f = sr.save(mols, dir.join("output"), format=["PRM7", "RST7"])

    # have to manually check that have written the two dihedrals between
    # atoms 2 and 5 (6 and 15) with only one having a 1-4 term, and the
    # other term ignored
    lines = open(f[0]).readlines()

    start = None

    for i, line in enumerate(lines):
        if line.find("%FLAG DIHEDRALS_WITHOUT_HYDROGEN") != -1:
            start = i + 2

    assert start is not None

    included = 0
    ignored = 0

    while True:
        line = lines[start]

        if line.startswith("%"):
            break

        start += 1

        atoms = line.split()
        [a, b, c, d] = [int(x) for x in atoms[0:4]]

        if (a == 6 and d == 15) or (a == 15 and d == 6):
            if c < 0:
                ignored += 1
            elif c > 0:
                included += 1

        try:
            [a, b, c, d] = [int(x) for x in atoms[5:9]]

            if (a == 6 and d == 15) or (a == 15 and d == 6):
                if c < 0:
                    ignored += 1
                elif c > 0:
                    included += 1
        except Exception as e:
            # this was the last line
            pass

    assert included == 1
    assert ignored == 1

    mols2 = sr.load(f)

    assert mols.energy().value() == pytest.approx(mols2.energy().value(), 1e-3)

    f = sr.save(mols, dir.join("output2"), format=["GroTop", "Gro87"])

    mols3 = sr.load(f, show_warnings=False)

    assert mols.energy().value() == pytest.approx(mols3.energy().value(), 1e-3)

    # manually check that the 2-5 pair appears only once
    lines = open(f[0]).readlines()

    start = None

    for i, line in enumerate(lines):
        if line.find("[ pairs ]") != -1:
            start = i + 2
            break

    assert start is not None

    included = 0

    while True:
        line = lines[start]
        start += 1

        atoms = line.split()
        if len(atoms) == 0:
            break

        [a, b] = [int(x) for x in atoms[0:2]]

        if (a == 2 and b == 5) or (a == 5 and b == 2):
            included += 1

    assert included == 1


def test_reorder_prmtop(reordered_protein):
    mols = reordered_protein

    # check that the atoms are in the correct order
    expected_number = 0

    for atom in mols.atoms():
        assert atom.number().value() == expected_number + 1
        expected_number += 1

    # check that the coordinates of the zinc atoms are correct
    zincs = mols.atoms("element Zn")

    assert len(zincs) == 4

    expected = [
        sr.v(7.646, 11.878, 28.867),
        sr.v(-0.332, 20.4, 19.036),
        sr.v(22.607, 10.037, 28.866),
        sr.v(32.863, 16.136, 38.415),
    ]

    for zinc, coord in zip(zincs, expected):
        z = zinc.coordinates()

        for i in range(3):
            assert z[i].value() == pytest.approx(coord[i].value(), 1e-3)

    # check that the zinc atoms are bonded to the proteins
    bonds = mols.bonds("element Zn")

    assert len(bonds) == 16

    zn = sr.mol.Element("Zn")

    for bond in bonds:
        assert bond[0].element() == zn or bond[1].element() == zn
        assert bond[0].molecule() == bond[1].molecule()

    assert mols[0]["element Zn"].num_atoms() == 2
    assert mols[1]["element Zn"].num_atoms() == 2
