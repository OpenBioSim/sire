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


def test_lj_1264_exceptions(apo_1264):
    mols = apo_1264

    # mols[3] is an ion with LJ 12-6-4 exceptions - check this is the case
    ex = mols[0].property("LJ").get_exceptions(mols[1].property("LJ"))
    assert len(ex) == 0

    ex = mols[0].property("LJ").get_exceptions(mols[2].property("LJ"))
    assert len(ex) == 0

    ex = mols[0].property("LJ").get_exceptions(mols[3].property("LJ"))
    assert len(ex) == 4647

    # for some reason, there aren't any 12-6-4 exceptions involving
    # amber type HO
    assert len(ex) == mols[0]["atom property ambertype != HO"].num_atoms()

    # there should be one exception with a water molecule - check this
    # for a couple of water molecules
    water_ex = mols["water[0]"].property("LJ").get_exceptions(mols[3].property("LJ"))
    assert len(water_ex) == 1

    water_ex2 = mols["water[10]"].property("LJ").get_exceptions(mols[3].property("LJ"))
    assert len(water_ex2) == 1

    assert water_ex == water_ex2

    assert water_ex[0][0] == 0
    assert water_ex[0][1] == 0
    assert water_ex[0][2] == sr.mm.LJ1264Parameter(
        "127672 kcal mol-1 Å^12", "1497.12 kcal mol-1 Å^6", "1 kcal mol-1 Å^4"
    )

    # there should be one exception with the ion - its self 12-6-4 term
    self_ex = mols[3].property("LJ").get_exceptions(mols[3].property("LJ"))
    assert len(self_ex) == 1


@pytest.mark.slow
def test_lj_1264_load_save(tmpdir, apo_1264):
    mols = apo_1264

    dir = tmpdir.mkdir("test_lj_1264_load_save")

    f = sr.save(mols, dir.join("output"), format=["prm7", "rst7"])

    mols2 = sr.load(f)

    assert mols[0:10].energy().value() == pytest.approx(
        mols2[0:10].energy().value(), 1e-3
    )

    # mols[3] is an ion with LJ 12-6-4 exceptions - check this is the case
    ex = mols2[0].property("LJ").get_exceptions(mols2[1].property("LJ"))
    assert len(ex) == 0

    ex = mols2[0].property("LJ").get_exceptions(mols2[2].property("LJ"))
    assert len(ex) == 0

    ex = mols2[0].property("LJ").get_exceptions(mols2[3].property("LJ"))
    assert len(ex) == 4647

    # for some reason, there aren't any 12-6-4 exceptions involving
    # amber type HO
    assert len(ex) == mols2[0]["atom property ambertype != HO"].num_atoms()

    # there should be one exception with a water molecule - check this
    # for a couple of water molecules
    water_ex = mols2["water[0]"].property("LJ").get_exceptions(mols2[3].property("LJ"))
    assert len(water_ex) == 1

    water_ex2 = (
        mols2["water[10]"].property("LJ").get_exceptions(mols2[3].property("LJ"))
    )
    assert len(water_ex2) == 1

    assert water_ex == water_ex2

    # there should be one exception with the ion - its self 12-6-4 term
    self_ex = mols2[3].property("LJ").get_exceptions(mols2[3].property("LJ"))
    assert len(self_ex) == 1
