import pytest

import sire as sr


def test_search_closest(ala_mols):
    mols = ala_mols

    C = mols[0][1]
    H = mols[0][2]

    assert mols["closest atom to (molidx 0 and atomidx 1)"] == H
    assert mols["closest atom to (molidx 0 and atomidx 2)"] == C

    five_waters = mols["closest 5 waters to resname ALA"]

    assert len(five_waters) == 5

    space = mols.space()

    mindists = [
        sr.minimum_distance(x, mols["resname ALA"], space) for x in five_waters
    ]

    # the five distances should be in order
    for i in range(0, 4):
        assert mindists[i].value() < mindists[i + 1].value()

    assert mindists[-1] < sr.u("3 angstrom")

    # now get the next 100 distances
    hundred_waters = mols["closest 100 waters to molidx 0"]

    assert len(hundred_waters) == 100

    mindists = [sr.minimum_distance(x, mols[0], space) for x in hundred_waters]

    for i in range(0, 99):
        assert mindists[i].value() < mindists[i + 1].value()

    assert mindists[-1] < sr.u("6.5 A")

    center = mols[0].coordinates()

    assert mols[f"closest molecule to {center:f}"] == mols[0]


def test_search_furthest(ala_mols):
    mols = ala_mols

    C = mols[0][1]

    F = mols["furthest atom from (molidx 0 and atomidx 1)"]

    assert sr.minimum_distance(C, F) > sr.u("20 A")

    five_waters = mols["furthest 5 waters from resname ALA"]

    assert len(five_waters) == 5

    mindists = [
        sr.minimum_distance(x, mols["resname ALA"]) for x in five_waters
    ]

    for i in range(0, 4):
        assert mindists[i].value() > mindists[i + 1].value()

    assert mindists[0] > sr.u("20 A")

    hundred_waters = mols["furthest 100 waters from molidx 0"]

    assert len(hundred_waters) == 100

    space = mols.space()

    mindists = [sr.minimum_distance(x, mols[0], space) for x in hundred_waters]

    for i in range(0, 99):
        assert mindists[i].value() > mindists[i + 1].value()

    assert mindists[0] > sr.u("20 A")

    center = mols[0][0].coordinates()

    assert (
        mols[f"furthest molecule from {center:f}"]
        == mols["furthest molecule from (molidx 0 and atomidx 0)"]
    )
