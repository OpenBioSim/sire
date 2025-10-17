import sire as sr

import pytest


def test_amber_cmap(tmpdir, amber_cmap):
    mols = amber_cmap

    cmap = mols[0].property("cmap")

    params = cmap.parameters()

    # there are 18 cmap parameters in this file
    assert len(params) == 18

    unique_params = []

    for param in cmap.parameters():
        if param.parameter() not in unique_params:
            unique_params.append(param.parameter())

    # there are 9 different terms from this file
    assert len(unique_params) == 9

    # write the file and re-read it
    dir = tmpdir.mkdir("test_amber_cmap")

    f = sr.save(mols, dir.join("output"), format="prm7")

    mols2 = sr.load(f)

    cmap2 = mols2[0].property("cmap")

    params2 = cmap2.parameters()

    assert len(params2) == 18

    unique_params2 = []

    for param in cmap2.parameters():
        if param.parameter() not in unique_params2:
            unique_params2.append(param.parameter())

    assert len(unique_params2) == 9

    # check that the parameters are the same
    for param in params:
        assert param in params2

    # write again - this should produce an identical file,
    # as we made sure to sort indexes and parameters in
    # a consistent way
    f2 = sr.save(mols2, dir.join("output2"), format="prm7")

    with open(f[0], "r") as f1:
        with open(f2[0], "r") as f2:
            for i, (line1, line2) in enumerate(zip(f1, f2)):
                # skip the first line since it includes a timestamp
                if i > 0:
                    assert line1 == line2


def test_amber_cmap_grotop(tmpdir, amber_cmap):
    """Testing reading and writing from amber to gromacs and back to amber."""
    mols = amber_cmap.clone()

    dir = tmpdir.mkdir("test_amber_cmap_grotop")

    # Save to a temporary file in GroTop format.
    f = sr.save(mols, dir.join("output"), format="GroTop")

    # Load the saved file.
    mols2 = sr.load(f, show_warnings=False)

    # Save back to prmtop format.
    f2 = sr.save(mols2, dir.join("output"), format="prmtop")

    # Load the saved file again.
    mols3 = sr.load(f2, show_warnings=False)

    # check that all of the cmap parameters are the same
    cmaps = mols[0].property("cmap").parameters()
    cmaps2 = mols2[0].property("cmap").parameters()
    cmaps3 = mols3[0].property("cmap").parameters()

    for cmap in cmaps:
        # find the same cmap in the other two molecules
        found = False

        for cmap2 in cmaps2:
            if (
                cmap.atom0() == cmap2.atom0()
                and cmap.atom1() == cmap2.atom1()
                and cmap.atom2() == cmap2.atom2()
                and cmap.atom3() == cmap2.atom3()
                and cmap.atom4() == cmap2.atom4()
            ):
                # make sure that the parameter is the same
                for val, val2 in zip(
                    cmap.parameter().values(), cmap2.parameter().values()
                ):
                    assert val == pytest.approx(val2, 1e-3)
                found = True

        assert found

        found = False

        for cmap3 in cmaps3:
            if (
                cmap.atom0() == cmap3.atom0()
                and cmap.atom1() == cmap3.atom1()
                and cmap.atom2() == cmap3.atom2()
                and cmap.atom3() == cmap3.atom3()
                and cmap.atom4() == cmap3.atom4()
            ):
                # make sure that the parameter is the same
                for val, val3 in zip(
                    cmap.parameter().values(), cmap3.parameter().values()
                ):
                    assert val == pytest.approx(val3, 1e-3)
                found = True

        assert found
