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
            assert f1.read() == f2.read()
