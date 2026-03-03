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


def test_amber_multichain_cmap(tmpdir, multichain_cmap):
    """Regression test for multi-chain CMAP write bug.

    When a topology has more than one molecule with CMAP terms, the write path
    previously applied a spurious '/= 3' to the per-molecule atom offset,
    producing wrong global atom indices for molecules 2, 3, … .  The re-parse
    of the generated text then raised:
      "there is a cmap between more than one different molecule"
    """
    mols = multichain_cmap

    # Collect per-molecule CMAP term counts indexed by position in the
    # molecule list. There must be at least two molecules with CMAP
    # to trigger the multi-chain bug.
    cmap_counts = {}
    for i, mol in enumerate(mols.molecules()):
        if mol.has_property("cmap"):
            cmap_counts[i] = len(mol.property("cmap").parameters())

    assert (
        len(cmap_counts) >= 2
    ), "Expected at least two molecules with CMAP terms in this topology"

    dir = tmpdir.mkdir("test_amber_multichain_cmap")

    # Write the topology back to file. (This is where the bug used to trigger.)
    f = sr.save(mols, dir.join("output"), format="prm7")

    # Reead the file back in.
    mols2 = sr.load(f)

    # Each molecule must still have the same number of CMAP terms after the
    # roundtrip.
    for i, count in cmap_counts.items():
        mol2 = mols2[i]
        assert mol2.has_property(
            "cmap"
        ), f"Molecule at index {i} lost its cmap property after roundtrip"
        count2 = len(mol2.property("cmap").parameters())
        assert (
            count2 == count
        ), f"Molecule at index {i}: CMAP count changed from {count} to {count2}"

    # Verify a second write also succeeds without error.
    sr.save(mols2, dir.join("output2"), format="prm7")


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
