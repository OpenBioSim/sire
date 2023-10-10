import sire as sr

import pytest
import sys


@pytest.mark.skipif(sys.platform == "win32", reason="Does not run on windows")
def test_mapping(ala_mols, ala_traj):
    mols = ala_mols.clone()
    mols2 = ala_traj.clone()

    mol = mols[0]
    mol2 = mols2[0]

    mapping = sr.match_atoms(mol, mol2, match_light_atoms=True)

    for i, atom in enumerate(mapping[mol.atoms()]):
        assert atom.element() == mol[i].element()

        if atom.element().num_protons() != 1:
            # the mapping can match different hydrogens
            assert atom.name() == mol[i].name()

    # should be the same order
    assert mapping[mol["not element H"]] == mol2["not element H"]

    assert mol["not element H"] != mol2["not element H"]

    # now map just 5 atoms
    assert mapping[mol[0:5]]["not element H"] == mol2[0:5]["not element H"]

    # now map the molecule in the wrong order
    assert (
        mapping[mol[-1::-1]]["not element H"] == mol2[-1::-1]["not element H"]
    )

    # now try to find the matching atoms in the container
    assert mapping.find(mol, mols2)["not element H"] == mol2["not element H"]

    assert (
        mapping.find(mol[0:5], mols2)["not element H"]
        == mol2[0:5]["not element H"]
    )

    assert (
        mapping.find(mol[-1::-1], mols2)["not element H"]
        == mol2[-1::-1]["not element H"]
    )

    # checked that swapping the match works
    assert mapping.swap()[mol2["not element H"]] == mol["not element H"]
