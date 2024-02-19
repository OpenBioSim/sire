import sire as sr

import pytest


def test_match(kigaki_mols):
    mols = kigaki_mols

    mol = mols.molecule("protein")

    ala = mol.residues("ALA")[1]
    phe = mol.residues("LYS")[1]

    m = sr.morph.match(ala, phe, match_light_atoms=True)

    assert m is not None

    assert len(m) == 9

    assert len(m.mapped_atoms0()) == 9
    assert len(m.mapped_atoms1()) == 9

    assert len(m.unmapped_atoms0().atoms()) == 1
    assert len(m.unmapped_atoms1().atoms()) == 13

    assert len(m.atoms0()) == len(ala.atoms())
    assert len(m.atoms1()) == len(phe.atoms())

    assert m.atoms0() == ala.atoms()
    assert m.atoms1() == phe.atoms()

    assert m.is_single_molecule()

    assert not m.is_empty()

    num_unmapped = 0

    for atom in ala.atoms():
        assert m.contains(atom)

        if atom.element().num_protons() != 1:
            assert m.is_mapped(atom)

        if not m.is_mapped(atom):
            num_unmapped += 1

    assert num_unmapped == 1

    m = m.swap()

    assert len(m) == 9

    assert len(m.mapped_atoms0()) == 9
    assert len(m.mapped_atoms1()) == 9

    assert len(m.unmapped_atoms0().atoms()) == 13
    assert len(m.unmapped_atoms1().atoms()) == 1

    assert len(m.atoms0()) == len(phe.atoms())
    assert len(m.atoms1()) == len(ala.atoms())

    assert m.atoms0() == phe.atoms()
    assert m.atoms1() == ala.atoms()

    assert m.is_single_molecule()

    assert not m.is_empty()

    num_unmapped = 0

    for atom in phe.atoms():
        assert m.contains(atom)

        if not m.is_mapped(atom):
            num_unmapped += 1

    assert m.atoms0().to_single_molecule() == phe.atoms().to_single_molecule()
    assert m.atoms1().to_single_molecule() == ala.atoms().to_single_molecule()
