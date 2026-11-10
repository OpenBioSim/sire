import pytest
import sire as sr

from sire.legacy.Mol import AtomIdx, AtomIdxMatcher


@pytest.mark.skipif(
    "rdkit" not in sr.convert.supported_formats(),
    reason="rdkit support is not available",
)
def test_make_compatible_removes_bonds_to_removed_atoms():
    mol = sr.smiles("C1CCCCC1")
    conn = mol.property("connectivity")

    last = AtomIdx(mol.num_atoms() - 1)
    neighbour = next(iter(conn.connections_to(last)))

    info = mol.edit().remove(last).commit().info()

    # AtomIdxMatcher never reports a change of order, even when atoms are removed.
    new = conn.make_compatible_with(info, AtomIdxMatcher())

    assert last not in new.connections_to(neighbour)
    assert new.in_ring(neighbour)


@pytest.mark.skipif(
    "rdkit" not in sr.convert.supported_formats(),
    reason="rdkit support is not available",
)
def test_rename_keeps_internal_properties():
    mol = sr.smiles("CCCC")

    cursor = mol.cursor()
    cursor.bonds()[0]["note"] = "bond"
    cursor.angles()[0]["note"] = "angle"
    cursor.dihedrals()[0]["note"] = "dihedral"
    mol = cursor.commit()

    mol = mol.edit().atom(0).rename("X").molecule().commit()

    conn = mol.property("connectivity")

    assert conn.property(mol.bonds()[0].id(), "note").value() == "bond"
    assert conn.property(mol.angles()[0].id(), "note").value() == "angle"
    assert conn.property(mol.dihedrals()[0].id(), "note").value() == "dihedral"
