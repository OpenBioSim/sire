import pytest
import sire as sr

from sire.legacy.Mol import AtomIdx


@pytest.mark.skipif(
    "rdkit" not in sr.convert.supported_formats(),
    reason="rdkit support is not available",
)
@pytest.mark.parametrize("size", [3, 6, 12, 13, 14, 30])
def test_in_ring_no_size_limit(size):
    mol = sr.smiles("C1" + "C" * (size - 2) + "C1")
    conn = mol.property("connectivity")

    for i in range(size):
        assert conn.in_ring(AtomIdx(i))
        assert conn.in_ring(AtomIdx(i), AtomIdx((i + 1) % size))

    # Hydrogens follow the heavy atoms.
    hydrogen = AtomIdx(mol.num_atoms() - 1)
    assert not conn.in_ring(hydrogen)
    assert not conn.in_ring(hydrogen, next(iter(conn.connections_to(hydrogen))))


@pytest.mark.skipif(
    "rdkit" not in sr.convert.supported_formats(),
    reason="rdkit support is not available",
)
def test_in_ring_separate_rings():
    # Biphenyl: ring A is atoms 0-5, ring B is atoms 6-11, joined by 5-6.
    conn = sr.smiles("c1ccccc1-c1ccccc1").property("connectivity")

    assert conn.in_ring(AtomIdx(5))
    assert conn.in_ring(AtomIdx(6))
    assert conn.in_ring(AtomIdx(1), AtomIdx(3))
    assert not conn.in_ring(AtomIdx(5), AtomIdx(6))
    assert not conn.in_ring(AtomIdx(1), AtomIdx(8))

    # Spiro[4.4]nonane: atom 3 is shared by rings 0-4 and 3,5-8.
    conn = sr.smiles("C1CCC2(C1)CCCC2").property("connectivity")

    assert conn.in_ring(AtomIdx(3), AtomIdx(0))
    assert conn.in_ring(AtomIdx(3), AtomIdx(6))
    assert not conn.in_ring(AtomIdx(0), AtomIdx(6))


@pytest.mark.skipif(
    "rdkit" not in sr.convert.supported_formats(),
    reason="rdkit support is not available",
)
def test_in_ring_after_edit():
    conn = sr.smiles("C1CCCCC1").property("connectivity")

    # Populate the cached ring systems before editing a copy.
    assert conn.in_ring(AtomIdx(0))

    broken = conn.edit().disconnect(AtomIdx(0), AtomIdx(1)).commit()

    assert conn.in_ring(AtomIdx(0))
    assert not broken.in_ring(AtomIdx(0))
    assert not broken.in_ring(AtomIdx(2), AtomIdx(3))

    closed = broken.edit().connect(AtomIdx(0), AtomIdx(1)).commit()

    assert closed.in_ring(AtomIdx(0))
    assert closed.in_ring(AtomIdx(2), AtomIdx(3))


@pytest.mark.skipif(
    "rdkit" not in sr.convert.supported_formats(),
    reason="rdkit support is not available",
)
@pytest.mark.parametrize("size", [3, 4, 6, 12, 13, 14])
def test_in_ring_max_ring_size(size):
    conn = sr.smiles("C1" + "C" * (size - 2) + "C1").property("connectivity")

    atoms = [AtomIdx(i) for i in range(min(size, 4))]

    assert conn.in_ring(*atoms[:3], max_ring_size=size)
    assert not conn.in_ring(*atoms[:3], max_ring_size=size - 1)
    assert conn.in_ring(*atoms[:3]) == (size <= 12)

    if size >= 4:
        assert conn.in_ring(*atoms, max_ring_size=size)
        assert not conn.in_ring(*atoms, max_ring_size=size - 1)
        assert conn.in_ring(*atoms) == (size <= 12)
