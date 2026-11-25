from sire.legacy.IO import patchIntrascale
from sire.legacy.MM import CLJNBPairs, CLJScaleFactor
from sire.legacy.Mol import AtomIdx


def _same(a, b):
    return a.coulomb() == b.coulomb() and a.lj() == b.lj()


def test_patch_intrascale(ala_mols):
    mol = ala_mols[0]
    n = mol.num_atoms()

    # Pairs scaled in only one of nb and base, so that both scans are needed.
    nb = mol.property("intrascale")
    pair14 = next(
        (AtomIdx(i), AtomIdx(j))
        for i in range(n)
        for j in range(i + 1, n)
        if _same(nb.get(AtomIdx(i), AtomIdx(j)), CLJScaleFactor(0.5 / 0.6, 0.5))
    )
    nb.set(*pair14, CLJScaleFactor(1, 1))

    base = CLJNBPairs(mol.property("connectivity"), CLJScaleFactor(0.5 / 0.6, 0.5))
    far = (AtomIdx(0), AtomIdx(n - 1))
    assert _same(base.get(*far), CLJScaleFactor(1, 1))
    base.set(*far, CLJScaleFactor(0.3, 0.7))

    # Leave the last few atoms unmapped.
    nmapped = n - 3

    far_mapped = (AtomIdx(0), AtomIdx(nmapped - 1))
    assert _same(nb.get(*far_mapped), CLJScaleFactor(1, 1))
    assert _same(base.get(*far_mapped), CLJScaleFactor(1, 1))
    nb.set(*far_mapped, CLJScaleFactor(0.4, 0.6))
    mapping = {AtomIdx(i): AtomIdx(i) for i in range(nmapped)}

    full = {AtomIdx(i): AtomIdx(i) for i in range(n)}
    patched, patched_full = patchIntrascale(nb, nb, base, base, mapping, full)

    for i in range(n):
        for j in range(i, n):
            ai, aj = AtomIdx(i), AtomIdx(j)

            assert _same(patched_full.get(ai, aj), nb.get(ai, aj))

            if i < nmapped and j < nmapped:
                assert _same(patched.get(ai, aj), nb.get(ai, aj))
            else:
                assert _same(patched.get(ai, aj), base.get(ai, aj))

    assert _same(patched.get(*pair14), CLJScaleFactor(1, 1))
    assert _same(patched.get(*far_mapped), CLJScaleFactor(0.4, 0.6))
