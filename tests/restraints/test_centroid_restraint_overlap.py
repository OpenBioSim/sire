import pytest
import sire as sr


@pytest.mark.parametrize(
    "build",
    [
        lambda a0, a1: sr.mm.BondRestraint(
            a0, a1, sr.u("1 kcal mol-1 A-2"), sr.u("1 A")
        ),
        lambda a0, a1: sr.mm.InverseBondRestraint(
            a0, a1, sr.u("1 kcal mol-1 A-2"), sr.u("1 A")
        ),
        lambda a0, a1: sr.mm.MorsePotentialRestraint(
            a0, a1, sr.u("1 kcal mol-1 A-2"), sr.u("1 A"), sr.u("1 kcal mol-1")
        ),
    ],
)
def test_centroid_restraint_overlap_names_atom(build):
    restraint = build([1, 2], [3, 4])
    assert restraint.is_centroid_restraint()

    with pytest.raises(Exception, match="Atom 7 appears in both"):
        build([5, 7], [3, 7])
