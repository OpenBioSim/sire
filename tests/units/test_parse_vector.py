import pytest
import sire as sr


@pytest.mark.parametrize(
    "args, expect",
    [
        ([0], sr.maths.Vector(0)),
        ([1, 2, 3], sr.maths.Vector(1, 2, 3)),
        (["1 A", "2 A", "3 A"], sr.maths.Vector(1, 2, 3)),
        ([1, 2, 3, "A"], sr.maths.Vector(1, 2, 3)),
        ([sr.u("1 A"), sr.u("2 A"), sr.u("3 A")], sr.maths.Vector(1, 2, 3)),
        (
            [3, 4, 5, "A ps-1"],
            sr.legacy.Mol.Velocity3D("3 A ps-1", "4 A ps-1", "5 A ps-1"),
        ),
        (
            [3, 4, 5, "newton"],
            sr.legacy.Mol.Force3D("3 newton", "4 newton", "5 newton"),
        ),
        (
            [1, 2, 3, "oC"],
            (sr.u("1 oC"), sr.u("2 oC"), sr.u("3 oC")),
        ),
    ],
)
def test_pase_vector(args, expect):
    assert sr.v(*args) == expect
