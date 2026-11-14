import copy

import pytest


@pytest.mark.parametrize(
    "create",
    [
        lambda sr: sr.u("1 ps"),
        lambda sr: sr.legacy.Mol.AtomIdx(3),
        lambda sr: sr.cas.LambdaSchedule.standard_morph(),
    ],
)
def test_deepcopy(create):
    import sire as sr

    obj = create(sr)
    other = copy.deepcopy(obj)
    assert other == obj
    assert other is not obj


def test_deepcopy_no_memo():
    import sire as sr

    obj = sr.u("1 ps")
    other = obj.__deepcopy__()
    assert other == obj
    assert other is not obj


def test_deepcopy_openmm():
    from sire.legacy.Convert import _SireOpenMM

    lever = _SireOpenMM.LambdaLever()
    assert copy.deepcopy(lever) is not lever


def test_deepcopy_container():
    import sire as sr

    d = {"time": sr.u("2 ps"), "atoms": [sr.legacy.Mol.AtomIdx(1)]}
    assert copy.deepcopy(d) == d
