import copy

import pytest
import sire as sr


@pytest.mark.parametrize(
    "obj",
    [
        sr.u("1 ps"),
        sr.legacy.Mol.AtomIdx(3),
        sr.cas.LambdaSchedule.standard_morph(),
    ],
)
def test_deepcopy(obj):
    other = copy.deepcopy(obj)
    assert other == obj
    assert other is not obj


def test_deepcopy_openmm():
    from sire.legacy.Convert import _SireOpenMM

    lever = _SireOpenMM.LambdaLever()
    assert copy.deepcopy(lever) is not lever


def test_deepcopy_container():
    d = {"time": sr.u("2 ps"), "atoms": [sr.legacy.Mol.AtomIdx(1)]}
    assert copy.deepcopy(d) == d
