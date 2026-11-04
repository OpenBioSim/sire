import pytest
import sire as sr


@pytest.mark.parametrize(
    "build",
    [
        lambda m: sr.restraints.distance(m, atoms0=m[0][0], atoms1=m[0][1]),
        lambda m: sr.restraints.inverse_distance(m, atoms0=m[0][0], atoms1=m[0][1]),
        lambda m: sr.restraints.positional(m, atoms=m[0][0]),
        lambda m: sr.restraints.angle(m, atoms=m[0][0:3]),
        lambda m: sr.restraints.dihedral(m, atoms=m[0][0:4]),
    ],
)
def test_restraints_add(ala_mols, build):
    restraints = build(ala_mols)

    assert len(restraints) == 1
    assert len(restraints + restraints) == 2
    assert len(restraints + restraints[0]) == 2
