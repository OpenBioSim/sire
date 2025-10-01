import pytest
import sire as sr


def test_morse_potential_restraints_setup(cyclopentane_cyclohexane):
    """Tests that morse_potential restraints can be set up correctly with custom parameters."""
    mols = cyclopentane_cyclohexane.clone()
    restraints = sr.restraints.morse_potential(
        mols,
        atoms0=mols["molecule property is_perturbable and atomidx 0"],
        atoms1=mols["molecule property is_perturbable and atomidx 4"],
        k="100 kcal mol-1 A-2",
        r0="1.5 A",
        de="50 kcal mol-1",
        )
    assert restraints.num_restraints() == 1
    assert restraints[0].atom0() == 0
    assert restraints[0].atom1() == 4
    assert restraints[0].k().value() == 100.0
    assert restraints[0].r0().value() == 1.5
    assert restraints[0].de().value() == 50.0


def test_morse_potential_restraint_auto_param(cyclopentane_cyclohexane):
    """Tests that morse_potential restraints can be set up correctly with automatic parametrisation."""
    mols = cyclopentane_cyclohexane.clone()
    restraints = sr.restraints.morse_potential(
        mols,
        de="25 kcal mol-1",
        auto_parametrise=True,
        )
    assert restraints.num_restraints() == 1
    assert restraints[0].atom0() == 0
    assert restraints[0].atom1() == 4
    assert restraints[0].k().value() == pytest.approx(228.89, rel=0.1)
    assert restraints[0].de().value() == 25.0

def test_morse_potential_restraint_auto_param_override(cyclopentane_cyclohexane):
    """Tests that morse_potential restraints can be set up correctly with automatic parametrisation and some parameters can be overwritten."""
    mols = cyclopentane_cyclohexane.clone()
    restraints = sr.restraints.morse_potential(
        mols,
        de="25 kcal mol-1",
        k="100 kcal mol-1 A-2",
        auto_parametrise=True,
        )
    assert restraints.num_restraints() == 1
    assert restraints[0].atom0() == 0
    assert restraints[0].atom1() == 4
    assert restraints[0].k().value() == 100.0
    assert restraints[0].de().value() == 25.0


def test_multiple_morse_potential_restraints(cyclopentane_cyclohexane):
    """Tests that multiple morse_potential restraints can be set up correctly."""
    mols = cyclopentane_cyclohexane.clone()
    restraints = sr.restraints.morse_potential(
        mols,
        de="25 kcal mol-1",
        auto_parametrise=True,
        )
    restraint1 = sr.restraints.morse_potential(
    mols,
    atoms0=mols["molecule property is_perturbable and atomidx 0"],
    atoms1=mols["molecule property is_perturbable and atomidx 1"],
    k="100 kcal mol-1 A-2",
    r0="1.5 A",
    de="50 kcal mol-1",
    )
    restraints.add(restraint1)
    assert restraints.num_restraints() == 2