import pytest
import sire as sr


def test_morse_potential_restraints_setup(cyclopentane_cyclohexane):
    """Tests that morse_potential restraints can be set up correctly with custom parameters."""
    mols = cyclopentane_cyclohexane.clone()
    restraints, mols = sr.restraints.morse_potential(
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


def test_morse_potential_restraint_annihiliation_auto_param(cyclopentane_cyclohexane):
    """Tests that morse_potential restraints can be set up correctly with automatic parametrisation
    when a bond is to be annihilated."""
    mols = cyclopentane_cyclohexane.clone()
    restraints, mols = sr.restraints.morse_potential(
        mols,
        de="25 kcal mol-1",
        auto_parametrise=True,
    )
    assert restraints.num_restraints() == 1
    assert restraints[0].atom0() == 0
    assert restraints[0].atom1() == 4
    assert restraints[0].k().value() == pytest.approx(457.78, rel=0.1)
    assert restraints[0].de().value() == 25.0


def test_morse_potential_restraint_creation_auto_param(propane_cyclopropane):
    """Tests that morse_potential restraints can be set up correctly with automatic parametrisation
    when a bond is to be created."""
    mols = propane_cyclopropane.clone()
    restraints, mols = sr.restraints.morse_potential(
        mols,
        de="25 kcal mol-1",
        auto_parametrise=True,
    )
    assert restraints.num_restraints() == 1
    assert restraints[0].atom0() == 0
    assert restraints[0].atom1() == 2


def test_morse_potential_restraint_auto_param_override(cyclopentane_cyclohexane):
    """Tests that morse_potential restraints can be set up correctly with automatic parametrisation and some parameters can be overwritten."""
    mols = cyclopentane_cyclohexane.clone()
    restraints, mols = sr.restraints.morse_potential(
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
    restraints, mols = sr.restraints.morse_potential(
        mols,
        de="25 kcal mol-1",
        auto_parametrise=True,
        direct_morse_replacement=False,
    )
    restraint1, mols = sr.restraints.morse_potential(
        mols,
        atoms0=mols["molecule property is_perturbable and atomidx 0"],
        atoms1=mols["molecule property is_perturbable and atomidx 1"],
        k="100 kcal mol-1 A-2",
        r0="1.5 A",
        de="50 kcal mol-1",
        direct_morse_replacement=False,
    )
    restraints.add(restraint1)
    assert restraints.num_restraints() == 2

def test_morse_potential_direct_morse_replacement(cyclopentane_cyclohexane):
    """Tests that morse_potential restraints by default will remove the annihilated harmonic bond."""
    mols = cyclopentane_cyclohexane.clone()
    bonds0_org_mol = mols[0].property("bond0")
    bonds1_org_mol = mols[0].property("bond1")
    num_bonds0_org_mol = len(bonds0_org_mol.potentials())
    num_bonds1_org_mol = len(bonds1_org_mol.potentials())

    restraints, mols = sr.restraints.morse_potential(
        mols,
        de="25 kcal mol-1",
        auto_parametrise=True,
    )
   
    bonds0_mod_mol = mols[0].property("bond0")
    bonds1_mod_mol = mols[0].property("bond1")
    num_bonds0_mod_mol = len(bonds0_mod_mol.potentials())
    num_bonds1_mod_mol = len(bonds1_mod_mol.potentials())

    # New bonds0 should be smaller by 1 than the original bonds0 if the function removed the bond
    assert num_bonds0_mod_mol == num_bonds0_org_mol - 1

    # New bonds1 should be same as the original bonds1, as the bonds1 will already be smaller
    # by 1 than the original bonds0 which is introduced during the merge.
    assert num_bonds1_mod_mol == num_bonds1_org_mol
    assert num_bonds0_org_mol == num_bonds1_org_mol + 1

def test_morse_potential_direct_morse_replacement_retain_harmonic(cyclopentane_cyclohexane):
    """Tests that morse_potential restraints can retain the annihilated harmonic bond, if specified."""
    mols = cyclopentane_cyclohexane.clone()
    bonds0_org_mol = mols[0].property("bond0")
    bonds1_org_mol = mols[0].property("bond1")
    num_bonds0_org_mol = len(bonds0_org_mol.potentials())
    num_bonds1_org_mol = len(bonds1_org_mol.potentials())

    restraints, mols = sr.restraints.morse_potential(
        mols,
        de="25 kcal mol-1",
        auto_parametrise=True,
        retain_harmonic_bond=True,
    )
    bonds0_mod_mol = mols[0].property("bond0")
    bonds1_mod_mol = mols[0].property("bond1")
    num_bonds0_mod_mol = len(bonds0_mod_mol.potentials())
    num_bonds1_mod_mol = len(bonds1_mod_mol.potentials())

    # If the harmonic bond is retained, the number of bonds in bonds0 and bonds1
    # should be the same as the original molecule
    assert num_bonds0_mod_mol == num_bonds0_org_mol
    assert num_bonds1_mod_mol == num_bonds1_org_mol
