import pytest
import sire as sr


def test_default_rmsd_restraints_setup(ala_mols):
    """Tests that rmsd restraints can be set up correctly with default parameters."""
    mols = ala_mols.clone()
    restraints = sr.restraints.rmsd(mols=mols, atoms=[0, 1, 2, 3, 4])
    assert restraints.num_restraints() == 1
    assert restraints[0].atoms() == [0, 1, 2, 3, 4]
    assert restraints[0].k().value() == 150.0
    assert restraints[0].r0().value() == 0.0

def test_rmsd_restraint_custom_k(ala_mols):
    """Tests that rmsd restraints can be set up correctly with custom k."""
    mols = ala_mols.clone()
    restraints = sr.restraints.rmsd(
        mols=mols, atoms=[0, 1, 2, 3, 4], k="10 kcal mol-1 A-2"
    )
    assert restraints.num_restraints() == 1
    assert restraints[0].atoms() == [0, 1, 2, 3, 4]
    assert restraints[0].k().value() == 10.0
    assert restraints[0].r0().value() == 0.0

def test_rmsd_restraint_custom_k_and_r0(ala_mols):
    """Tests that rmsd restraints can be set up correctly with custom k and r0."""
    mols = ala_mols.clone()
    restraints = sr.restraints.rmsd(
        mols=mols, atoms="atomname CA, C, N", k="10 kcal mol-1 A-2", r0="1 A"
    )
    assert restraints.num_restraints() == 1
    assert restraints[0].k().value() == 10.0
    assert restraints[0].r0().value() == 1.0

def test_rmsd_restraint_custom_k_and_r0_and_ref(ala_mols):
    """Tests that rmsd restraints can be set up correctly with custom k, r0 and reference state."""
    mols = ala_mols.clone()
    ref = mols.clone().minimisation().run().commit()
    restraints = sr.restraints.rmsd(
        mols=mols, atoms="atomname CA, C, N", k="10 kcal mol-1 A-2", r0="1 A", ref=ref
    )
    assert restraints.num_restraints() == 1
    assert mols.count() == ref.count()
    assert restraints[0].k().value() == 10.0
    assert restraints[0].r0().value() == 1.0

def test_multiple_rmsd_restraints(ala_mols):
    """Tests that multiple RMSD restraints can be set up correctly."""
    mols = ala_mols.clone()
    restraints = sr.restraints.rmsd(mols=mols, atoms=[0, 1, 2])
    restraint1 = sr.restraints.rmsd(mols=mols, atoms=[3, 4, 5])
    restraints.add(restraint1)
    assert restraints.num_restraints() == 2