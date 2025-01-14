import pytest
import sire as sr


def test_default_angle_restraints_setup(ala_mols):
    """Tests that angle restraints can be set up correctly with default parameters."""
    mols = ala_mols.clone()
    restraints = sr.restraints.angle(mols=mols, atoms=[0, 1, 2])
    assert restraints.num_restraints() == 1
    assert restraints[0].atoms() == [0, 1, 2]
    assert restraints[0].ktheta().value() == 100.0


def test_angle_restraint_custom_ktheta(ala_mols):
    """Tests that angle restraints can be set up correctly with custom ktheta."""
    mols = ala_mols.clone()
    restraints = sr.restraints.angle(
        mols=mols, atoms=[0, 1, 2], ktheta="10 kcal mol-1 rad-2"
    )
    assert restraints.num_restraints() == 1
    assert restraints[0].atoms() == [0, 1, 2]
    assert restraints[0].ktheta().value() == 10.0


def test_angle_restraint_custom_ktheta_and_theta0(ala_mols):
    """Tests that angle restraints can be set up correctly with custom ktheta and theta0."""
    mols = ala_mols.clone()
    restraints = sr.restraints.angle(
        mols=mols, atoms=[0, 1, 2], ktheta="10 kcal mol-1 rad-2", theta0="2 rad"
    )
    assert restraints.num_restraints() == 1
    assert restraints[0].atoms() == [0, 1, 2]
    assert restraints[0].ktheta().value() == 10.0
    assert restraints[0].theta0().value() == 2.0


def test_angle_restraint_molecular_container(ala_mols):
    """Tests that angle restraints can be set up correctly with a molecular container (i.e. passing the angle of the molecule)."""
    mols = ala_mols.clone()
    ang = mols.angles()[0]
    restraints = sr.restraints.angle(mols=mols, atoms=ang.atoms())


def test_default_dihedral_restraints_setup(ala_mols):
    """Tests that dihedral restraints can be set up correctly with default parameters."""
    mols = ala_mols.clone()
    restraints = sr.restraints.dihedral(mols=mols, atoms=[0, 1, 2, 3])
    assert restraints.num_restraints() == 1
    assert restraints[0].atoms() == [0, 1, 2, 3]
    assert restraints[0].kphi().value() == 100.0


def test_dihedral_restraint_custom_kphi(ala_mols):
    """Tests that dihedral restraints can be set up correctly with custom kphi."""
    mols = ala_mols.clone()
    restraints = sr.restraints.dihedral(
        mols=mols, atoms=[0, 1, 2, 3], kphi="10 kcal mol-1 rad-2"
    )
    assert restraints.num_restraints() == 1
    assert restraints[0].atoms() == [0, 1, 2, 3]
    assert restraints[0].kphi().value() == 10.0


def test_dihedral_restraint_custom_kphi_and_phi0(ala_mols):
    """Tests that dihedral restraints can be set up correctly with custom kphi and phi0."""
    mols = ala_mols.clone()
    restraints = sr.restraints.dihedral(
        mols=mols, atoms=[0, 1, 2, 3], kphi="10 kcal mol-1 rad-2", phi0="2 rad"
    )
    assert restraints.num_restraints() == 1
    assert restraints[0].atoms() == [0, 1, 2, 3]
    assert restraints[0].kphi().value() == 10.0
    assert restraints[0].phi0().value() == 2.0


def test_dihedral_restraint_molecular_container(ala_mols):
    """Tests that dihedral restraints can be set up correctly with a molecular container (i.e. passing the dihedral of the molecule)."""
    mols = ala_mols.clone()
    dih = mols.dihedrals()[0]
    restraints = sr.restraints.dihedral(mols=mols, atoms=dih.atoms())


def test_multiple_angle_restraints(ala_mols):
    """Tests that multiple angle restraints can be set up correctly."""
    mols = ala_mols.clone()
    ang0 = mols.angles()[0]
    ang1 = mols.angles()[1]
    restraints = sr.restraints.angle(mols=mols, atoms=ang0.atoms())
    restraint1 = sr.restraints.angle(mols=mols, atoms=ang1.atoms())
    restraints.add(restraint1)
    assert restraints.num_restraints() == 2


def test_multiple_dihedral_restraints(ala_mols):
    """Tests that multiple dihedral restraints can be set up correctly."""
    mols = ala_mols.clone()
    dih0 = mols.dihedrals()[0]
    dih1 = mols.dihedrals()[1]
    dih2 = mols.dihedrals()[2]
    restraints = sr.restraints.dihedral(mols=mols, atoms=dih0.atoms())
    restraint1 = sr.restraints.dihedral(mols=mols, atoms=dih1.atoms())
    restraint2 = sr.restraints.dihedral(mols=mols, atoms=dih2.atoms())
    restraints.add(restraint1)
    restraints.add(restraint2)
    assert restraints.num_restraints() == 3
