import pytest
import sire as sr


def test_repartition_hydrogen_masses(ala_mols):
    mol = ala_mols.molecules()[0]

    newmol = sr.morph.repartition_hydrogen_masses(mol)

    for atom0, atom1 in zip(mol.atoms(), newmol.atoms()):
        if atom0.element().num_protons() == 1:
            assert atom1.mass() > atom0.mass()
            assert atom1.mass() >= sr.u("4 g mol-1")
        else:
            assert atom1.mass() <= atom0.mass()

    assert mol.mass().value() == pytest.approx(newmol.mass().value())
