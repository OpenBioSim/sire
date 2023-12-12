import sire as sr
import pytest


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_dynamics_return_type(ala_mols):
    mols = ala_mols

    m = mols.dynamics(platform="Reference").run("1 fs").commit()

    assert isinstance(m, type(mols))

    mol = mols[0]

    m = mol.dynamics(platform="Reference").run("1 fs").commit()

    assert isinstance(m, type(mol))

    atom = mol[0]

    a = atom.dynamics(platform="Reference").run("1 fs").commit()

    assert isinstance(a, type(atom))


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_minimisation_return_type(ala_mols):
    mols = ala_mols

    m = mols.minimisation(platform="Reference").run(1).commit()

    assert isinstance(m, type(mols))

    mol = mols[0]

    m = mol.minimisation(platform="Reference").run(1).commit()

    assert isinstance(m, type(mol))

    atom = mol[0]

    a = atom.minimisation(platform="Reference").run(1).commit()

    assert isinstance(a, type(atom))
