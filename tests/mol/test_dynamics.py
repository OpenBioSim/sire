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


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_cutoff_options(ala_mols):
    mols = ala_mols

    d = mols.dynamics(platform="Reference", cutoff="10 A")

    assert d.info().cutoff() == sr.u("10 A")

    d = mols[0].dynamics(platform="Reference", cutoff="infinite", vacuum=True)

    # OpenMM cannot have no cutoff, so sets it to a very large number
    assert d.info().cutoff() > sr.u("1000 A")

    d = mols[0].dynamics(platform="Reference", cutoff="none", vacuum=True)

    # OpenMM cannot have no cutoff, so sets it to a very large number
    assert d.info().cutoff() > sr.u("1000 A")

    d = mols.dynamics(platform="Reference", cutoff=sr.u("7.5A"), cutoff_type="PME")

    assert d.info().cutoff() == sr.u("7.5 A")
    assert d.info().cutoff_type() == "PME"
