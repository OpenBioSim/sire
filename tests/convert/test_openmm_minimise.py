import pytest
import sire as sr


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_openmm_simple_minimise(kigaki_mols):
    mols = kigaki_mols

    nrg = mols.energy()

    mols = mols.minimisation(platform="cpu").run(5).commit()

    nrg2 = mols.energy()

    assert nrg2.value() < nrg.value()


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_openmm_minimise_lambda(merged_ethane_methanol):
    mols = merged_ethane_methanol.clone()

    for mol in mols.molecules("molecule property is_perturbable"):
        mols.update(mol.perturbation().link_to_reference().commit())

    # this blows up because of incompatible exceptions/exclusions
    mol = (
        mols[0].minimisation(platform="cpu", lambda_value=1.0).run(5).commit()
    )

    mols = (
        mols[0].minimisation(platform="cpu", lambda_value=0.0).run(5).commit()
    )

    mols = (
        mols[0].minimisation(platform="cpu", lambda_value=0.5).run(5).commit()
    )
