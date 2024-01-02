import pytest
import sire as sr


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_openmm_simple_minimise(ala_mols, openmm_platform):
    mols = ala_mols

    nrg = mols.energy()

    mols = mols.minimisation(platform=openmm_platform).run(5).commit()

    nrg2 = mols.energy()

    assert nrg2.value() < nrg.value()


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_openmm_minimise_lambda(merged_ethane_methanol, openmm_platform):
    mols = merged_ethane_methanol.clone()

    for mol in mols.molecules("molecule property is_perturbable"):
        mols.update(mol.perturbation().link_to_reference().commit())

    # this blows up because of incompatible exceptions/exclusions
    mol = (
        mols[0].minimisation(platform=openmm_platform, lambda_value=1.0).run(5).commit()
    )

    mols = (
        mols[0].minimisation(platform=openmm_platform, lambda_value=0.0).run(5).commit()
    )

    mols = (
        mols[0].minimisation(platform=openmm_platform, lambda_value=0.5).run(5).commit()
    )


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_openmm_minimise_unbonded_water(kigaki_mols, openmm_platform):
    mols = kigaki_mols

    atoms = mols[100].atoms()

    # the water molecules have no internal bonds, so this tests
    # whether or not constraints have been added correctly
    mols = mols.minimisation(platform=openmm_platform).run(10).commit()

    new_atoms = mols[100].atoms()

    # make sure the water has not blown up... (low tolerance as
    # the water bond lengths will be constrained to a common shared
    # value)
    assert (
        sr.measure(atoms[0], atoms[1]).value()
        == pytest.approx(sr.measure(new_atoms[0], new_atoms[1]).value(), abs=1e-2)
    ) or (sr.measure(new_atoms[0], new_atoms[2]).value() == pytest.approx(0.9572))

    assert (
        sr.measure(atoms[0], atoms[2]).value()
        == pytest.approx(sr.measure(new_atoms[0], new_atoms[2]).value(), abs=1e-2)
    ) or (sr.measure(new_atoms[0], new_atoms[2]).value() == pytest.approx(0.9572))

    assert sr.measure(atoms[1], atoms[2]).value() == pytest.approx(
        sr.measure(new_atoms[1], new_atoms[2]).value(), abs=1e-2
    )


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_openmm_minimise_vacuum(kigaki_mols, openmm_platform):
    mols = kigaki_mols

    mols = mols.minimisation(platform=openmm_platform, vacuum=True).run(10).commit()

    assert not mols.property("space").is_periodic()

    mols = kigaki_mols.clone()

    mols.add_shared_property("space", sr.vol.Cartesian())

    mols = mols.minimisation(platform=openmm_platform, vacuum=True).run(10).commit()

    assert not mols.property("space").is_periodic()
