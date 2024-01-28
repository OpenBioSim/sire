import sire as sr

import pytest

# get the directory of this script file
import os

test_dir = os.path.dirname(os.path.realpath(__file__))

neopentane_methane_pert = os.path.join(test_dir, "input", "neopentane_methane.pert")


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_pertfile(neopentane_methane):
    mols = neopentane_methane.clone()

    mols = sr.morph.link_to_reference(mols)

    ref_mols = sr.morph.extract_reference(mols)

    mols2 = ref_mols.clone()
    mols2.update(sr.morph.create_from_pertfile(ref_mols[0], neopentane_methane_pert))

    assert len(mols2.molecules("property is_perturbable")) == 1

    d = mols.dynamics(lambda_value=0.0)
    d2 = mols2.dynamics(lambda_value=0.0)

    assert d.current_potential_energy().value() == pytest.approx(
        d2.current_potential_energy().value(), 1e-3
    )

    d = mols.dynamics(lambda_value=1.0)
    d2 = mols2.dynamics(lambda_value=1.0)

    assert d.current_potential_energy().value() == pytest.approx(
        d2.current_potential_energy().value(), 1e-3
    )

    d = mols.dynamics(lambda_value=0.5)
    d2 = mols2.dynamics(lambda_value=0.5)

    assert d.current_potential_energy().value() == pytest.approx(
        d2.current_potential_energy().value(), 1e-3
    )


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_extract_and_link(neopentane_methane):
    mols = neopentane_methane.clone()

    assert len(mols.molecules("property is_perturbable")) == 1

    ref_mols = sr.morph.extract_reference(mols)

    with pytest.raises(KeyError):
        assert len(ref_mols.molecules("property is_perturbable")) == 0

    mols = sr.morph.link_to_reference(mols)

    for key in ref_mols[0].property_keys():
        assert mols[0].property(key) == ref_mols[0].property(key)

    d = mols.dynamics(lambda_value=0.0)
    d2 = ref_mols.dynamics()

    assert d.current_potential_energy().value() == pytest.approx(
        d2.current_potential_energy().value()
    )

    pert_mols = sr.morph.extract_perturbed(mols)

    with pytest.raises(KeyError):
        assert len(pert_mols.molecules("property is_perturbable")) == 0

    mols = sr.morph.link_to_perturbed(mols)

    for key in ref_mols[0].property_keys():
        assert mols[0].property(key) == pert_mols[0].property(key)

    d = mols.dynamics(lambda_value=1.0)
    d2 = pert_mols.dynamics()

    assert d.current_potential_energy().value() == pytest.approx(
        d2.current_potential_energy().value()
    )


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_get_pert_info(neopentane_methane):
    mols = neopentane_methane.clone()

    mols = sr.morph.link_to_reference(mols)

    mols = sr.morph.zero_dummy_torsions(mols)

    p = mols[0].perturbation().to_openmm(constraint="bonds")

    atoms = p.changed_atoms(to_pandas=True)

    print(atoms)

    bonds = p.changed_bonds(to_pandas=True)

    print(bonds)

    angles = p.changed_angles(to_pandas=True)

    print(angles)

    torsions = p.changed_torsions(to_pandas=True)

    print(torsions)

    exceptions = p.changed_exceptions(to_pandas=True)

    print(exceptions)

    constraints = p.changed_constraints(to_pandas=True)

    print(constraints)

    assert False
