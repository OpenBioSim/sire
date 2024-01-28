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

    ref_mols = sr.morph.extract_reference(mols)

    mols2 = ref_mols.clone()
    mols2.update(sr.morph.create_from_pertfile(ref_mols[0], neopentane_methane_pert))

    print(mols2.molecules("property is_perturbable"))

    assert len(mols2.molecules("property is_perturbable")) == 1

    assert mols.energy().value() == pytest.approx(mols2.energy().value())


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

    assert mols.energy().value() == pytest.approx(ref_mols.energy().value())

    pert_mols = sr.morph.extract_perturbed(mols)

    with pytest.raises(KeyError):
        assert len(pert_mols.molecules("property is_perturbable")) == 0

    mols = sr.morph.link_to_perturbed(mols)

    assert mols.energy().value() == pytest.approx(pert_mols.energy().value())
