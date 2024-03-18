import sire as sr

import pytest


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
def test_decouple(ala_mols):
    mol = ala_mols[0]

    amol = sr.morph.annihilate(mol)

    omm_fwds = amol.perturbation().to_openmm(constraint="h-bonds")

    omm_bwds = amol.perturbation().to_openmm(constraint="h-bonds",
                                             swap_end_states=True)

    # use .changed_bonds() etc to validate that the forwards
    # and backwards perturbations are the same

    # also check that parameters are being set equal to zero correctly

    # raise an error until I have written this!
    assert False