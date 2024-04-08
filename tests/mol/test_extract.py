import sire as sr

import pytest


def test_extract(neura_mols):
    protein = neura_mols["protein"]

    # Extract the middle 10 residues
    middle_residues = protein.residues()[10:20]

    extracted_middle_residues = middle_residues.extract().molecule()

    assert middle_residues.num_atoms() == extracted_middle_residues.num_atoms()

    assert middle_residues.num_residues() == extracted_middle_residues.num_residues()

    # Extract the third of these residues
    res_0 = middle_residues[2].extract().molecule()
    res_1 = extracted_middle_residues.residues()[2].extract().molecule()

    assert res_0.num_atoms() == middle_residues[2].num_atoms()
    assert res_0.num_residues() == 1

    assert res_0.num_atoms() == res_1.num_atoms()
    assert res_0.num_residues() == res_1.num_residues()

    # this validates that most of the properties have been extracted
    # properly
    assert res_0.energy().value() == pytest.approx(res_1.energy().value())
    assert res_0.energy().value() == pytest.approx(middle_residues[2].energy().value())
    assert extracted_middle_residues.energy().value() == pytest.approx(
        middle_residues.energy().value()
    )
