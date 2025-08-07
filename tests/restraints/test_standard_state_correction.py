import pytest

# Boresch parameters from old test_boresch_corr SOMD test so we can compare
# to previous results. Note that the force constants are
BORESCH_PARAMS_DEFAULT = {
    "receptor_selection": "atomidx 1538 or atomidx 1518 or atomidx 1540",
    "ligand_selection": "atomidx 4 or atomidx 3 or atomidx 5",
    "kr": "25.49 kcal mol-1 A-2",
    "ktheta": ["66.74 kcal mol-1 rad-2", "38.39 kcal mol-1 rad-2"],
    "kphi": [
        "215.36 kcal mol-1 rad-2",
        "49.23 kcal mol-1 rad-2",
        "49.79 kcal mol-1 rad-2",
    ],
    "r0": "5.92 A",
    "theta0": ["1.85 rad", "1.59 rad"],
    "phi0": ["-0.30 rad", "-1.55 rad", "2.90 rad"],
    "name": "boresch_restraint_1",
}


@pytest.mark.slow
def test_standard_state_correction_boresch(thrombin_complex):
    """
    Check that the parameters of the created Boresch restraint are as expected.
    """
    import sire as sr
    from sire.restraints import boresch, get_standard_state_correction

    boresch_restraints = boresch(
        thrombin_complex,
        receptor=thrombin_complex["protein"][
            BORESCH_PARAMS_DEFAULT["receptor_selection"]
        ],
        ligand=thrombin_complex["resname LIG"][
            BORESCH_PARAMS_DEFAULT["ligand_selection"]
        ],
        kr=BORESCH_PARAMS_DEFAULT["kr"],
        ktheta=BORESCH_PARAMS_DEFAULT["ktheta"],
        kphi=BORESCH_PARAMS_DEFAULT["kphi"],
        r0=BORESCH_PARAMS_DEFAULT["r0"],
        theta0=BORESCH_PARAMS_DEFAULT["theta0"],
        phi0=BORESCH_PARAMS_DEFAULT["phi0"],
        name=BORESCH_PARAMS_DEFAULT["name"],
    )
    boresch_restraint = boresch_restraints[0]

    std_state_correction = get_standard_state_correction(
        boresch_restraint, temperature=298.15 * sr.units.kelvin
    )
    assert std_state_correction.value() == pytest.approx(-10.98, abs=1e-2)
