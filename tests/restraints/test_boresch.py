import pytest

import sire as sr
from sire.restraints import boresch

# Valid Boresch restraint parameters.
BORESCH_PARAMS_DEFAULT = {
    "receptor_selection": "atomidx 1538 or atomidx 1518 or atomidx 1540",
    "ligand_selection": "atomidx 4 or atomidx 3 or atomidx 5",
    "kr": "6.2012 kcal mol-1 A-2",
    "ktheta": ["28.7685 kcal mol-1 rad-2", "24.8204 kcal mol-1 rad-2"],
    "kphi": [
        "59.8626 kcal mol-1 rad-2",
        "0.7923 kcal mol-1 rad-2",
        "55.1775 kcal mol-1 rad-2",
    ],
    "r0": "7.687 A",
    "theta0": ["1.3031 rad", "1.4777 rad"],
    "phi0": ["2.5569 rad", "2.9359 rad", "1.4147 rad"],
    "name": "boresch_restraint_1",
}
BORESCH_PARAMS_UNNAMED = BORESCH_PARAMS_DEFAULT.copy()
BORESCH_PARAMS_UNNAMED["name"] = None
BORESCH_PARAMS_NON_EXPLICIT = BORESCH_PARAMS_DEFAULT.copy()
for param in ("kr", "ktheta", "kphi", "r0", "theta0", "phi0"):
    BORESCH_PARAMS_NON_EXPLICIT[param] = None
BORESCH_PARAMS_SINGLE_FORCE_CONSTS = BORESCH_PARAMS_DEFAULT.copy()
BORESCH_PARAMS_SINGLE_FORCE_CONSTS["ktheta"] = ["28.7685 kcal mol-1 rad-2"]
BORESCH_PARAMS_SINGLE_FORCE_CONSTS["kphi"] = ["59.8626 kcal mol-1 rad-2"]

# Invalid Boresch restraint parameters.
BORESCH_PARAMS_2_RECEPTOR_ATOMS = BORESCH_PARAMS_DEFAULT.copy()
BORESCH_PARAMS_2_RECEPTOR_ATOMS["receptor_selection"] = "atomidx 1538 or atomidx 1518"
BORESCH_PARAMS_2_LIGAND_ATOMS = BORESCH_PARAMS_DEFAULT.copy()
BORESCH_PARAMS_2_LIGAND_ATOMS["ligand_selection"] = "atomidx 4 or atomidx 3"
BORESCH_PARAMS_WRONG_NUM_KPHI = BORESCH_PARAMS_DEFAULT.copy()
BORESCH_PARAMS_WRONG_NUM_KPHI["kphi"] = [
    "28.7685 kcal mol-1 rad-2",
    "24.8204 kcal mol-1 rad-2",
]
BORESCH_PARAMS_KR_0 = BORESCH_PARAMS_DEFAULT.copy()
BORESCH_PARAMS_KR_0["kr"] = "0.0 kcal mol-1 A-2"
BORESCH_PARAMS_KTHETA_A_0 = BORESCH_PARAMS_DEFAULT.copy()
BORESCH_PARAMS_KTHETA_A_0["ktheta"] = [
    "0.0 kcal mol-1 rad-2",
    "24.8204 kcal mol-1 rad-2",
]
BORESCH_PARAMS_KTHETA_B_0 = BORESCH_PARAMS_DEFAULT.copy()
BORESCH_PARAMS_KTHETA_B_0["ktheta"] = [
    "28.7685 kcal mol-1 rad-2",
    "0.0 kcal mol-1 rad-2",
]

# Restraints likely to be unstable due to thetaA0 or thetaB0 being too close to 0 or pi
# or the associated force constants being too low.
BORESCH_PARAMS_UNSTABLE_FORCE_CONSTS_LOW = BORESCH_PARAMS_DEFAULT.copy()
BORESCH_PARAMS_UNSTABLE_FORCE_CONSTS_LOW["ktheta"] = [
    "0.1 kcal mol-1 rad-2",
    "24.8204 kcal mol-1 rad-2",
]
BORESCH_PARAMS_UNSTABLE_THETA_A_0 = BORESCH_PARAMS_DEFAULT.copy()
BORESCH_PARAMS_UNSTABLE_THETA_A_0["theta0"] = ["0.1 rad", "1.4777 rad"]
BORESCH_PARAMS_UNSTABLE_THETA_B_0 = BORESCH_PARAMS_DEFAULT.copy()
BORESCH_PARAMS_UNSTABLE_THETA_B_0["theta0"] = ["1.3031 rad", "3.1 rad"]


# Parameterise the test with many combinations of parameters which should pass.
@pytest.mark.parametrize(
    (
        "receptor_selection",
        "ligand_selection",
        "kr",
        "ktheta",
        "kphi",
        "r0",
        "theta0",
        "phi0",
        "name",
    ),
    [
        [val for val in BORESCH_PARAMS_DEFAULT.values()],
        [val for val in BORESCH_PARAMS_UNNAMED.values()],
        [val for val in BORESCH_PARAMS_NON_EXPLICIT.values()],
        [val for val in BORESCH_PARAMS_SINGLE_FORCE_CONSTS.values()],
    ],
)
def test_create_boresch_restraint(
    thrombin_complex,
    receptor_selection,
    ligand_selection,
    kr,
    ktheta,
    kphi,
    r0,
    theta0,
    phi0,
    name,
):
    """
    Check that we can create a Boresch restraint with valid parameters.
    """
    boresch_restraints = boresch(
        thrombin_complex,
        receptor=thrombin_complex["protein"][receptor_selection],
        ligand=thrombin_complex["resname LIG"][ligand_selection],
        kr=kr,
        ktheta=ktheta,
        kphi=kphi,
        r0=r0,
        theta0=theta0,
        phi0=phi0,
        name=name,
    )


# Parameterise the test with many combinations of parameters which should fail.
@pytest.mark.parametrize(
    (
        "receptor_selection",
        "ligand_selection",
        "kr",
        "ktheta",
        "kphi",
        "r0",
        "theta0",
        "phi0",
        "name",
    ),
    [
        [val for val in BORESCH_PARAMS_2_RECEPTOR_ATOMS.values()],
        [val for val in BORESCH_PARAMS_2_LIGAND_ATOMS.values()],
        [val for val in BORESCH_PARAMS_WRONG_NUM_KPHI.values()],
        [val for val in BORESCH_PARAMS_KR_0.values()],
        [val for val in BORESCH_PARAMS_KTHETA_A_0.values()],
        [val for val in BORESCH_PARAMS_KTHETA_B_0.values()],
    ],
)
def test_create_boresch_restraint_fails(
    thrombin_complex,
    receptor_selection,
    ligand_selection,
    kr,
    ktheta,
    kphi,
    r0,
    theta0,
    phi0,
    name,
):
    """
    Check that we can create a Boresch restraint with valid parameters.
    """
    with pytest.raises(ValueError):
        boresch_restraints = boresch(
            thrombin_complex,
            receptor=thrombin_complex["protein"][receptor_selection],
            ligand=thrombin_complex["resname LIG"][ligand_selection],
            kr=kr,
            ktheta=ktheta,
            kphi=kphi,
            r0=r0,
            theta0=theta0,
            phi0=phi0,
            name=name,
        )


def test_boresch_restraint_params(thrombin_complex):
    """
    Check that the parameters of the created Boresch restraint are as expected.
    """
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

    assert boresch_restraint.receptor_atoms() == [1574, 1554, 1576]
    assert boresch_restraint.ligand_atoms() == [4, 3, 5]
    assert boresch_restraint.kr().value() == 6.2012
    assert boresch_restraint.ktheta()[0].value() == 28.7685
    assert boresch_restraint.ktheta()[1].value() == 24.8204
    assert boresch_restraint.kphi()[0].value() == 59.8626
    assert boresch_restraint.kphi()[1].value() == 0.7923
    assert boresch_restraint.kphi()[2].value() == 55.1775
    assert boresch_restraint.r0().value() == 7.687
    assert boresch_restraint.theta0()[0].value() == 1.3031
    assert boresch_restraint.theta0()[1].value() == 1.4777
    assert boresch_restraint.phi0()[0].value() == 2.5569
    assert boresch_restraint.phi0()[1].value() == 2.9359
    assert boresch_restraint.phi0()[2].value() == 1.4147


@pytest.mark.parametrize(
    (
        "receptor_selection",
        "ligand_selection",
        "kr",
        "ktheta",
        "kphi",
        "r0",
        "theta0",
        "phi0",
        "name",
    ),
    [
        [val for val in BORESCH_PARAMS_UNSTABLE_FORCE_CONSTS_LOW.values()],
        [val for val in BORESCH_PARAMS_UNSTABLE_THETA_A_0.values()],
        [val for val in BORESCH_PARAMS_UNSTABLE_THETA_B_0.values()],
    ],
)
def test_boresch_restraint_params_unstable(
    thrombin_complex,
    receptor_selection,
    ligand_selection,
    kr,
    ktheta,
    kphi,
    r0,
    theta0,
    phi0,
    name,
):
    """
    Check that a warning is raised when creating a Boresch restraint with unstable parameters.
    """
    with pytest.warns(UserWarning):
        boresch_restraints = boresch(
            thrombin_complex,
            receptor=thrombin_complex["protein"][receptor_selection],
            ligand=thrombin_complex["resname LIG"][ligand_selection],
            kr=kr,
            ktheta=ktheta,
            kphi=kphi,
            r0=r0,
            theta0=theta0,
            phi0=phi0,
            name=name,
        )
