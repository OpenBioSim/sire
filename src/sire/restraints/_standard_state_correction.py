__all__ = ["get_standard_state_correction"]

import numpy as _np

import sire as sr

from .. import u as _u
from .. import units as _units


def get_standard_state_correction(restraint, temperature=_u("300 K")):
    """
    Get the entropic correction for releasing a given restraint to
    the standard state volume at a given temperature.

    Parameters
    ----------
    restraint : sire.legacy.MM._MM.BoreschRestraint
        The restraint for which to calculate the standard state correction.

    temperature : sire.legacy.Units._Units.GeneralUnit, optional
        The temperature at which to calculate the standard state correction.

    Returns
    -------
    correction : sire.legacy.Units._Units.GeneralUnit
        The standard state correction for the given restraint.

    Examples
    --------
    To create a Boresch restraint and calculate the standard state correction:

    >>> # Create the Boresch restraints object.
    >>> my_boresch_restraints = boresch(
    >>>     system,
    >>>     receptor=system["protein"]["atomidx 1538 or atomidx 1518 or atomidx 1540"],
    >>>     ligand=system["resname LIG"]["atomidx 4 or atomidx 3 or atomidx 5"],
    >>>     kr="6.2012 kcal mol-1 A-2",
    >>>     ktheta=["28.7685 kcal mol-1 rad-2", "24.8204 kcal mol-1 rad-2"],
    >>>     kphi=["59.8626 kcal mol-1 rad-2", "0.7923 kcal mol-1 rad-2", "55.1775 kcal mol-1 rad-2"],
    >>>     r0="7.687 A",
    >>>     theta0=["1.3031 rad", "1.4777 rad"],
    >>>     phi0=["2.5569 rad", "2.9359 rad", "1.4147 rad"],
    >>> )
    >>>
    >>> # Extract the single Boresch restraint from the Boresch restraints.
    >>> my_boresch_restraint = my_boresch_restraints[0]
    >>>
    >>> # Calculate the standard state correction for the Boresch restraint.
    >>> standard_state_correction = get_standard_state_correction(my_boresch_restraint)
    >>> print(standard_state_correction)
    """
    if isinstance(restraint, sr.legacy.MM._MM.BoreschRestraint):
        return _get_boresch_standard_state_correction(restraint, temperature)
    else:
        raise NotImplementedError(
            f"Standard state correction for restraint type {type(restraint)} is not implemented. "
            "This function currently only supports Boresch restraints."
        )


def _get_boresch_standard_state_correction(restraint, temperature):
    """
    Get the entropic correction for releasing a given Boresch restraint to
    the standard state volume at a given temperature.

    Parameters
    ----------
    restraint : sire.legacy.MM._MM.BoreschRestraint
        The Boresch restraint for which to calculate the standard state correction.
    temperature : float
        The temperature at which to calculate the standard state correction.

    Returns
    -------
    correction : sire.legacy.Units._Units.GeneralUnit
        The standard state correction for the given Boresch restraint.

    Examples
    --------
    To create a Boresch restraint and calculate the standard state correction:
    ```
    # Create the Boresch restraints object.
    my_boresch_restraints = boresch(
        system,
        receptor=system["protein"]["atomidx 1538 or atomidx 1518 or atomidx 1540"],
        ligand=system["resname LIG"]["atomidx 4 or atomidx 3 or atomidx 5"],
        kr="6.2012 kcal mol-1 A-2",
        ktheta=["28.7685 kcal mol-1 rad-2", "24.8204 kcal mol-1 rad-2"],
        kphi=["59.8626 kcal mol-1 rad-2", "0.7923 kcal mol-1 rad-2", "55.1775 kcal mol-1 rad-2"],
        r0="7.687 A",
        theta0=["1.3031 rad", "1.4777 rad"],
        phi0=["2.5569 rad", "2.9359 rad", "1.4147 rad"],
    )

    # Extract the single Boresch restraint from the Boresch restraints.
    my_boresch_restraint = my_boresch_restraints[0]

    # Calculate the standard state correction for the Boresch restraint.
    standard_state_correction = get_standard_state_correction(my_boresch_restraint)
    print(standard_state_correction)
    ```
    """
    # Remove units so that we can raise to non-integer powers and take sine.

    # Params.
    T = (temperature / _units.kelvin).value()  # K
    r0 = (restraint.r0() / _units.angstrom).value()  # A
    thetaA0 = (
        restraint.theta0()[0] / _units.radians
    ).value()  # Remove units so we can take sin.
    thetaB0 = (
        restraint.theta0()[1] / _units.radians
    ).value()  # Remove units so we can take sin.

    # Constants.
    v0 = (
        (_units.meter3 / 1000) / _units.mole.value()
    ).value()  # A^3, the standard state volume.
    R = _units.gasr.value()  # kcal mol-1, the molar gas constant.
    prefactor = 8 * (_np.pi**2) * v0  # Divide this to account for force constants of 0.
    force_constants = []

    # Correct for force constants of zero which break the analytical correction.
    # kr
    kr = (
        restraint.kr() / (_units.kcal * _units.mole.pow(-1) * _units.angstrom.pow(-2))
    ).value()
    force_constants.append(kr)

    # kthetas
    for i in range(2):
        if restraint.ktheta()[i] == 0:
            prefactor /= 2 / _np.sin((restraint.theta0()[i] / _units.radians).value())
        else:
            force_const = (
                restraint.ktheta()[i]
                / (_units.kcal * _units.mole.pow(-1) * _units.radians.pow(-2))
            ).value()
            force_constants.append(force_const)

    # kphis
    for i in range(3):
        if restraint.kphi()[i] == 0:
            prefactor /= 2 * _np.pi
        else:
            force_const = (
                restraint.kphi()[i]
                / (_units.kcal * _units.mole.pow(-1) * _units.radians.pow(-2))
            ).value()
            force_constants.append(force_const)

    n_nonzero_k = len(force_constants)

    # Calculation
    numerator = prefactor * _np.sqrt(_np.prod(force_constants))
    denominator = (
        r0**2
        * _np.sin(thetaA0)
        * _np.sin(thetaB0)
        * (2 * _np.pi * R * T) ** (n_nonzero_k / 2)
    )

    # Use values with units to return a sire.legacy.Units._Units.GeneralUnit object.
    dg = -_units.gasr * temperature * _np.log(numerator / denominator)
    return dg
