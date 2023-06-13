__all__ = [
    "angle",
    "convert",
    "length",
    "clear_default_units",
    "set_default_unit",
    "set_default_units",
    "set_energy_unit",
    "set_internal_units",
    "set_length_unit",
    "set_mass_unit",
    "set_quantity_unit",
    "set_si_units",
    "set_time_unit",
    "GeneralUnit",
    "akma_time",
    "akma_velocity",
    "amp",
    "angle_minute",
    "angle_minutes",
    "angle_second",
    "angle_seconds",
    "angstrom",
    "angstroms",
    "angstrom2",
    "angstrom3",
    "angstroms_per_fs",
    "atm",
    "atomic_mass_constant",
    "bar",
    "bohr_radii",
    "c",
    "cal",
    "cal_per_mol",
    "celsius",
    "centimeter",
    "centimeters",
    "circumference",
    "convert",
    "coulomb",
    "coulomb_per_mol",
    "day",
    "degree",
    "degrees",
    "e_charge",
    "electron_mass",
    "epsilon0",
    "fahrenheit",
    "farad",
    "faraday",
    "feet",
    "femtogram",
    "femtosecond",
    "fg_per_mol",
    "fortnight",
    "four_pi_eps0",
    "g_accel",
    "g_per_mol",
    "gasr",
    "grad",
    "gradian",
    "gradians",
    "gram",
    "gon",
    "G_newton",
    "h_planck",
    "h_slash",
    "hartree",
    "hour",
    "hundredweight",
    "inch",
    "int_cal",
    "int_cal_per_mol",
    "int_kcal",
    "int_kcal_per_mol",
    "joule",
    "J_per_mol",
    "k_boltz",
    "kcal",
    "kcal_per_mol",
    "kelvin",
    "kilogram",
    "kilojoule",
    "kilometer",
    "kilometers",
    "kilometers_per_hour",
    "kg_per_mol",
    "kJ_per_mol",
    "kph",
    "megajoule",
    "meter",
    "meter2",
    "meter3",
    "meters",
    "meters_per_second",
    "mg_per_mol",
    "microgram",
    "micrometer",
    "microsecond",
    "miles",
    "miles_per_hour",
    "milligram",
    "millimeter",
    "millimeters",
    "millisecond",
    "minute",
    "mmHg",
    "mod_electron",
    "molar_volume",
    "mole",
    "mph",
    "mu0",
    "MJ_per_mol",
    "nanogram",
    "nanometer",
    "nanometers",
    "nanometer2",
    "nanosecond",
    "neutron_mass",
    "newton",
    "ng_per_mol",
    "octant",
    "octants",
    "one_over_four_pi_eps0",
    "ounce",
    "pascal",
    "pg_per_mol",
    "picogram",
    "picometer",
    "picometers",
    "picometer2",
    "picometer3",
    "picosecond",
    "pound",
    "proton_mass",
    "psi",
    "quadrant",
    "quadrants",
    "radian",
    "radians",
    "revolution",
    "revolutions",
    "revs",
    "second",
    "sextant",
    "sextants",
    "stone",
    "tonne",
    "tonne_per_mol",
    "ug_per_mol",
    "volt",
    "watt",
    "watt_per_mol",
    "week",
    "yards",
]

from ..legacy import Units as _Units

from .. import use_new_api as _use_new_api

from ..legacy.Units import (
    akma_time,
    akma_velocity,
    amp,
    angle_minute,
    angle_minutes,
    angle_second,
    angle_seconds,
    angstrom,
    angstroms,
    angstrom2,
    angstrom3,
    angstroms_per_fs,
    atm,
    atomic_mass_constant,
    bar,
    bohr_radii,
    c,
    cal,
    cal_per_mol,
    celsius,
    centimeter,
    centimeters,
    circumference,
    convert,
    coulomb,
    coulomb_per_mol,
    day,
    degree,
    degrees,
    e_charge,
    electron_mass,
    epsilon0,
    fahrenheit,
    farad,
    faraday,
    feet,
    femtogram,
    femtosecond,
    fg_per_mol,
    fortnight,
    four_pi_eps0,
    g_accel,
    g_per_mol,
    gasr,
    grad,
    gradian,
    gradians,
    gram,
    gon,
    G_newton,
    h_planck,
    h_slash,
    hartree,
    hour,
    hundredweight,
    inch,
    int_cal,
    int_cal_per_mol,
    int_kcal,
    int_kcal_per_mol,
    joule,
    J_per_mol,
    k_boltz,
    kcal,
    kcal_per_mol,
    kelvin,
    kilogram,
    kilojoule,
    kilometer,
    kilometers,
    kilometers_per_hour,
    kg_per_mol,
    kJ_per_mol,
    kph,
    megajoule,
    meter,
    meter2,
    meter3,
    meters,
    meters_per_second,
    mg_per_mol,
    microgram,
    micrometer,
    microsecond,
    miles,
    miles_per_hour,
    milligram,
    millimeter,
    millimeters,
    millisecond,
    minute,
    mmHg,
    mod_electron,
    molar_volume,
    mole,
    mph,
    mu0,
    MJ_per_mol,
    nanogram,
    nanometer,
    nanometers,
    nanometer2,
    nanosecond,
    neutron_mass,
    newton,
    ng_per_mol,
    octant,
    octants,
    one_over_four_pi_eps0,
    ounce,
    pascal,
    pg_per_mol,
    picogram,
    picometer,
    picometers,
    picometer2,
    picometer3,
    picosecond,
    pound,
    proton_mass,
    psi,
    quadrant,
    quadrants,
    radian,
    radians,
    revolution,
    revolutions,
    revs,
    second,
    sextant,
    sextants,
    stone,
    tonne,
    tonne_per_mol,
    ug_per_mol,
    volt,
    watt,
    watt_per_mol,
    week,
    yards,
    GeneralUnit,
)

_use_new_api()


def _fix_generalunit():
    def _generalunit_approx_equal(u, v):
        if not hasattr(u, "what"):
            u = GeneralUnit(u)

        if not hasattr(v, "what"):
            v = GeneralUnit(v)

        if u.what().endswith("Property"):
            return _generalunit_approx_equal(u.value(), v)
        elif v.what().endswith("Property"):
            return _generalunit_approx_equal(u, v.value())

        # make sure that the units are the same
        if u.has_same_units(v):
            from ..search import approx_equal

            return approx_equal(u.value(), v.value())
        else:
            return False

    def _generalunit_to_default(obj):
        """Return a floating point value that represents
        this value in default units for this dimension

        Example
        -------

        >>> import sire as sr
        >>> l = 5 * sr.units.angstrom
        >>> sr.units.set_length_unit(sr.units.picometer)
        >>> print(l.to_default())
            500.0
        """
        return obj.value() / obj.get_default().value()

    def __generalunit__getitem__(obj, key):
        return obj.get_component(key)

    def __generalunit__setitem__(obj, key, value):
        obj.set_component(key, value)
        return obj

    GeneralUnit.approx_equal = _generalunit_approx_equal
    GeneralUnit.to_default = _generalunit_to_default
    GeneralUnit.__setitem__ = __generalunit__setitem__
    GeneralUnit.__getitem__ = __generalunit__getitem__

    def __generalunit__bool__(obj):
        return not obj.is_zero()

    def __generalunit__float__(obj):
        if not obj.is_dimensionless():
            raise TypeError(
                f"You cannot convert the dimensioned value {obj} "
                "to a dimensionless floating point number."
            )

        return obj.value()

    def __generalunit__int__(obj):
        if not obj.is_dimensionless():
            raise TypeError(
                f"You cannot convert the dimensioned value {obj} "
                "to a dimensionless whole number integer."
            )

        return int(obj.value())

    GeneralUnit.__bool__ = __generalunit__bool__
    GeneralUnit.__float__ = __generalunit__float__
    GeneralUnit.__int__ = __generalunit__int__


if not hasattr(GeneralUnit, "to_default"):
    _fix_generalunit()


_names = None


def _get_unit_name(unit):
    global _names

    if _names is None:
        _names = {
            mole: "mol",
            cal: "cal",
            joule: "joule",
            int_cal: "int_cal",
            kcal: "kcal",
            kilojoule: "kJ",
            int_kcal: "int_kcal",
            angstrom: "Å",
            picometer: "pm",
            nanometer: "nm",
            gram: "g",
            kilogram: "kg",
            picosecond: "ps",
            nanosecond: "ns",
            femtosecond: "fs",
        }

    try:
        return _names[unit]
    except Exception:
        pass

    return unit._to_cpp_type()


def _incompatible_units(a, b, typ: str):
    raise TypeError(
        f"Unit {a.unit_string()} is not a {typ}, and so is "
        f"incompatible with {b.unit_string()}."
    )


def set_quantity_unit(q, name: str = None):
    """
    Set the default quantity unit to be used for
    output and default conversions
    """

    if not q.has_same_units(mole):
        _incompatible_units(q, mole, "quantity")

    if name is None:
        name = _get_unit_name(q)

    q.set_as_default(name)


def set_energy_unit(energy, name: str = None):
    """
    Set the default energy unit to be used for
    output and default conversions
    """

    if not energy.has_same_units(kcal):
        _incompatible_units(energy, kcal, "energy")

    if name is None:
        name = _get_unit_name(energy)

    energy.set_as_default(name)


def set_length_unit(length, name: str = None):
    """
    Set the default length unit to be used for
    output and default conversions
    """

    if not length.has_same_units(angstrom):
        _incompatible_units(length, angstrom, "length")

    if name is None:
        name = _get_unit_name(length)

    length.set_as_default(name)


def set_mass_unit(mass, name: str = None):
    """
    Set the default mass unit to be used for
    output and default conversions
    """

    if not mass.has_same_units(gram):
        _incompatible_units(mass, gram, "mass")

    if name is None:
        name = _get_unit_name(mass)

    mass.set_as_default(name)


def set_time_unit(time, name: str = None):
    """
    Set the default time unit to be used for
    output and default conversions
    """

    if not time.has_same_units(picosecond):
        _incompatible_units(time, picosecond, "time")

    if name is None:
        name = _get_unit_name(time)

    time.set_as_default(name)


def clear_default_units():
    """
    Clear the current set of default units. Note that you will need
    to set some units, or else the code will complain.
    """
    GeneralUnit.clear_defaults()


def set_default_units(units):
    """
    Set the default unit names for the units named in the passed
    list of strings.

    For example, set_default_units(["m"]) would set the default length
    unit to "meters", and will print it out using "m" as the string.

    Calling set_default_unit(["kcal", "kcal.mol-1"]) would set the default
    molar energy unit to kcal_per_mol, and will print it out
    using "kcal.mol-1", and the default energy unit to "kcal", printed
    out using "kcal".

    All of the strings must be something that can be parsed by
    the units grammar in GeneralUnit, and should be units only -
    they cannot contain values.
    """
    if type(units) is not list:
        units = [units]

    unit_strings = {}

    # parse all the units first
    for unit in units:
        unit_strings[GeneralUnit(1.0, unit)] = unit

    # we have fully parsed units, and only have one of each type
    for key in unit_strings.keys():
        if hasattr(type(key), "set_as_default"):
            key.set_as_default(unit_strings[key])
        else:
            key.setAsDefault(unit_strings[key])


def set_default_unit(unit: str):
    """
    Set the default unit name for the unit named in the passed
    string.

    For example, set_default_unit("m") would set the default length
    unit to "meters", and will print it out using "m" as the string.

    Calling set_default_unit("kcal.mol-1") would set the default
    molar energy unit to kcal_per_mol, and will print it out
    using "kcal.mol-1"

    The string must be something that can be parsed by
    the units grammar in GeneralUnit, and should be units only -
    they cannot contain values.
    """
    if type(unit) is list:
        set_default_units(unit)
    else:
        set_default_units([unit])


def set_si_units(scaled: bool = True):
    """
    Switch over to using SI units for output and default conversions

    If 'scaled' is true, then units scaled appropriately
    for atomic-level simulations are chosen.

    In this case, we would use

        mass:        gram (g)
        energy:      kilojoule (J)
        length:      nanometer
        time:        picosecond
        quantity:    mole
        charge:      coulomb (C)
        angle:       radian
        temperature: kelvin (K)

    Unscaled, we would use

        mass:        kilogram (g)
        energy:      joule (J)
        length:      meter
        time:        second
        quantity:    mole
        charge:      coulomb (C)
        angle:       radian
        temperature: kelvin (K)
    """
    try:
        GeneralUnit.clear_defaults()
    except AttributeError:
        # this happens on load before the new api has been activated
        pass

    if scaled:
        set_default_units(
            [
                "g",
                "kJ",
                "nm",
                "ps",
                "mol",
                "C",
                "rad",
                "K",
                "kJ mol-1",
                "kJ mol-1 nm-1",
                "kJ mol-1 nm-2",
                "kJ mol-1 nm-3",
                "kJ mol-1 rad-1",
                "kJ mol-1 rad-2",
                "kJ mol-1 K-1",
                "kJ mol-1 K-2",
                "J s",
                "W",
                "J mol-1 s",
                "W mol-1",
                "N",
            ]
        )
    else:
        set_default_units(
            [
                "kg",
                "J",
                "m",
                "s",
                "mol",
                "C",
                "rad",
                "K",
                "J mol-1",
                "J mol-1 nm-1",
                "J mol-1 nm-2",
                "J mol-1 nm-3",
                "J mol-1 rad-1",
                "J mol-1 rad-2",
                "J mol-1 K-1",
                "J mol-1 K-2",
                "J s",
                "W",
                "J mol-1 s",
                "W mol-1",
                "N",
            ]
        )


def set_internal_units():
    """
    Switch over to using default (AKMA-style) units for output
    and default conversions.

    This uses:

        mass:     gram (g)
        energy:   kilocalorie (kcal)
        length:   angstrom
        time:     picosecond
        quantity: mole
        charge:   mod_electron charges
    """
    try:
        GeneralUnit.clear_defaults()
    except AttributeError:
        # this happens on load before the new api has been activated
        pass

    set_default_units(
        [
            "g",
            "kcal",
            "Å",
            "ps",
            "|e|",
            "°",
            "K",
            "kcal mol-1",
            "kcal mol-1 Å-1",
            "kcal mol-1 Å-2",
            "kcal mol-1 Å-3",
            "kcal mol-1 °-1",
            "kcal mol-1 °-2",
            "kcal mol-1 K-1",
            "kcal mol-1 K-2",
            "kcal s",
            "kcal s-1",
            "kcal mol-1 s",
            "kcal mol-1 s-1",
            "kcal Å-1",
        ]
    )


def length(length):
    """Convert the passed argument into a length. This will
    return `length` if it is already a length, or it will
    convert it into a length with default units if this
    is a float or something that can be converted to
    a float
    """
    try:
        return float(length) * angstrom.get_default()
    except Exception:
        pass

    if type(length) != type(angstrom):
        raise TypeError(
            f"The value '{length}' of type {type(length)} is "
            "not a type that is compatible with a GeneralUnit Length"
        )

    return length


def angle(angle):
    """Convert the passed argument into an angle. This will
    return 'angle' if it is already an angle, or will
    convert it into an angle with default units if this
    is a float or something that can be converted to
    a float
    """
    try:
        return float(angle) * degrees.get_default()
    except Exception:
        pass

    if type(angle) != type(degrees):
        raise TypeError(
            f"The value '{angle}' of type {type(angle)} is "
            "not a type that is compatible with a GeneralUnit angle"
        )

    return angle


# set internal units on load
set_internal_units()
