import pytest
import sire as sr

try:
    import pint

    ureg = pint.UnitRegistry()
except Exception:
    ureg = None

U = sr.units.GeneralUnit


@pytest.mark.skipif(ureg is None, reason="pint is not installed")
@pytest.mark.parametrize(
    "text, expect",
    [
        ("10 A", "10 amp"),
        ("10 Å", "10 A"),
        ("3.1e5 kcal mol**-1", "3.1e5 kcal mol-1"),
        ("5 m ** 2 / s ** 2", "5 m2.s-2"),
        ("-3 * meter/second**2", "-3 m s-2"),
    ],
)
def test_parse_pint(text, expect):
    p = sr.u(ureg(text))
    s = sr.u(expect)

    assert p.has_same_units(s)
    assert p.value() == pytest.approx(s.value())


@pytest.mark.parametrize(
    "text, expect",
    [
        ("0", U(0)),
        ("5", U(5)),
        ("-3.2", U(-3.2)),
        ("-3.141e-10", U(-3.141e-10)),
        ("1 angstrom", 1 * sr.units.angstrom),
        ("1 angstrom /* block comment */", 1 * sr.units.angstrom),
        ("1 angstrom // line comment", 1 * sr.units.angstrom),
        ("10 Å", 10 * sr.units.angstrom),
        ("-10.5e5 A", -10.5e5 * sr.units.angstrom),
        ("10 s-1", 10 / sr.units.second),
        ("6.5 m s-1", 6.5 * sr.units.meter / sr.units.second),
        ("1.2 m.s-1", 1.2 * sr.units.meter / sr.units.second),
        (
            "3.2m2.s-2",
            3.2 * sr.units.meter * sr.units.meter / (sr.units.second * sr.units.second),
        ),
        ("55 mm", 55 * sr.units.millimeter),
        ("15 ps", 15 * sr.units.picosecond),
        ("20*angstroms", 20 * sr.units.angstrom),
        ("20 angstrom", 20 * sr.units.angstrom),
        ("120 µm", 120 * sr.units.micrometer),
        (
            "120 μm",
            120 * sr.units.micrometer,
        ),  # this is a different micro symbol
        ("15 * kcal.mol-1", 15 * sr.units.kcal_per_mol),
        (
            "1 kJ/(mol nm**2)",
            1 * sr.units.kJ_per_mol / (sr.units.nanometer * sr.units.nanometer),
        ),
        ("-34.5 kilocalorie mole^-1", -34.5 * sr.units.kcal_per_mol),
        (
            "1.2 kcal mol-1 A-2",
            1.2 * sr.units.kcal_per_mol / (sr.units.angstrom * sr.units.angstrom),
        ),
        (
            "8 (kcal per mol) / angstrom**2",
            8 * sr.units.kcal_per_mol / (sr.units.angstrom * sr.units.angstrom),
        ),
        (
            "8 kcal per mol / angstrom**2",
            8 * sr.units.kcal_per_mol / (sr.units.angstrom * sr.units.angstrom),
        ),
        ("1 fs", 1 * sr.units.femtosecond),
        ("1e-3 femtoseconds", 1e-3 * sr.units.femtosecond),
        (
            "5 (kJ per mol) / nm**2",
            5 * sr.units.kJ_per_mol / (sr.units.nanometer * sr.units.nanometer),
        ),
        (
            "5 kJ per mol / nm**2",
            5 * sr.units.kJ_per_mol / (sr.units.nanometer * sr.units.nanometer),
        ),
        ("5m", 5 * sr.units.meter),
        ("5mm", 5 * sr.units.millimeter),
        ("2 meter / second", 2 * sr.units.meter / sr.units.second),
        ("2.54 centimeter", 2.54 * sr.units.centimeter),
        ("2.54cm", 2.54 * sr.units.centimeter),
        ("15 s2", 15 * sr.units.second * sr.units.second),
        ("-0.45e10s-2", -0.45e10 / (sr.units.second * sr.units.second)),
        (
            "9.8 m.s-1",
            9.8 * sr.units.meter / (sr.units.second),
        ),
        (
            "9.8 meters/second**1",
            9.8 * sr.units.meter / (sr.units.second),
        ),
        (
            "9.8 m.s-2",
            9.8 * sr.units.meter / (sr.units.second * sr.units.second),
        ),
        (
            "9.8 meters/second**2",
            9.8 * sr.units.meter / (sr.units.second * sr.units.second),
        ),
        (
            "5 m**2 / s**2",
            5 * sr.units.meter * sr.units.meter / (sr.units.second * sr.units.second),
        ),
        (
            "5 (m / s)**2",
            5 * sr.units.meter * sr.units.meter / (sr.units.second * sr.units.second),
        ),
        (
            "10 m ** 3 / (s ** 2 * kg)",
            10
            * sr.units.meter
            * sr.units.meter
            * sr.units.meter
            / (sr.units.second * sr.units.second * sr.units.kilogram),
        ),
        ("35.1°", 35.1 * sr.units.degree),
        ("90 degrees", 90 * sr.units.degrees),
        ("3.141 radians", 3.141 * sr.units.radians),
        ('12"', 12 * sr.units.inch),
        ("10 |e|", 10 * sr.units.mod_electron),
        ("25°C", U(25 * sr.units.celsius)),
        ("100 fahrenheit", U(100 * sr.units.fahrenheit)),
        ("0 celsius", U(0 * sr.units.celsius)),
        ("13 J K-1", 13 * sr.units.joule / sr.units.kelvin),
        ("60 miles per hour", 60 * sr.units.miles / sr.units.hour),
        ("70mph", 70 * sr.units.miles / sr.units.hour),
        ("120 kph", 120 * sr.units.kilometer / sr.units.hour),
        (
            "radian * degree**2 / radian^2",
            sr.units.degree * sr.units.degree / sr.units.radian,
        ),
        (
            "angstrom**3 / nanometer",
            sr.units.angstrom
            * sr.units.angstrom
            * sr.units.angstrom
            / (sr.units.nanometer),
        ),
        (
            "coulombs * angstrom**-2 * nanometer**2",
            sr.units.coulomb
            * sr.units.nanometer
            * sr.units.nanometer
            / (sr.units.angstrom * sr.units.angstrom),
        ),
        (
            "kcal_per_mol / angstrom**2 * nanometer**2",
            sr.units.kcal_per_mol
            / (
                sr.units.angstrom
                * sr.units.angstrom
                * sr.units.nanometer
                * sr.units.nanometer
            ),
        ),
        (
            "angstrom**3 * nanometer^-1 / picometer",
            sr.units.angstrom
            * sr.units.angstrom
            * sr.units.angstrom
            / (sr.units.nanometer * sr.units.picometer),
        ),
        (
            "bar * kJ_per_mol**2 / (kcal_per_mol * kJ_per_mol)",
            sr.units.bar * sr.units.kJ_per_mol / (sr.units.kcal_per_mol),
        ),
        (
            "coulomb * kelvin^-2 * kelvin^3 / e_charge",
            sr.units.coulomb * sr.units.kelvin / sr.units.e_charge,
        ),
        (
            "nanoseconds^3 * kelvin^-3 * kelvin**3 / milliseconds**2",
            sr.units.nanosecond
            * sr.units.nanosecond
            * sr.units.nanosecond
            / (sr.units.millisecond * sr.units.millisecond),
        ),
        (
            "angstroms3 * atm^-3 * bar**3",
            sr.units.angstrom
            * sr.units.angstrom
            * sr.units.angstrom
            * sr.units.bar
            * sr.units.bar
            * sr.units.bar
            / (sr.units.atm * sr.units.atm * sr.units.atm),
        ),
        ("degree", sr.units.degree),
        ("meters2", sr.units.meter * sr.units.meter),
        ("coulombs", sr.units.coulomb),
        ("kJ_per_mol", sr.units.kJ_per_mol),
        ("nanometer", sr.units.nanometer),
        ("bar", sr.units.bar),
        ("fahrenheit", U(sr.units.fahrenheit)),
        ("days", sr.units.day),
        (
            "picometers**3",
            sr.units.picometer * sr.units.picometer * sr.units.picometer,
        ),
        (
            "kcal per mol / angstrom**2",
            sr.units.kcal_per_mol / (sr.units.angstrom * sr.units.angstrom),
        ),
    ],
)
def test_parse_unit(text, expect):
    u = U(text)

    assert u.has_same_units(expect)
    assert u.value() == pytest.approx(expect.value())


@pytest.mark.parametrize(
    "initial, final, expect",
    [
        (
            "10 angstrom",
            "nanometer",
            (10 * sr.units.angstrom).to(sr.units.nanometer),
        ),
        (
            "15 kcal per mol",
            "kJ.mol-1",
            (15 * sr.units.kcal_per_mol).to(sr.units.kJ_per_mol),
        ),
        ("0 celsius", "kelvin", 273.15),
        ("273.15 K", "celsius", 0),
    ],
)
def test_convert_unit(initial, final, expect):
    result = U(initial).to(final)

    assert result == pytest.approx(expect)


@pytest.mark.parametrize(
    "unit, expect",
    [
        ("5 meter", U("5 meter")),
        ("0.5 kcal mol-1", 0.5 * sr.units.kcal_per_mol),
    ],
)
def test_create_unit(unit, expect):
    u = sr.u(unit)

    assert u.has_same_units(expect)
    assert u.value() == pytest.approx(expect.value())


@pytest.mark.parametrize(
    "unit1, unit2, tol",
    [
        ("m", "s", 1e-6),
        ("cal", "mol", 1e-6),
        ("J", "mol", 1e-3),
        ("g", "m", 1e-6),
        ("m", "g", 1e-6),
    ],
)
def test_default_units(unit1, unit2, tol):
    unit1_u = sr.u(unit1)
    unit2_u = sr.u(unit2)

    expect = sr.u(f"{unit1} {unit2}-1")

    old_unit1 = sr.u(unit1).get_default().unit_string()
    old_unit2 = sr.u(unit2).get_default().unit_string()

    sr.units.set_default_units([unit1, unit2])

    try:
        u = sr.u(f"{unit1} {unit2}-1")

        assert u.has_same_units(expect)
        assert u.value() == pytest.approx(expect.value(), tol)

        u2 = sr.u(str(u))
        assert u2.has_same_units(u)
        assert u2.value() == pytest.approx(u.value(), 1e-3)

        u = sr.u(f"m{unit1} {unit2}-1")

        u2 = sr.u(str(u))
        assert u2.has_same_units(u)
        assert u2.value() == pytest.approx(u.value())

        assert u.has_same_units(expect)
        assert u.value() == pytest.approx(1e-3 * expect.value(), tol)

        u = sr.u(f"{unit1} m{unit2}-1")

        u2 = sr.u(str(u))
        assert u2.has_same_units(u)
        assert u2.value() == pytest.approx(u.value(), tol)

        assert u.has_same_units(expect)
        assert u.value() == pytest.approx(1e3 * expect.value(), tol)

        u = sr.u(f"m{unit1} m{unit2}-1")

        u2 = sr.u(str(u))
        assert u2.has_same_units(u)
        assert u2.value() == pytest.approx(u.value(), tol)

        assert u.has_same_units(expect)
        assert u.value() == pytest.approx(expect.value(), tol)

        u = sr.u(f"{unit1} {unit2}-3")

        u2 = sr.u(str(u))
        assert u2.has_same_units(u)
        assert u2.value() == pytest.approx(u.value(), tol)

        assert u.has_same_units(expect / (unit2_u * unit2_u))
        assert u.value() == pytest.approx(expect.value() / (unit2_u.value() ** 2), tol)

        u = sr.u(f"{unit1}3 {unit2}-1")

        u2 = sr.u(str(u))
        assert u2.has_same_units(u)
        assert u2.value() == pytest.approx(u.value(), tol)

        assert u.has_same_units(expect * unit1_u * unit1_u)
        assert u.value() == pytest.approx(expect.value() * (unit1_u.value() ** 2), tol)

        u = sr.u(f"{unit1}3 {unit2}-3")

        u2 = sr.u(str(u))
        assert u2.has_same_units(u)
        assert u2.value() == pytest.approx(u.value(), tol)

        assert u.has_same_units(expect.pow(3))
        assert u.value() == pytest.approx(expect.pow(3).value(), tol)

    except Exception as e:
        sr.units.set_default_units([old_unit1, old_unit2])
        raise e

    sr.units.set_default_units([old_unit1, old_unit2])


def test_temperature():
    assert sr.u("1 kelvin").to("celsius") == pytest.approx(-272.15)
    assert sr.u("298.15 K").to("celsius") == pytest.approx(25)
    assert sr.u("25 celsius").to("kelvin") == pytest.approx(298.15)
    assert sr.u("0 celsius").to("kelvin") == pytest.approx(273.15)
    assert sr.u("0 kelvin").to("celsius") == pytest.approx(-273.15)
