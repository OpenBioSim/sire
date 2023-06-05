import pytest
import sire as sr

U = sr.units.GeneralUnit


@pytest.mark.parametrize(
    "text, expect",
    [
        ("5", U(5)),
        ("-3.2", U(-3.2)),
        ("-3.141e-10", U(-3.141e-10)),
        ("1 angstrom", 1 * sr.units.angstrom),
        ("1 angstrom /* block comment */", 1 * sr.units.angstrom),
        ("1 angstrom // line comment", 1 * sr.units.angstrom),
        ("-10.5e5 A", -10.5e5 * sr.units.angstrom),
        ("10 s-1", 10 / sr.units.second),
        ("6.5 m s-1", 6.5 * sr.units.meter / sr.units.second),
        ("1.2 m.s-1", 1.2 * sr.units.meter / sr.units.second),
        ("55 mm", 55 * sr.units.millimeter),
        ("15 ps", 15 * sr.units.picosecond),
        ("20 angstroms", 20 * sr.units.angstrom),
        ("20 angstrom", 20 * sr.units.angstrom),
        ("15 kcal.mol-1", 15 * sr.units.kcal_per_mol),
        (
            "1 kJ/(mol nm**2)",
            1
            * sr.units.kJ_per_mol
            / (sr.units.nanometer * sr.units.nanometer),
        ),
        ("-34.5 kilocalorie mole^-1", -34.5 * sr.units.kcal_per_mol),
        (
            "1.2 kcal mol-1 A-2",
            1.2
            * sr.units.kcal_per_mol
            / (sr.units.angstrom * sr.units.angstrom),
        ),
        (
            "8 (kcal per mol) / angstrom**2",
            8
            * sr.units.kcal_per_mol
            / (sr.units.angstrom * sr.units.angstrom),
        ),
        (
            "8 kcal per mol / angstrom**2",
            8
            * sr.units.kcal_per_mol
            / (sr.units.angstrom * sr.units.angstrom),
        ),
        ("1 fs", 1 * sr.units.femtosecond),
        ("1e-3 femtoseconds", 1e-3 * sr.units.femtosecond),
        (
            "5 (kJ per mol) / nm**2",
            5
            * sr.units.kJ_per_mol
            / (sr.units.nanometer * sr.units.nanometer),
        ),
        (
            "5 kJ per mol / nm**2",
            5
            * sr.units.kJ_per_mol
            / (sr.units.nanometer * sr.units.nanometer),
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
            5
            * sr.units.meter
            * sr.units.meter
            / (sr.units.second * sr.units.second),
        ),
        (
            "5 (m / s)**2",
            5
            * sr.units.meter
            * sr.units.meter
            / (sr.units.second * sr.units.second),
        ),
        (
            "10 m ** 3 / (s ** 2 * kg)",
            10
            * sr.units.meter
            * sr.units.meter
            * sr.units.meter
            / (sr.units.second * sr.units.second * sr.units.kilogram),
        ),
    ],
)
def test_parse_unit(text, expect):
    u = U(text)

    assert u.has_same_units(expect)
    assert u.value() == pytest.approx(expect.value())
