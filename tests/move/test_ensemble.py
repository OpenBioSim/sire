import sire as sr


def test_ensemble():
    from sire.move import Ensemble

    e = Ensemble()

    assert e.is_micro_canonical()

    e = Ensemble(map={"temperature": 25 * sr.units.celsius})

    assert e.is_canonical()

    assert e.temperature().value() == (25 * sr.units.celsius).value()

    e = Ensemble(
        map={"temperature": 53 * sr.units.kelvin, "pressure": 2 * sr.units.atm}
    )

    assert e.is_isothermal_isobaric()

    assert e.temperature().value() == (53 * sr.units.kelvin).value()
    assert e.pressure().value() == (2 * sr.units.atm).value()

    e = Ensemble("nve")

    assert e.is_nve()

    e = Ensemble("nvt")

    assert e.is_nvt()
    assert e.temperature().value() == (25 * sr.units.celsius).value()

    e = Ensemble("npt")

    assert e.is_npt()
    assert e.temperature().value() == (25 * sr.units.celsius).value()
    assert e.pressure().value() == (1 * sr.units.atm).value()

    assert Ensemble("nve") == Ensemble(map={"ensemble": "nve"})
    assert Ensemble("nvt") == Ensemble(map={"ensemble": "nvt"})
    assert Ensemble("npt") == Ensemble(map={"ensemble": "npt"})
