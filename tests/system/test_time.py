def test_set_zero_time(ala_mols):
    import sire as sr

    mols = ala_mols.clone()

    mols.set_time(0 * sr.units.picosecond)
    mols.set_time(mols.time())
    assert mols.time().value() == 0

    mols.set_time(sr.u("5 ps"))
    mols.set_time(mols.time())
    assert mols.time().to(sr.units.picosecond) == 5
