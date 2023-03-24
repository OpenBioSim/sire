import sire as sr
import pytest


def test_forcefieldinfo(kigaki_mols):
    mols = kigaki_mols

    ffinfo = sr.system.ForceFieldInfo(mols)

    assert ffinfo.space() == mols.property("space")

    assert ffinfo.has_cutoff()

    assert ffinfo.cutoff_type() == "CUTOFF"

    assert ffinfo.cutoff() > 0
    assert ffinfo.cutoff() < mols.property("space").maximum_cutoff()
    assert (
        mols.property("space").maximum_cutoff() - ffinfo.cutoff()
        < 1.5 * sr.units.angstrom
    )

    ffinfo = sr.system.ForceFieldInfo(mols, {"cutoff": 5 * sr.units.angstrom})

    assert ffinfo.space() == mols.property("space")

    assert ffinfo.has_cutoff()

    assert ffinfo.cutoff() == 5 * sr.units.angstrom

    with pytest.raises(ValueError):
        ffinfo = sr.system.ForceFieldInfo(
            mols, {"cutoff": 2 * mols.property("space").maximum_cutoff()}
        )

    ffinfo = sr.system.ForceFieldInfo(
        mols,
        {
            "cutoff": 5 * sr.units.angstrom,
            "cutoff_type": "PME",
            "tolerance": 0.5,
        },
    )

    assert ffinfo.space() == mols.property("space")

    assert ffinfo.has_cutoff()
    assert ffinfo.cutoff() == 5 * sr.units.angstrom

    assert ffinfo.cutoff_type() == "PME"
    assert ffinfo.get_parameter("tolerance") == 0.5
