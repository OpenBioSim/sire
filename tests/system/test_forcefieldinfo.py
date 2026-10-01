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


def test_forcefieldinfo_string_options(kigaki_mols):
    m = {
        "cutoff": 5 * sr.units.angstrom,
        "cutoff_type": "PME",
        "tolerance": "5e-4",
        "pme_alpha": "3.47",
        "pme_grid": "64",
    }

    ffinfo = sr.system.ForceFieldInfo(kigaki_mols, m)

    assert ffinfo.get_parameter("tolerance") == pytest.approx(5e-4)
    assert ffinfo.get_parameter("pme_alpha") == pytest.approx(3.47)
    assert ffinfo.get_parameter("pme_grid_x") == 64

    m = {"cutoff": 5 * sr.units.angstrom, "cutoff_type": "RF", "dielectric": "80"}

    ffinfo = sr.system.ForceFieldInfo(kigaki_mols, m)

    assert ffinfo.get_parameter("dielectric") == pytest.approx(80)


@pytest.mark.parametrize(
    "options",
    [
        {"tolerance": "five"},
        {"pme_alpha": "big", "pme_grid": 64},
        {"pme_spacing": "0.12 nm"},
    ],
)
def test_forcefieldinfo_bad_string_options(kigaki_mols, options):
    m = {"cutoff": 5 * sr.units.angstrom, "cutoff_type": "PME"}
    m.update(options)

    with pytest.raises(RuntimeError, match="must be a"):
        sr.system.ForceFieldInfo(kigaki_mols, m)
