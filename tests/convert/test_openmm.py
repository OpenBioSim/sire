import sire as sr
import pytest


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_openmm_single_energy(kigaki_mols):
    mols = kigaki_mols

    mol = mols[0]

    map = {"space": sr.vol.Cartesian()}

    omm = sr.convert.to(mol, "openmm", map=map)

    state = omm.getState(getPositions=True, getEnergy=True)

    positions = state.getPositions()

    assert len(positions) == mol.num_atoms()

    for i, atom in enumerate(mol.atoms()):
        c = atom.coords()
        p = positions[i]

        assert c.x().to(sr.units.nanometer) == pytest.approx(p.x)
        assert c.y().to(sr.units.nanometer) == pytest.approx(p.y)
        assert c.z().to(sr.units.nanometer) == pytest.approx(p.z)

    energy = state.getPotentialEnergy()

    # get this as a float in kJ mol-1
    energy = energy.value_in_unit(energy.unit)

    # these won't be exactly the same - this is 5227 +/- 0.1 kJ mol-1
    assert mol.energy(map=map).to(sr.units.kJ_per_mol) == pytest.approx(
        energy, abs=0.1
    )


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_openmm_multi_energy_small_cart(kigaki_mols):
    # first, try just 50 molecules in a cartesian space
    mols = kigaki_mols[0:50]

    map = {"space": sr.vol.Cartesian()}

    omm = sr.convert.to(mols, "openmm", map=map)

    state = omm.getState(getEnergy=True)

    energy = state.getPotentialEnergy()

    # get this as a float in kJ mol-1
    energy = energy.value_in_unit(energy.unit)

    # these won't be exactly the same - this is 4865.82 +/- 0.04
    assert mols.energy(map=map).to(sr.units.kJ_per_mol) == pytest.approx(
        energy, abs=0.5
    )


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_openmm_multi_energy_all_cart(kigaki_mols):
    # use all of the molecules
    mols = kigaki_mols

    map = {
        "space": sr.vol.Cartesian(),
        "cutoff": 10000 * sr.units.angstrom,
        "cutoff_type": "REACTION_FIELD",
        "dielectric": 1.0,
    }

    omm = sr.convert.to(mols, "openmm", map=map)

    state = omm.getState(getEnergy=True)

    energy = state.getPotentialEnergy()

    # get this as a float in kJ mol-1
    energy = energy.value_in_unit(energy.unit)

    # -127882 +/- 0.5
    assert mols.energy(map=map).to(sr.units.kJ_per_mol) == pytest.approx(
        energy, abs=1.0
    )


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_openmm_multi_energy_all_cart_cutoff(kigaki_mols):
    # use all of the molecules
    mols = kigaki_mols

    map = {
        "space": sr.vol.Cartesian(),
        "cutoff": 10 * sr.units.angstrom,
        "cutoff_type": "REACTION_FIELD",
        "dielectric": 78.0,
    }

    omm = sr.convert.to(mols, "openmm", map=map)

    state = omm.getState(getEnergy=True)

    energy = state.getPotentialEnergy()

    # get this as a float in kJ mol-1
    energy = energy.value_in_unit(energy.unit)

    # -125622.2 +/- 0.05
    assert mols.energy(map=map).to(sr.units.kJ_per_mol) == pytest.approx(
        energy, abs=0.5
    )


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_openmm_multi_energy_all_periodic_cutoff(kigaki_mols):
    # use all of the molecules
    mols = kigaki_mols

    # GET DISAGREEMNT FROM FIRST MOLECULE, LIKELY BECAUSE OF NO
    # SPACE IN THE FIRST MOLECULE!

    map = {
        "space": mols.property("space"),
        "cutoff": 10 * sr.units.angstrom,
        "cutoff_type": "REACTION_FIELD",
        "dielectric": 78.0,
    }

    omm = sr.convert.to(mols, "openmm", map=map)

    state = omm.getState(getEnergy=True)

    energy = state.getPotentialEnergy()

    # get this as a float in kJ mol-1
    energy = energy.value_in_unit(energy.unit)

    # -74975.9 +/- 0.5
    assert mols.energy(map=map).to(sr.units.kJ_per_mol) == pytest.approx(
        energy, abs=0.5
    )
