import sire as sr
import pytest


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_openmm_single_energy(kigaki_mols):
    mols = kigaki_mols

    mol = mols[0]

    map = {
        "space": sr.vol.Cartesian(),
        "platform": "Reference",
        "constraint": "bonds-h-angles",
    }

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
    assert mol.energy(map=map).to(sr.units.kJ_per_mol) == pytest.approx(energy, abs=0.1)


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_openmm_multi_energy_small_cart(kigaki_mols):
    # first, try just 50 molecules in a cartesian space
    mols = kigaki_mols[0:50]

    map = {
        "space": sr.vol.Cartesian(),
        "platform": "Reference",
        "constraint": "bonds-h-angles",
    }

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
@pytest.mark.slow
def test_openmm_multi_energy_all_cart(kigaki_mols):
    # use all of the molecules
    mols = kigaki_mols

    map = {
        "space": sr.vol.Cartesian(),
        "cutoff": 10000 * sr.units.angstrom,
        "cutoff_type": "REACTION_FIELD",
        "dielectric": 1.0,
        "platform": "cpu",
        "constraint": "bonds-h-angles",
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
        "platform": "cpu",
        "constraint": "bonds-h-angles",
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

    map = {
        "cutoff": 10 * sr.units.angstrom,
        "cutoff_type": "REACTION_FIELD",
        "dielectric": 78.0,
        "platform": "Reference",
        "constraint": "bonds-h-angles",
    }

    omm = sr.convert.to(mols, "openmm", map=map)

    state = omm.getState(getEnergy=True, enforcePeriodicBox=True)

    energy = state.getPotentialEnergy()

    # get this as a float in kJ mol-1
    energy = energy.value_in_unit(energy.unit)

    # -74975.9 +/- 0.5
    assert mols.energy(map=map).to(sr.units.kJ_per_mol) == pytest.approx(
        energy, abs=0.5
    )


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
@pytest.mark.slow
def test_openmm_dynamics(ala_mols):
    mols = ala_mols

    map = {
        "cutoff": 10 * sr.units.angstrom,
        "cutoff_type": "REACTION_FIELD",
        "dielectric": 78.0,
        "temperature": 25 * sr.units.celsius,
        "platform": "Reference",
        "constraint": "bonds-h-angles",
        # "pressure": 1 * sr.units.atm,   # currently disagree with energies for NPT...
    }

    sire_nrg = mols.energy(map=map)

    # Need to set constraints to `none` so that we
    # get energy agreement with sire - without this
    # we will be missing the bond energies
    # (and some angle energies)

    d = mols.dynamics(
        timestep=1 * sr.units.femtosecond,
        save_frequency=1 * sr.units.picosecond,
        map=map,
        constraint="none",
    )

    omm_nrg = d.current_potential_energy()

    assert sire_nrg.value() == pytest.approx(omm_nrg.value(), abs=0.5)

    assert d.ensemble().is_canonical()
    assert d.ensemble().temperature() == 25 * sr.units.celsius

    assert d.timestep() == 1 * sr.units.femtosecond

    d.run(0.1 * sr.units.picosecond, 0.01 * sr.units.picosecond)

    assert d.current_step() == 100
    # the molecules started from 6 ns
    assert d.current_time().to(sr.units.picosecond) == pytest.approx(6000.1)

    mols = d.commit()

    assert mols.num_frames() == 10

    sire_nrg = mols.energy(map=map)

    omm_nrg = d.current_potential_energy()

    assert sire_nrg.value() == pytest.approx(omm_nrg.value(), abs=0.5)


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_openmm_options(ala_mols):
    mols = ala_mols

    mol = mols[0]

    m = {
        "integrator": "langevin_middle",
        "temperature": 25 * sr.units.celsius,
        "pressure": 1 * sr.units.atm,
        "friction": 5 / sr.units.picosecond,
        "platform": "Reference",
        "constraint": "bonds-h-angles",
    }

    omm = sr.convert.to(mol, "openmm", map=m)

    for platform in ["CPU", "OpenCL", "CUDA"]:
        m["platform"] = platform

        try:
            omm = sr.convert.to(mol, "openmm", map=m)
        except ValueError:
            # maybe OpenCL or CUDA are not supported
            pass


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_openmm_ignore_constrained(ala_mols):
    mols = ala_mols

    mol = mols[0]

    d = mol.dynamics(
        constraint="bonds-h-angles",
        include_constrained_energies=True,
        platform="Reference",
    )

    nrg1 = d.current_potential_energy()

    d = mol.dynamics(
        constraint="bonds-h-angles",
        include_constrained_energies=False,
        platform="Reference",
    )

    nrg2 = d.current_potential_energy()

    # these two energies should be different, because
    # we should be ignoring the constrained bonds and angles
    assert abs(nrg2.value() - nrg1.value()) > 1.0


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_openmm_no_zero_sigmas(zero_lj_mols):
    mols = zero_lj_mols

    omm = sr.convert.to(mols, "openmm", 
                        map={"constraint": "h-bonds",
                             "platform": "Reference"})

    from openmm import XmlSerializer

    xml = XmlSerializer.serialize(omm.getSystem())

    for line in xml.split("\n"):
        assert 'sig="0"' not in line


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_openmm_skipped_constrained_bonds(zero_lj_mols):
    mols = zero_lj_mols

    omm1 = sr.convert.to(
        mols,
        "openmm",
        map={"constraint": "h-bonds", 
             "include_constrained_energies": True,
             "platform": "Reference"},
    )

    omm2 = sr.convert.to(
        mols,
        "openmm",
        map={"constraint": "h-bonds", 
             "include_constrained_energies": False,
             "platform": "Reference"},
    )

    nrg1 = omm1.get_potential_energy().to(sr.units.kcal_per_mol)
    nrg2 = omm2.get_potential_energy().to(sr.units.kcal_per_mol)

    # Check the energies haven't changed
    # (regression check - here are the current values)
    assert nrg1 == pytest.approx(-447.44, 1e-3)
    assert nrg2 == pytest.approx(-3279.87, 1e-3)

    from openmm import XmlSerializer

    xml1 = XmlSerializer.serialize(omm1.getSystem())
    xml2 = XmlSerializer.serialize(omm2.getSystem())

    assert xml1 != xml2

    lines1 = xml1.split("\n")
    lines2 = xml2.split("\n")

    i = 0

    for j in range(0, len(lines2)):
        line1 = lines1[i]
        line2 = lines2[j]
        i += 1

        while line1 != line2:
            line1 = lines1[i]
            assert "Bond" in line1
            i += 1
