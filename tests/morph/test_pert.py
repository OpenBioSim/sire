import sire as sr

import pytest

# get the directory of this script file
import os

test_dir = os.path.dirname(os.path.realpath(__file__))

neopentane_methane_pert = os.path.join(test_dir, "input", "neopentane_methane.pert")


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_pertfile(neopentane_methane, openmm_platform):
    mols = neopentane_methane.clone()

    map = {
        "space": sr.vol.Cartesian(),
        "platform": openmm_platform,
        "constraint": "none",
        "cutoff": "none",
    }

    mols = sr.morph.link_to_reference(mols)

    # (the pert file has zeroed ghost torsions)
    mols = sr.morph.zero_ghost_torsions(mols)

    mols2 = sr.morph.extract_reference(mols)

    with pytest.raises(KeyError):
        assert len(mols2.molecules("property is_perturbable")) == 0

    mols2.update(sr.morph.create_from_pertfile(mols2[0], neopentane_methane_pert))

    assert len(mols2.molecules("property is_perturbable")) == 1

    p1 = mols[0].perturbation().to_openmm(constraint="bonds")
    p2 = mols2[0].perturbation().to_openmm(constraint="bonds")

    c_atm1 = p1.changed_atoms(to_pandas=False)
    c_atm2 = p2.changed_atoms(to_pandas=False)

    assert len(c_atm1) == len(c_atm2)

    for atm1 in c_atm1:
        found = False
        for atm2 in c_atm2:
            found = True
            if atm1[0].name() == atm2[0].name():
                for i in range(1, len(atm1)):
                    assert atm1[i] == pytest.approx(atm2[i], abs=1e-3)

                break

        assert found

    c_bnd1 = p1.changed_bonds(to_pandas=False)
    c_bnd2 = p2.changed_bonds(to_pandas=False)

    assert len(c_bnd1) == len(c_bnd2)

    for bnd1 in c_bnd1:
        found = False
        for bnd2 in c_bnd2:
            found = True
            if bnd1[0].id() == bnd2[0].id():
                for i in range(1, len(bnd1)):
                    assert bnd1[i] == pytest.approx(bnd2[i], abs=1e-3)

                break

        assert found

    c_ang1 = p1.changed_angles(to_pandas=False)
    c_ang2 = p2.changed_angles(to_pandas=False)

    assert len(c_ang1) == len(c_ang2)

    for ang1 in c_ang1:
        found = False
        for ang2 in c_ang2:
            found = True
            if ang1[0].id() == ang2[0].id():
                for i in range(1, len(ang1)):
                    assert ang1[i] == pytest.approx(ang2[i], abs=1e-3)

                break

        assert found

    c_tor1 = p1.changed_torsions(to_pandas=False)
    c_tor2 = p2.changed_torsions(to_pandas=False)

    assert len(c_tor1) == len(c_tor2)

    for tor1 in c_tor1:
        found = False
        for tor2 in c_tor2:
            found = True
            if tor1[0].id() == tor2[0].id():
                for i in range(1, len(tor1)):
                    assert tor1[i] == pytest.approx(tor2[i], abs=1e-3)

                break

        assert found

    c_con1 = p1.changed_constraints(to_pandas=False)
    c_con2 = p2.changed_constraints(to_pandas=False)

    assert len(c_con1) == len(c_con2)

    for con1 in c_con1:
        found = False
        for con2 in c_con2:
            found = True
            if con1[0][0].name() == con2[0][0].name():
                for i in range(1, len(con1)):
                    assert con1[i] == pytest.approx(con2[i], abs=1e-3)

                break

        assert found

    d = mols.dynamics(lambda_value=0.0, map=map)
    d2 = mols2.dynamics(lambda_value=0.0, map=map)

    assert d.current_potential_energy().value() == pytest.approx(
        d2.current_potential_energy().value(), abs=1e-3
    )

    d = mols.dynamics(lambda_value=1.0, map=map)
    d2 = mols2.dynamics(lambda_value=1.0, map=map)

    assert d.current_potential_energy().value() == pytest.approx(
        d2.current_potential_energy().value(), abs=1e-3
    )

    d = mols.dynamics(lambda_value=0.5, map=map)
    d2 = mols2.dynamics(lambda_value=0.5, map=map)

    assert d.current_potential_energy().value() == pytest.approx(
        d2.current_potential_energy().value(), abs=1e-3
    )


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_extract_and_link(neopentane_methane, openmm_platform):
    mols = neopentane_methane.clone()

    map = {
        "space": sr.vol.Cartesian(),
        "platform": openmm_platform,
        "constraint": "none",
        "cutoff": "none",
    }

    assert len(mols.molecules("property is_perturbable")) == 1

    ref_mols = sr.morph.extract_reference(mols, remove_ghosts=False)
    pert_mols = sr.morph.extract_perturbed(mols, remove_ghosts=False)

    with pytest.raises(KeyError):
        assert len(ref_mols.molecules("property is_perturbable")) == 0

    mols = sr.morph.link_to_reference(mols)

    nrg_0_0 = (
        mols.dynamics(ignore_perturbations=True, map=map)
        .current_potential_energy()
        .value()
    )

    nrg_1_0 = (
        mols.dynamics(swap_end_states=True, ignore_perturbations=True, map=map)
        .current_potential_energy()
        .value()
    )

    mols = sr.morph.link_to_perturbed(mols)

    nrg_0_1 = (
        mols.dynamics(ignore_perturbations=True, map=map)
        .current_potential_energy()
        .value()
    )

    nrg_1_1 = (
        mols.dynamics(swap_end_states=True, ignore_perturbations=True, map=map)
        .current_potential_energy()
        .value()
    )

    assert nrg_0_0 != nrg_1_0
    assert nrg_0_1 != nrg_1_1

    assert ref_mols.energy().value() == pytest.approx(nrg_0_0, abs=1e-3)
    assert pert_mols.energy().value() == pytest.approx(nrg_1_1, abs=1e-3)

    mols = sr.morph.link_to_reference(mols)

    assert nrg_0_0 == pytest.approx(
        mols.dynamics(lambda_value=0.0, map=map).current_potential_energy().value(),
        1e-3,
    )

    assert nrg_1_0 == pytest.approx(
        mols.dynamics(lambda_value=1.0, map=map).current_potential_energy().value(),
        1e-3,
    )

    for key in ref_mols[0].property_keys():
        assert mols[0].property(key) == ref_mols[0].property(key)

    nrg_ref = ref_mols.dynamics(map=map).current_potential_energy().value()

    assert nrg_0_0 == pytest.approx(nrg_ref, 1e-3)

    with pytest.raises(KeyError):
        assert len(pert_mols.molecules("property is_perturbable")) == 0

    mols = sr.morph.link_to_perturbed(mols)

    for key in ref_mols[0].property_keys():
        assert mols[0].property(key) == pert_mols[0].property(key)

    nrg_pert = pert_mols.dynamics(map=map).current_potential_energy().value()

    assert nrg_1_1 == pytest.approx(nrg_pert, 1e-3)
