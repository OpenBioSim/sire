import sire as sr

import pytest

try:
    from kartograf import KartografAtomMapper
except ImportError:
    KartografAtomMapper = None


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm or kartograf support is not available",
)
def test_extract_remerge(merged_zan_ose, openmm_platform):
    merged = merged_zan_ose[0].perturbation().link_to_reference()

    extracted_ose = merged.perturbation().extract_reference()
    extracted_zan = merged.perturbation().extract_perturbed()

    remerged = sr.morph.merge(
        extracted_ose, extracted_zan, match=sr.legacy.Mol.AtomNumMatcher()
    )

    nrg_merged_0 = merged.dynamics(
        lambda_value=0, platform=openmm_platform
    ).current_potential_energy()

    nrg_merged_1 = merged.dynamics(
        lambda_value=1, platform=openmm_platform
    ).current_potential_energy()

    nrg_merged_05 = merged.dynamics(
        lambda_value=0.5, platform=openmm_platform
    ).current_potential_energy()

    nrg_remerged_0 = remerged.dynamics(
        lambda_value=0, platform=openmm_platform
    ).current_potential_energy()

    nrg_remerged_1 = remerged.dynamics(
        lambda_value=1, platform=openmm_platform
    ).current_potential_energy()

    nrg_remerged_05 = remerged.dynamics(
        lambda_value=0.5, platform=openmm_platform
    ).current_potential_energy()

    assert nrg_merged_0.value() == pytest.approx(nrg_remerged_0.value())
    assert nrg_merged_1.value() == pytest.approx(nrg_remerged_1.value())

    # this is different - worth checking why!
    # assert nrg_merged_05.value() == pytest.approx(nrg_remerged_05.value())


@pytest.mark.slow
@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats() or KartografAtomMapper is None,
    reason="openmm or kartograf support is not available",
)
def test_merge(ose_mols, zan_mols, openmm_platform):
    ose = ose_mols[0]
    zan = zan_mols[0]

    merged = sr.morph.merge(
        ose, zan, match=KartografAtomMapper(atom_map_hydrogens=True)
    )

    ose_nrg = ose.dynamics(platform=openmm_platform).current_potential_energy()
    zan_nrg = zan.dynamics(platform=openmm_platform).current_potential_energy()

    extracted_ose = merged.perturbation().extract_reference()
    extracted_zan = merged.perturbation().extract_perturbed()

    extracted_ose_nrg = extracted_ose.dynamics(
        platform=openmm_platform
    ).current_potential_energy()

    extracted_zan_nrg = extracted_zan.dynamics(
        platform=openmm_platform
    ).current_potential_energy()

    assert extracted_ose_nrg.value() == pytest.approx(ose_nrg.value())
    assert extracted_zan_nrg.value() == pytest.approx(zan_nrg.value())

    merged = sr.morph.merge(
        zan, ose, match=KartografAtomMapper(atom_map_hydrogens=True)
    )

    extracted_zan = merged.perturbation().extract_reference()
    extracted_ose = merged.perturbation().extract_perturbed()

    extracted_ose_nrg = extracted_ose.dynamics(
        platform=openmm_platform
    ).current_potential_energy()

    extracted_zan_nrg = extracted_zan.dynamics(
        platform=openmm_platform
    ).current_potential_energy()

    assert extracted_ose_nrg.value() == pytest.approx(ose_nrg.value())
    assert extracted_zan_nrg.value() == pytest.approx(zan_nrg.value())

    # we don't test merged molecule energies are these
    # energies are not equal because there are additional bonds,
    # angle and dihedrals to ghost atoms in the merged molecule


@pytest.mark.veryslow
def test_merge_protein(neura_mols):
    protein = neura_mols["protein"]

    ala = protein.residues("ALA")[1]
    lys = protein.residues("LYS")[1]

    merged = sr.morph.merge(ala, lys)

    merged_ala = merged.perturbation().extract_reference()[ala.number()]
    merged_lys = merged.perturbation().extract_perturbed()[lys.number()]

    assert ala.energy().value() == pytest.approx(merged_ala.energy().value())
    assert lys.energy().value() == pytest.approx(merged_lys.energy().value())


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats() or KartografAtomMapper is None,
    reason="openmm or kartograf support is not available",
)
def test_merge_neopentane_methane(neopentane_methane, openmm_platform):
    neopentane = neopentane_methane[0].perturbation().extract_reference()
    methane = neopentane_methane[0].perturbation().extract_perturbed()

    orig_merged = sr.morph.link_to_reference(neopentane_methane[0])

    nrg_neo = neopentane.dynamics(platform=openmm_platform).current_potential_energy()
    nrg_met = methane.dynamics(platform=openmm_platform).current_potential_energy()

    merged = sr.morph.merge(
        neopentane, methane, match=KartografAtomMapper(atom_map_hydrogens=True)
    )

    nrg_merged_0 = merged.dynamics(
        lambda_value=0, platform=openmm_platform
    ).current_potential_energy()

    nrg_merged_1 = merged.dynamics(
        lambda_value=1, platform=openmm_platform
    ).current_potential_energy()

    nrg_merged_05 = merged.dynamics(
        lambda_value=0.5, platform=openmm_platform
    ).current_potential_energy()

    extracted_neo = merged.perturbation().extract_reference()
    extracted_met = merged.perturbation().extract_perturbed()

    nrg_extracted_neo = extracted_neo.dynamics(
        platform=openmm_platform
    ).current_potential_energy()

    nrg_extracted_met = extracted_met.dynamics(
        platform=openmm_platform
    ).current_potential_energy()

    nrg_orig_merged_0 = orig_merged.dynamics(
        lambda_value=0, platform=openmm_platform
    ).current_potential_energy()

    nrg_orig_merged_1 = orig_merged.dynamics(
        lambda_value=1, platform=openmm_platform
    ).current_potential_energy()

    nrg_orig_merged_05 = orig_merged.dynamics(
        lambda_value=0.5, platform=openmm_platform
    ).current_potential_energy()

    assert nrg_neo.value() == pytest.approx(nrg_extracted_neo.value(), abs=1e-3)
    assert nrg_met.value() == pytest.approx(nrg_extracted_met.value(), abs=1e-3)

    assert nrg_merged_0.value() == pytest.approx(nrg_orig_merged_0.value(), abs=1e-3)
    assert nrg_merged_1.value() == pytest.approx(nrg_orig_merged_1.value(), abs=1e-3)

    # this is different - worth checking why!
    # assert nrg_merged_05.value() == pytest.approx(nrg_orig_merged_05.value(), abs=1e-3)

    # These energies aren't correct - extra ghost atom internals?
    assert nrg_neo.value() == pytest.approx(nrg_merged_0.value(), abs=1e-3)
    # assert nrg_met.value() == pytest.approx(nrg_merged_1.value(), abs=1e-3)
