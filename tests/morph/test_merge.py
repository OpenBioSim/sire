import sire as sr

import pytest

try:
    from kartograf import KartografAtomMapper
except ImportError:
    KartografAtomMapper = None


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
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

    nrg_remerged_0 = remerged.dynamics(
        lambda_value=0, platform=openmm_platform
    ).current_potential_energy()

    nrg_remerged_1 = remerged.dynamics(
        lambda_value=1, platform=openmm_platform
    ).current_potential_energy()

    assert nrg_merged_0.value() == pytest.approx(nrg_remerged_0.value())
    assert nrg_merged_1.value() == pytest.approx(nrg_remerged_1.value())


@pytest.mark.slow
@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats() or KartografAtomMapper is None,
    reason="openmm support is not available",
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

    merged = merged.perturbation().link_to_reference()

    nrg_merged_0 = merged.dynamics(
        lambda_value=0, platform=openmm_platform
    ).current_potential_energy()

    nrg_merged_1 = merged.dynamics(
        lambda_value=1, platform=openmm_platform
    ).current_potential_energy()

    print(ose_nrg, zan_nrg)
    print(extracted_ose_nrg, extracted_zan_nrg)
    print(nrg_merged_0, nrg_merged_1)

    merged = merged.perturbation().link_to_reference()

    nrg_merged_0 = merged.dynamics(
        lambda_value=0, platform=openmm_platform
    ).current_potential_energy()

    nrg_merged_1 = merged.dynamics(
        lambda_value=1, platform=openmm_platform
    ).current_potential_energy()

    print(nrg_merged_0, nrg_merged_1)

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

    merged = merged.perturbation().link_to_reference()

    nrg_merged_0 = merged.dynamics(
        lambda_value=0, platform=openmm_platform
    ).current_potential_energy()

    nrg_merged_1 = merged.dynamics(
        lambda_value=1, platform=openmm_platform
    ).current_potential_energy()

    print(zan_nrg, ose_nrg)
    print(extracted_zan_nrg, extracted_ose_nrg)
    print(nrg_merged_0, nrg_merged_1)

    merged = merged.perturbation().link_to_reference()

    nrg_merged_0 = merged.dynamics(
        lambda_value=0, platform=openmm_platform
    ).current_potential_energy()

    nrg_merged_1 = merged.dynamics(
        lambda_value=1, platform=openmm_platform
    ).current_potential_energy()

    print(nrg_merged_0, nrg_merged_1)

    assert ose_nrg.value() == pytest.approx(nrg_merged_0.value())
    assert zan_nrg.value() == pytest.approx(nrg_merged_1.value())


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
