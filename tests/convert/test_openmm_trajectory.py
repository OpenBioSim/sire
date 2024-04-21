import sire as sr
import pytest


@pytest.mark.veryslow
@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_openmm_trajectory(ala_mols, openmm_platform, tmpdir):
    mols = ala_mols

    mols = mols.minimisation(platform=openmm_platform).run().commit()

    d = mols.dynamics(timestep="1fs", temperature="25oC", platform=openmm_platform)

    d.run("0.1ps", save_frequency="0.01ps")

    mols = d.commit()

    assert mols.num_frames() == 10

    rmsd = mols.trajectory().rmsd()

    assert len(rmsd) == 10

    dir = tmpdir.mkdir("test_openmm_trajectory")

    f = sr.save(mols.trajectory(), dir.join("test"), format=["PRMTOP", "RST"])

    mols2 = sr.load(f)

    rmsd2 = mols2.trajectory().rmsd()

    assert len(rmsd2) == 10

    for r1, r2 in zip(rmsd, rmsd2):
        assert r1.value() == pytest.approx(r2.value(), abs=1e-3)
