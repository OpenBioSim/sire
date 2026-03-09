import sire as sr

import pytest


def test_glycam(tmpdir):
    """Test that a topology using the GLYCAM force field (SCEE=1.0, SCNB=1.0)
    round-trips correctly through AMBER prm7 format.

    GLYCAM uses full 1-4 interactions (no scaling), so SCEE=1.0 and SCNB=1.0
    for glycan dihedrals. The protein dihedrals use standard AMBER scaling
    (SCEE=1.2, SCNB=2.0). Before the fix, CLJScaleFactor(1.0, 1.0) pairs were
    silently dropped when building AmberParams, so glycan dihedrals were written
    with SCEE=0 and SCNB=0, giving zero 1-4 interactions.
    """

    # Load the GLYCAM topology and coordinates.
    mols = sr.load_test_files("glycam.top", "glycam.gro")

    # Write to AMBER prm7 + rst7 format.
    d = tmpdir.mkdir("test_glycam_amber")
    f = sr.save(mols, d.join("glycam_out"), format=["PRM7", "RST7"])

    # Parse SCEE_SCALE_FACTOR and SCNB_SCALE_FACTOR from the written prm7.
    # Format: 5 values per line, 16 chars each (AmberFormat FLOAT 5 16 8).
    scee_values = []
    scnb_values = []
    reading = None

    with open(f[0], "r") as fh:
        for line in fh:
            if line.startswith("%FLAG SCEE_SCALE_FACTOR"):
                reading = "scee"
                continue
            elif line.startswith("%FLAG SCNB_SCALE_FACTOR"):
                reading = "scnb"
                continue
            elif line.startswith("%FLAG") or line.startswith("%FORMAT"):
                if line.startswith("%FLAG"):
                    reading = None
                continue
            if reading == "scee":
                scee_values.extend(float(x) for x in line.split())
            elif reading == "scnb":
                scnb_values.extend(float(x) for x in line.split())

    assert len(scee_values) > 0, "No SCEE_SCALE_FACTOR entries found"
    assert len(scnb_values) > 0, "No SCNB_SCALE_FACTOR entries found"

    # The system has both glycan (SCEE=1.0, SCNB=1.0) and protein
    # (SCEE=1.2, SCNB=2.0) dihedrals. Both values must be present.
    # Before the fix, all glycan entries would be 0.0.
    assert any(
        v == pytest.approx(1.0) for v in scee_values
    ), "No SCEE=1.0 entries found; GLYCAM glycan dihedrals were not written correctly"
    assert any(
        v == pytest.approx(1.2, rel=1e-3) for v in scee_values
    ), "No SCEE=1.2 entries found; standard AMBER protein dihedrals are missing"
    assert any(
        v == pytest.approx(1.0) for v in scnb_values
    ), "No SCNB=1.0 entries found; GLYCAM glycan dihedrals were not written correctly"
    assert any(
        v == pytest.approx(2.0, rel=1e-3) for v in scnb_values
    ), "No SCNB=2.0 entries found; standard AMBER protein dihedrals are missing"

    # SCEE/SCNB=0.0 is valid and expected — it marks dihedrals that share
    # terminal atoms with another dihedral and should not contribute a 1-4 term.
    # The pre-fix bug was that ALL glycan dihedral entries were 0.0 because
    # CLJScaleFactor(1.0, 1.0) pairs were silently dropped. Having both 1.0 and
    # 1.2 present (checked above) confirms the fix is working correctly.

    # Reload and verify the energy is self-consistent after the AMBER roundtrip.
    # Before the fix, glycan 1-4 pairs had SCEE=0/SCNB=0, giving zero 1-4
    # interactions and a different energy.
    mols2 = sr.load(f)
    assert mols2.energy().value() == pytest.approx(mols.energy().value(), rel=1e-3)
