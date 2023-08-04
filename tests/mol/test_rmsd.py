import sire as sr
import pytest


def test_rmsd(ala_trr, ala_mols):
    mols = ala_trr.clone()
    mols2 = ala_mols.clone()

    rmsds = mols.trajectory().rmsd(mols[0])

    expected = [
        0,
        0.1537,
        0.1417,
        0.1288,
        0.1189,
        0.1372,
        0.1188,
        0.1607,
        0.1097,
        0.1160,
        0.1504,
    ]

    for r, e in zip(rmsds, expected):
        rmsd = r.to(sr.units.nanometer)
        assert rmsd == pytest.approx(e, 1e-3)

    # try frame 3
    expected = [
        0.1288,
        0.1078,
        0.0945,
        0.0000,
        0.1142,
        0.1187,
        0.0927,
        0.1010,
        0.1098,
        0.1126,
        0.1261,
    ]

    rmsds = mols.trajectory().rmsd(mols[0], frame=3)

    for r, e in zip(rmsds, expected):
        rmsd = r.to(sr.units.nanometer)
        assert rmsd == pytest.approx(e, 1e-3)

    mapping = sr.match_atoms(mols2[0], mols[0])

    rmsds = mols.trajectory().rmsd(mols2[0], mapping=mapping, match_all=False)

    for rmsd in rmsds:
        rmsd = r.to(sr.units.angstrom)

        assert rmsd < 1.5
