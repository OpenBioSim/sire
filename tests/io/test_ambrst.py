import sire as sr
import pytest


@pytest.mark.slow
def test_ambrst_velocities(tmpdir, kigaki_mols):
    mols = kigaki_mols.clone()

    mol = mols[10]

    a_per_ps = sr.units.angstrom / sr.units.picosecond

    v = sr.legacy.Mol.Velocity3D(1 * a_per_ps, 2 * a_per_ps, 3 * a_per_ps)

    c = mol.cursor().atoms()

    c["velocity"] = v

    mol = c.molecule().commit()

    mols.update(mol)

    d = tmpdir.mkdir("test_ambrst_velocities")

    f = sr.save(mols, d.join("test"), format=["RST", "PRM7"])

    mols = sr.load(f, show_warnings=False)

    def compare(v1, v2):
        assert v1.x().value() == pytest.approx(v2.x().value())
        assert v1.y().value() == pytest.approx(v2.y().value())
        assert v1.z().value() == pytest.approx(v2.z().value())

    for i, mol in enumerate(mols[0:20]):
        if i == 10:
            for atom in mol.atoms():
                compare(atom.property("velocity"), v)
        else:
            for atom in mol.atoms():
                compare(atom.property("velocity"), sr.legacy.Mol.Velocity3D())

    for i, mol in enumerate(mols[30::100]):
        for atom in mol.atoms():
            compare(atom.property("velocity"), sr.legacy.Mol.Velocity3D())
