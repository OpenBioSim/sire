import sire as sr
import pytest


def test_ambrst_vels_and_forces(tmpdir, kigaki_mols):
    mols = kigaki_mols.clone()

    mols = mols[5:10]

    mol = mols[2]

    a_per_ps = sr.units.angstrom / sr.units.picosecond

    v = sr.legacy.Mol.Velocity3D(1 * a_per_ps, 2 * a_per_ps, 3 * a_per_ps)

    c = mol.cursor().atoms()

    c["velocity"] = v

    mol = c.molecule().commit()

    mols.update(mol)

    mol = mols[1]

    kcal_per_a = sr.units.kcal / sr.units.angstrom

    force = sr.legacy.Mol.Force3D(1 * kcal_per_a, 2 * kcal_per_a, 3 * kcal_per_a)

    c = mol.cursor().atoms()

    c["force"] = force

    mol = c.molecule().commit()

    mols.update(mol)

    d = tmpdir.mkdir("test_ambrst_vels_and_forces")

    f = sr.save(mols, d.join("test"), format=["RST", "PRM7"])

    mols = sr.load(f, show_warnings=False)

    def compare(v1, v2):
        assert v1.x().value() == pytest.approx(v2.x().value())
        assert v1.y().value() == pytest.approx(v2.y().value())
        assert v1.z().value() == pytest.approx(v2.z().value())

    for i, mol in enumerate(mols):
        if i == 1:
            for atom in mol.atoms():
                compare(atom.property("velocity"), sr.legacy.Mol.Velocity3D())
                compare(atom.property("force"), force)
        elif i == 2:
            for atom in mol.atoms():
                compare(atom.property("velocity"), v)
                compare(atom.property("force"), sr.legacy.Mol.Force3D())
        else:
            for atom in mol.atoms():
                compare(atom.property("velocity"), sr.legacy.Mol.Velocity3D())
                compare(atom.property("force"), sr.legacy.Mol.Force3D())
