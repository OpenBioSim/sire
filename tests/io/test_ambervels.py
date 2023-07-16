import pytest
import sire as sr


def _assert_same_vel(v1, v2, precision):
    for i in range(0, 3):
        x1 = v1[i]
        x2 = v2[i]

        assert x1.has_same_units(x2)
        assert x1.value() == pytest.approx(x2.value(), precision)


@pytest.mark.veryslow
def test_amber_vels_compare(tmpdir, triclinic_protein, triclinic_protein_rst7):
    mols = triclinic_protein
    mols7 = triclinic_protein_rst7

    # check that the velocities are the same between these two files
    for mol, mol7 in zip(mols[0::100], mols7[0::100]):
        for atom, atom7 in zip(mol.atoms()[0::20], mol7.atoms()[0::20]):
            _assert_same_vel(
                atom.property("velocity"),
                atom7.property("velocity"),
                precision=1e-3,
            )

    # now write out to a rst7 file and rst file, then compare those
    d = tmpdir.mkdir("test_amber_vels_compare")

    f = sr.save(mols, d.join("test"), format=["PRMTOP", "RST7"])
    f7 = sr.save(mols7, d.join("test"), format=["RST"])

    check = sr.load(f[0], f[1])
    check7 = sr.load(f[0], f7[0])

    for mol, mol7 in zip(check[0::100], mols7[0::100]):
        for atom, atom7 in zip(mol.atoms()[0::20], mol7.atoms()[0::20]):
            _assert_same_vel(
                atom.property("velocity"),
                atom7.property("velocity"),
                precision=1e-3,
            )

    for mol, mol7 in zip(mols[0::100], check7[0::100]):
        for atom, atom7 in zip(mol.atoms()[0::20], mol7.atoms()[0::20]):
            _assert_same_vel(
                atom.property("velocity"),
                atom7.property("velocity"),
                precision=1e-3,
            )
