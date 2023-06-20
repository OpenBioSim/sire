import sire as sr
import pytest


def test_alignment(ala_mols):
    mols = ala_mols

    mol = mols[0]

    coords = mol.property("coordinates").to_vector()

    mol = (
        mol.cursor()
        .translate((10, 9, 8))
        .rotate(30 * sr.units.degrees, (1, 2, 3))
        .commit()
    )

    moved_coords = mol.property("coordinates").to_vector()

    transform = sr.maths.get_alignment(coords, moved_coords)

    def _assert_equal(c0, c1):
        assert len(c0) == len(c1)

        for a, b in zip(c0, c1):
            assert a.x().value() == pytest.approx(b.x().value())
            assert a.y().value() == pytest.approx(b.y().value())
            assert a.z().value() == pytest.approx(b.z().value())

    _assert_equal(coords, transform.apply(moved_coords))

    aligned_coords = sr.maths.align(coords, moved_coords)

    _assert_equal(coords, aligned_coords)

    _assert_equal(moved_coords, transform.inverse().apply(coords))

    _assert_equal(moved_coords, transform.reverse(coords))


def test_trajectory_alignment(ala_traj):
    mols = ala_traj.clone()

    mol = mols[0]

    t = mols.trajectory()

    unaligned_rmsds = []
    aligned_rmsds = []
    transforms = []

    # align the molecules manually
    for frame in t[1:11]:
        # calculate the first RMSD, just so we can be sure we have
        # reduced it!
        rmsd = mol.evaluate().rmsd(frame[0])

        unaligned_rmsds.append(rmsd.value())

        # calculate the transformation needed to align
        transform = sr.mol.get_alignment(mol, frame[0])

        transforms.append(transform)

        aligned_mol = frame[0].cursor().transform(transform).commit()

        aligned_rmsd = mol.evaluate().rmsd(aligned_mol)

        assert aligned_rmsd.value() < rmsd.value()

        aligned_rmsds.append(aligned_rmsd.value())

    # Now ask for the trajectory to be aligned against the first molecule
    t = mols.trajectory().align("molidx 0")

    for i, frame in enumerate(t[1:11]):
        aligned_mol = frame[0]

        # check the RMSD for the trajectory aligned molecule
        # is close to the manually aligned molecule
        rmsd = mol.evaluate().rmsd(aligned_mol)
        assert rmsd.value() == pytest.approx(aligned_rmsds[i])

        # now transform the trajectory aligned molecule back
        # to its original position
        mol2 = aligned_mol.cursor().transform(transforms[i].inverse()).commit()

        rmsd = mol.evaluate().rmsd(mol2)

        assert unaligned_rmsds[i] == pytest.approx(rmsd.value())


def test_smooth_alignment(ala_traj):
    mols = ala_traj

    smoothed = mols.trajectory(align="element C", smooth=20)
    aligned = mols.trajectory(align="element C")

    smoothed_mol = smoothed[0].current()[0]
    aligned_mol = aligned[0].current()[0]

    for s, a in zip(smoothed[1:11], aligned[1:11]):
        rmsd_s = smoothed_mol.evaluate().rmsd(s[0])
        rmsd_a = aligned_mol.evaluate().rmsd(a[0])

        # the smoothed RMSDs will be significantly less than the
        # aligned RMSDs
        assert rmsd_s.value() < 0.5 * rmsd_a.value()
