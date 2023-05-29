import pytest

import sire as sr


def test_trajectories(tmpdir, ala_traj):
    # ala_traj is a AmberRst file - it contains time and space info
    mols = ala_traj

    d = tmpdir.mkdir("test_trajectories")

    f = sr.save(
        mols.trajectory()[0:5],
        d.join("test"),
        format=["PRMTOP", "TRAJ", "DCD", "TRR", "XTC", "RST"],
    )

    def assert_space_equal(box0, box1, precision):
        dims0 = box0.dimensions()
        dims1 = box1.dimensions()

        assert dims0.x().value() == pytest.approx(dims1.x().value(), precision)
        assert dims0.y().value() == pytest.approx(dims1.y().value(), precision)
        assert dims0.z().value() == pytest.approx(dims1.z().value(), precision)

    def assert_coords_equal(mol0, mol1, precision):
        coords0 = mol0.atoms()[3].coordinates()
        coords1 = mol1.atoms()[3].coordinates()

        assert coords0.x().value() == pytest.approx(
            coords1.x().value(), precision
        )
        assert coords0.y().value() == pytest.approx(
            coords1.y().value(), precision
        )
        assert coords0.z().value() == pytest.approx(
            coords1.z().value(), precision
        )

        coords0 = mol0.atoms()[-1].coordinates()
        coords1 = mol1.atoms()[-1].coordinates()

        assert coords0.x().value() == pytest.approx(
            coords1.x().value(), precision
        )
        assert coords0.y().value() == pytest.approx(
            coords1.y().value(), precision
        )
        assert coords0.z().value() == pytest.approx(
            coords1.z().value(), precision
        )

    for traj in f[1:]:
        print(traj)

        precision = 0.001
        check_time = True
        time_delta = 0

        if traj.endswith(".xtc"):
            precision = 0.01

        if traj.endswith(".traj"):
            # this format doesn't store the simulation time
            check_time = False
        elif traj.endswith(".dcd"):
            # this format only stores a single timestep, and the
            # first frame is 0 ps (the original trajectory starts
            # at 0.2 ps)
            time_delta = 0.2 * sr.units.picosecond

        mols2 = sr.load(f[0], traj)

        assert mols2.num_frames() == 5

        m = mols.trajectory()[3].current()
        m2 = mols2.trajectory()[3].current()

        assert m.num_molecules() == m2.num_molecules()

        if check_time:
            assert m.property("time").value() == pytest.approx(
                m2.property("time").value() + time_delta, precision
            )
        else:
            assert m2.property("time").value() == 0

        assert_space_equal(
            m.property("space"), m2.property("space"), precision
        )

        assert_coords_equal(m[0], m2[0], precision)
