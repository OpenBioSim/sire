import pytest

import sire as sr


def assert_space_equal(box0, box1, precision):
    dims0 = box0.dimensions()
    dims1 = box1.dimensions()

    assert dims0.x().value() == pytest.approx(dims1.x().value(), precision)
    assert dims0.y().value() == pytest.approx(dims1.y().value(), precision)
    assert dims0.z().value() == pytest.approx(dims1.z().value(), precision)


def assert_coords_equal(mol0, mol1, precision):
    coords0 = mol0.atoms()[3].coordinates()
    coords1 = mol1.atoms()[3].coordinates()

    assert coords0.x().value() == pytest.approx(coords1.x().value(), precision)
    assert coords0.y().value() == pytest.approx(coords1.y().value(), precision)
    assert coords0.z().value() == pytest.approx(coords1.z().value(), precision)

    coords0 = mol0.atoms()[-1].coordinates()
    coords1 = mol1.atoms()[-1].coordinates()

    assert coords0.x().value() == pytest.approx(coords1.x().value(), precision)
    assert coords0.y().value() == pytest.approx(coords1.y().value(), precision)
    assert coords0.z().value() == pytest.approx(coords1.z().value(), precision)


def test_ignore_topology_frame(tmpdir, ala_traj):
    mols = ala_traj.clone()

    assert mols.num_frames() == 500

    mol = mols[0]

    mols.load_frame(100)

    d = tmpdir.mkdir("test_ignore_topology_frame")

    pdb = sr.save(mols, d.join("test"), format=["PDB"])

    rst = sr.save(mols.trajectory()[0:10], d.join("test"), format=["RST"])

    # load including the PDB frame... - the should be the same
    # as the 100th frame (as this was saved as the PDB)
    mols2 = sr.load(pdb + rst)

    assert mols2.num_frames() == 11
    assert_coords_equal(mols[0], mols2[0], 0.001)

    # load ignoring the PDB frame - this should be the same as the loaded
    # molecule that was frame 0 of the TRAJ file
    mols2 = sr.load(pdb + rst, ignore_topology_frame=True)

    assert mols2.num_frames() == 10
    assert_coords_equal(mol, mols2[0], 0.001)


def test_trajectories(tmpdir, ala_traj):
    # ala_traj is a AmberRst file - it contains time and space info
    mols = ala_traj.clone()

    d = tmpdir.mkdir("test_trajectories")

    f = []

    for format in ["PRMTOP", "TRAJ", "DCD", "TRR", "XTC", "RST"]:
        f.append(sr.save(mols.trajectory()[0:5], d.join("test"), format=[format])[0])

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

        print(f"LOAD {f[0]} + {traj}")
        mols2 = sr.load(f[0], traj)

        assert mols2.num_frames() == 5

        m = mols.trajectory()[3].current()
        m2 = mols2.trajectory()[3].current()

        print(m.property("time"))
        print(m2.property("time"))

        assert m.num_molecules() == m2.num_molecules()

        if check_time:
            correct_time = m.property("time").value() == pytest.approx(
                m2.property("time").value() + time_delta, precision
            )

            if not correct_time:
                print("NOT CORRECT TIME!")
                print(m.property("time"))
                print(m2.property("time"))
                m2 = mols2.trajectory()[3].current()
                print(m2.property("time"))
                mols2 = sr.load(f[0], traj)
                m2 = mols2.trajectory()[3].current()

                for frame in mols2.trajectory():
                    print(frame.property("time"))

                assert m.property("time").value() == pytest.approx(
                    m2.property("time").value() + time_delta, precision
                )

        else:
            assert m2.property("time").value() == 0

        assert_space_equal(m.property("space"), m2.property("space"), precision)

        assert_coords_equal(m[0], m2[0], precision)


def test_frame_names(tmpdir, ala_traj):
    mols = ala_traj.clone()

    d = tmpdir.mkdir("test_frame_names")

    # First, write two trajectory frames to named Gro87 files.
    sr.save(mols.trajectory()[0:2], [d.join("cat.gro"), d.join("dog.gro")])

    # Make sure the output files exist.
    cat_path = d.join("cat.gro")
    dog_path = d.join("dog.gro")
    assert cat_path.check()
    assert dog_path.check()

    # Try writing to mixed formats. This should raise a ValueError.
    with pytest.raises(ValueError):
        sr.save(
            mols.trajectory()[0:2],
            [d.join("cat.gro"), d.join("dog.rst7")],
        )

    # Try writing with mismatched number of frames and filenames.
    with pytest.raises(ValueError):
        sr.save(
            mols.trajectory()[0:3],
            [d.join("cat.gro"), d.join("dog.gro")],
        )
