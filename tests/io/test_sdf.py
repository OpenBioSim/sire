import sire as sr


def test_cis_or_trans():
    """
    Make sure that we can read an SDF file containing a cis or trans double bond.
    """
    mols = sr.load_test_files("cis_trans_double_bond.sdf")


def test_charge():
    """
    Make sure that the SDF formal charge field is parsed correctly.
    """
    mol = sr.load_test_files("charge.sdf")[0]

    from math import isclose

    import tempfile

    # SDF format has a weird mapping for formal charge. Here the keys are
    # the SDF value, and the values are the formal charge.
    mapping = {
        7: -3,
        6: -2,
        5: -1,
        0: 0,
        3: 1,
        2: 2,
        1: 3,
    }

    # Make sure the charges are parsed correctly.
    for c0, c1 in zip(mol.property("formal_charge").to_list(), mapping.values()):
        assert isclose(c0.value(), c1)

    # Write back to file.
    with tempfile.NamedTemporaryFile(suffix=".sdf") as f:
        sr.save(mol, f.name)

        # Read back in and check that the charges are still correct.
        for c0, c1 in zip(mol.property("formal_charge").to_list(), mapping.values()):
            assert isclose(c0.value(), c1)


def test_name(tagged_sdf):
    """
    Make sure that the molecule takes its name from the SDF title field.
    """
    assert tagged_sdf.name().value() == tagged_sdf.property("name").value()
