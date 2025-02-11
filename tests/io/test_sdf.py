import sire as sr


def test_cis_or_trans():
    """
    Make sure that we can read an SDF file containing a cis or trans double bond.
    """
    mols = sr.load_test_files("cis_trans_double_bond.sdf")


def test_positive_charge():
    """
    Make sure that we correctly set postive atomic formal charges on
    read, then re-convert to the correct SDF format on write.
    """
    mol = sr.load_test_files("positive_charge.sdf")[0]

    from math import isclose

    import tempfile

    # Make sure that the charge is set correctly.
    # A value of 3 in the file should be converted to a charge of 1.
    assert isclose(mol.charge().value(), 1.0)

    # Write back to file.
    with tempfile.NamedTemporaryFile(suffix=".sdf") as f:
        sr.save(mol, f.name)

        # Read back in and check that the charge is still 1.
        mol = sr.load(f.name)[0]
        assert isclose(mol.charge().value(), 1.0)
