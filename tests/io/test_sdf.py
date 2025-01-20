import sire as sr


def test_cis_or_trans():
    """
    Make sure that we can read an SDF file containing a cis or trans double bond.
    """
    mols = sr.load_test_files("cis_trans_double_bond.sdf")
