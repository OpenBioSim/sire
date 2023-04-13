import sire as sr


def test_no_label():
    # try to load these files - they segfault if the label is missing
    # Need to fix this bug :-)
    mols = sr.load_test_files("no_label.rst7", "no_label.prm7")
