import sire as sr

import pytest


def test_grospace(tmpdir, kigaki_mols):
    mols = kigaki_mols.clone()

    energy = mols.energy()

    mol = mols[0]

    mol = mol.edit().rename(" \t\t  NAME \t\t WITH     SPACE    ").commit()

    mols.update(mol)

    d = tmpdir.mkdir("test_grospace")

    f = sr.save(mols, d.join("test"), format=["Gro87", "GroTop"])

    lines = open(f[1], "r").readlines()

    found_correct_line = False

    for line in lines:
        # look for the name with underscores replacing spaces
        if line.find("NAME_WITH_SPACE") != -1:
            found_correct_line = True
            break

    assert found_correct_line

    mols = sr.load(f, show_warnings=False)

    assert energy.value() == pytest.approx(mols.energy().value())


def test_posre():
    # Make sure we can parse a file with BioSimSpace position restraint include
    # directives.
    mols = sr.load_test_files("posre.top")


def test_fep_atoms():
    """
    Test that GROMACS FEP atomtypes and atoms are created correctly.
    """

    import os
    import tempfile

    atomtypes = """[ atomtypes ]
    ; name      at.num        mass      charge   ptype       sigma     epsilon
        C1           6   12.010700    0.000000       A    0.348065    0.363503
    C1_du           0    0.000000    0.000000       A    0.348065    0.000000
        C2           6   12.010700    0.000000       A    0.337953    0.455389
        H1           1    1.007940    0.000000       A    0.110343    0.058956
    H1_du           0    0.000000    0.000000       A    0.110343    0.000000
        H2           1    1.007940    0.000000       A    0.257258    0.065318
    H2_du           0    0.000000    0.000000       A    0.264454    0.000000
    H2_dux           0    0.000000    0.000000       A    0.257258    0.000000
        H3           1    1.007940    0.000000       A    0.264454    0.066021
    H3_du           0    0.000000    0.000000       A    0.258323    0.000000
        H4           1    1.007940    0.000000       A    0.258323    0.068656
        H5           1    1.007940    0.000000       A    0.245363    0.054840
        N1           7   14.006700    0.000000       A    0.320688    0.701621
        O1           8   15.999400    0.000000       A    0.303981    0.879502
    O1_du           0    0.000000    0.000000       A    0.303981    0.000000
        O2           8   15.999400    0.000000       A    0.302511    0.704858
    """

    atoms = """[ atoms ]
    ;   nr   type0  resnr residue  atom   cgnr    charge0        mass0   type1    charge1        mass1
        1      O1      1     LIG   O1x      1  -0.803190    15.999430   O1_du   0.000000    15.999430
        2      C1      1     LIG   C1x      2   0.916780    12.010780      N1  -0.662150    14.006720
        3      O1      1     LIG   O2x      3  -0.803190    15.999430   O1_du   0.000000    15.999430
        4      C1      1     LIG   C2x      4   0.082410    12.010780      C1  -0.105710    12.010780
        5      N1      1     LIG   N1x      5  -0.564860    14.006720      N1  -0.484790    14.006720
        6      H1      1     LIG   H1x      6   0.449360     1.007947      H1   0.407460     1.007947
        7      C1      1     LIG   C3x      7   0.061040    12.010780      C1   0.001920    12.010780
        8      C1      1     LIG   C4x      8  -0.087870    12.010780      C1  -0.087450    12.010780
        9      C1      1     LIG   C5x      9  -0.101450    12.010780      N1  -0.015090    14.006720
        10      C1      1     LIG   C6x     10  -0.165020    12.010780      C1  -0.045730    12.010780
        11      H2      1     LIG   H2x     11   0.106400     1.007947      H5   0.196160     1.007947
        12      C1      1     LIG   C7x     12  -0.107040    12.010780   C1_du   0.000000    12.010780
        13      H2      1     LIG   H3x     13   0.120980     1.007947   H2_dux   0.000000     1.007947
        14      C1      1     LIG   C8x     14  -0.057540    12.010780      C1  -0.204760    12.010780
        15      C1      1     LIG   C9x     15  -0.203310    12.010780      C1   0.114810    12.010780
        16      C2      1     LIG  C10x     16   0.013180    12.010780      C2  -0.054450    12.010780
        17      H3      1     LIG   H4x     17   0.044890     1.007947      H3   0.058540     1.007947
        18      H3      1     LIG   H5x     18   0.044890     1.007947      H3   0.058540     1.007947
        19      C2      1     LIG  C11x     19  -0.084960    12.010780      C2  -0.080770    12.010780
        20      H3      1     LIG   H6x     20   0.058790     1.007947      H3   0.071370     1.007947
        21      H3      1     LIG   H7x     21   0.058790     1.007947      H3   0.071370     1.007947
        22      C2      1     LIG  C12x     22   0.126550    12.010780      C2   0.130690    12.010780
        23      H4      1     LIG   H8x     23   0.037860     1.007947      H4   0.038330     1.007947
        24      H4      1     LIG   H9x     24   0.037860     1.007947      H4   0.038330     1.007947
        25      O2      1     LIG   O3x     25  -0.332540    15.999430      O2  -0.334080    15.999430
        26      C1      1     LIG  C13x     26   0.142880    12.010780      C1   0.108960    12.010780
        27      C1      1     LIG  C14x     27  -0.181170    12.010780      C1  -0.171770    12.010780
        28      H2      1     LIG  H10x     28   0.144890     1.007947      H2   0.138510     1.007947
        29      C1      1     LIG  C15x     29  -0.052890    12.010780      C1  -0.042000    12.010780
        30      C2      1     LIG  C16x     30  -0.051460    12.010780      C2  -0.059900    12.010780
        31      H3      1     LIG  H11x     31   0.040320     1.007947      H3   0.050460     1.007947
        32      H3      1     LIG  H12x     32   0.040320     1.007947      H3   0.050460     1.007947
        33      H3      1     LIG  H13x     33   0.040320     1.007947      H3   0.050460     1.007947
        34      C1      1     LIG  C17x     34  -0.175410    12.010780      C1  -0.147890    12.010780
        35      H2      1     LIG  H14x     35   0.122060     1.007947      H2   0.142460     1.007947
        36      C1      1     LIG  C18x     36  -0.101660    12.010780      C1  -0.094400    12.010780
        37      H2      1     LIG  H15x     37   0.123170     1.007947      H2   0.139870     1.007947
        38      C1      1     LIG  C19x     38  -0.184410    12.010780      C1  -0.174770    12.010780
        39      H2      1     LIG  H16x     39   0.144680     1.007947      H2   0.138300     1.007947
        40      C2      1     LIG  C20x     40  -0.038120    12.010780      C2  -0.032020    12.010780
        41      H3      1     LIG  H17x     41   0.030840     1.007947   H2_du   0.000000     1.007947
        42      H3      1     LIG  H18x     42   0.030840     1.007947   H2_du   0.000000     1.007947
        43      H3      1     LIG  H19x     43   0.030840     1.007947   H2_du   0.000000     1.007947
        44      C2      1     LIG  C21x     44  -0.035100    12.010780      C2   0.001470    12.010780
        45      H3      1     LIG  H20x     45   0.026750     1.007947   H2_du   0.000000     1.007947
        46      H3      1     LIG  H21x     46   0.026750     1.007947   H2_du   0.000000     1.007947
        47      H3      1     LIG  H22x     47   0.026750     1.007947   H2_du   0.000000     1.007947
        48   H1_du      1     LIG   H1x     48   0.000000     1.007947      H1   0.459340     1.007947
        49   H1_du      1     LIG   H2x     49   0.000000     1.007947      H1   0.459340     1.007947
        50   H1_du      1     LIG   H3x     50   0.000000     1.007947      H1   0.459340     1.007947
        51   H2_du      1     LIG  H19x     51   0.000000     1.007947      H3   0.062430     1.007947
        52   H2_du      1     LIG  H20x     52   0.000000     1.007947      H3   0.062430     1.007947
        53   H2_du      1     LIG  H21x     53   0.000000     1.007947      H3   0.062430     1.007947
        54   H3_du      1     LIG  H22x     54   0.000000     1.007947      H4   0.074650     1.007947
        55   H3_du      1     LIG  H23x     55   0.000000     1.007947      H4   0.074650     1.007947
        56   H3_du      1     LIG  H24x     56   0.000000     1.007947      H4   0.074650     1.007947
    """

    # Load the molecule.
    merged_mol = sr.load_test_files("merged_molecule_grotop.s3")

    # Save to a temporary file.
    with tempfile.TemporaryDirectory() as tmpdir:
        sr.save(merged_mol, os.path.join(tmpdir, "merged_molecule"), format="GroTop")

        # Read the [ atomtypes ] and [ atoms ] sections from the file.
        with open(os.path.join(tmpdir, "merged_molecule.grotop"), "r") as f:
            atomtypes_lines = []
            atoms_lines = []
            found_atomtypes = False
            found_atoms = False
            for line in f:
                if line.startswith("[ atomtypes ]"):
                    found_atomtypes = True
                elif line.startswith("[ atoms ]"):
                    found_atoms = True

                if found_atomtypes and not found_atoms:
                    if line != "\n":
                        atomtypes_lines.append(line)
                    else:
                        found_atomtypes = False

                if found_atoms:
                    if line != "\n":
                        atoms_lines.append(line)
                    else:
                        found_atoms = False

        # Check that the atomtypes and atoms are as expected.
        atomtypes_list = atomtypes.split("\n")
        for a, b in zip(atomtypes_list, atomtypes_lines):
            assert a.strip() == b.strip()
        atoms_list = atoms.split("\n")
        for a, b in zip(atoms_list, atoms_lines):
            assert a.strip() == b.strip()


def test_grotop_cmap(tmpdir, gromacs_cmap):
    """
    Test that the GROMACS cmap is correctly parsed and saved.
    """
    mols = gromacs_cmap.clone()

    dir = tmpdir.mkdir("test_grotop_cmap")

    mol = mols[0]
    cmaps = mol.property("cmap").parameters()

    assert len(cmaps) == 127

    unique_cmaps = {}

    # these parameters are all between atoms called "C", "N", "CA", "C", "N"
    for cmap in cmaps:
        assert mol.atom(cmap.atom0()).name().value() == "C"
        assert mol.atom(cmap.atom1()).name().value() == "N"
        assert mol.atom(cmap.atom2()).name().value() == "CA"
        assert mol.atom(cmap.atom3()).name().value() == "C"
        assert mol.atom(cmap.atom4()).name().value() == "N"

        # there should be 24 rows and 24 columns in the cmap
        assert cmap.parameter().num_rows() == 24
        assert cmap.parameter().num_columns() == 24

        unique_cmaps[cmap.parameter().to_string()] = 1

    assert len(unique_cmaps) == 3

    orig_cmaps = cmaps

    # now other molecule should have any cmap parameters
    for mol in mols[1:]:
        try:
            cmaps = mol.property("cmap")
        except Exception:
            cmaps = None

        if cmaps is not None:
            assert len(cmaps.parameters()) == 0

    # Save to a temporary file.
    f = sr.save(mols, dir.join("output"), format="GroTop")

    # Load the saved file.
    mols2 = sr.load(f, show_warnings=False)

    # assert that we have the same cmap parameters
    cmaps = mols2[0].property("cmap").parameters()

    assert len(cmaps) == 127

    unique_cmaps2 = {}

    # these parameters are all between atoms called "C", "N", "CA", "C", "N"
    for cmap in cmaps:
        assert mols2[0].atom(cmap.atom0()).name().value() == "C"
        assert mols2[0].atom(cmap.atom1()).name().value() == "N"
        assert mols2[0].atom(cmap.atom2()).name().value() == "CA"
        assert mols2[0].atom(cmap.atom3()).name().value() == "C"
        assert mols2[0].atom(cmap.atom4()).name().value() == "N"

        # there should be 24 rows and 24 columns in the cmap
        assert cmap.parameter().num_rows() == 24
        assert cmap.parameter().num_columns() == 24

        unique_cmaps2[cmap.parameter().to_string()] = 1

    assert len(unique_cmaps2) == 3

    assert cmaps == orig_cmaps

    # make sure that no other molecule have any cmap parameters
    for mol in mols2[1:]:
        try:
            cmaps = mol.property("cmap")
        except Exception:
            cmaps = None

        if cmaps is not None:
            assert len(cmaps.parameters()) == 0
