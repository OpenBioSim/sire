import sire as sr
import pytest


def pytest_addoption(parser):
    parser.addoption(
        "--runslow", action="store_true", default=False, help="run slow tests"
    )
    parser.addoption(
        "--runveryslow",
        action="store_true",
        default=False,
        help="run slow and veryslow tests",
    )


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "slow: mark tests as slow to run (takes more than a couple of seconds",
    )
    config.addinivalue_line(
        "markers",
        "veryslow: mark tests as veryslow to run (takes more than 5-10 seconds to run",
    )


def pytest_collection_modifyitems(config, items):
    runslow = False
    runveryslow = False

    if config.getoption("--runslow"):
        runslow = True

    if config.getoption("--runveryslow"):
        runslow = True
        runveryslow = True

    if not runslow:
        skip_slow = pytest.mark.skip(reason="need --runslow option to run")
        for item in items:
            if "slow" in item.keywords:
                item.add_marker(skip_slow)

    if not runveryslow:
        skip_veryslow = pytest.mark.skip(reason="need --runveryslow option to run")
        for item in items:
            if "veryslow" in item.keywords:
                item.add_marker(skip_veryslow)


@pytest.fixture(scope="session")
def ala_mols():
    return sr.load_test_files("ala.top", "ala.crd")


@pytest.fixture(scope="session")
def ose_mols():
    return sr.load_test_files("ose.top", "ose.crd")


@pytest.fixture(scope="session")
def h7n9_mols():
    return sr.load_test_files("h7n9.pdb", "h7n9.dcd")


@pytest.fixture(scope="session")
def chol_mols():
    return sr.load_test_files("cholesterol.sdf")


@pytest.fixture(scope="session")
def ala_traj():
    return sr.load_test_files("ala.top", "ala.traj")


@pytest.fixture(scope="session")
def ala_trr():
    return sr.load_test_files("ala.top", "ala.trr")


@pytest.fixture(scope="session")
def p38_mols():
    return sr.load_test_files("p38.pdb")


@pytest.fixture(scope="session")
def alanin_mols():
    return sr.load_test_files("alanin.psf")


@pytest.fixture(scope="session")
def kigaki_mols():
    return sr.load_test_files("kigaki.gro", "kigaki.top")


@pytest.fixture(scope="session")
def neura_mols():
    return sr.load_test_files("proteinbox.crd", "proteinbox.top")


@pytest.fixture(scope="session")
def excluded_mols():
    return sr.load_test_files("excluded.rst7", "excluded.prm7")


@pytest.fixture(scope="session")
def wrapped_mols():
    return sr.load_test_files("wrapped.rst7", "wrapped.prm7")


@pytest.fixture(scope="session")
def openmm_interchange_mols():
    return sr.load_test_files("openmm_interchange.rst7", "openmm_interchange.prm7")


@pytest.fixture(scope="session")
def triclinic_protein():
    return sr.load_test_files("triclinic_protein.prm7", "triclinic_protein.crd")


@pytest.fixture(scope="session")
def triclinic_protein_rst7():
    return sr.load_test_files("triclinic_protein.prm7", "triclinic_protein.rst")


@pytest.fixture(scope="session")
def merged_ethane_methanol():
    return sr.load_test_files("merged_molecule.s3")


@pytest.fixture(scope="session")
def merged_zan_ose():
    return sr.load_test_files("merged_ligand.s3")


@pytest.fixture(scope="session")
def ethane_12dichloroethane():
    return sr.load_test_files("ethane_12dichloroethane.bss")


@pytest.fixture(scope="session")
def pentane_cyclopentane():
    return sr.load_test_files("pentane_cyclopentane.bss")


@pytest.fixture(scope="session")
def neopentane_methane():
    return sr.load_test_files("neo_meth_scratch.bss")


@pytest.fixture(scope="session")
def zero_lj_mols():
    return sr.load_test_files("zero_lj.prm7", "zero_lj.rst7")
