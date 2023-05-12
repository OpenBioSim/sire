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
def h7n9_mols():
    return sr.load_test_files("h7n9.pdb", "h7n9.dcd")


@pytest.fixture(scope="session")
def chol_mols():
    return sr.load_test_files("cholesterol.sdf")


@pytest.fixture(scope="session")
def ala_traj():
    return sr.load_test_files("ala.top", "ala.traj")


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
def openmm_interchange_mols():
    return sr.load_test_files("openmm_interchange.rst7", "openmm_interchange.prm7")
