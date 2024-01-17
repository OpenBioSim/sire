import math
import pytest
import tempfile

from sire.legacy.Convert._SireOpenMM import EMLECallback, EMLEEngine

from sire.qm import emle

try:
    from emle.calculator import EMLECalculator

    has_emle = True
except:
    has_emle = False


def test_callback():
    """Makes sure that a callback method works correctly"""

    class Test:
        def callback(self, a, b, c, d):
            return (42, d, c)

    # Instantiate the class.
    test = Test()

    # Create a callback object.
    cb = EMLECallback(test, "callback")

    # Create some lists to hold test data.
    a = [1, 2]
    b = [3, 4]
    c = [a, b]
    d = [b, a]

    # Call the callback.
    result = cb.call(a, b, c, d)

    # Make sure the result is correct.
    assert result == (42, d, c) == test.callback(a, b, c, d)


@pytest.mark.skipif(not has_emle, reason="emle-engine is not installed")
def test_interpolate(ala_mols):
    """
    Make sure that lambda interpolation between pure MM and EMLE potentials works.
    """

    # Create a local copy of the test system.
    mols = ala_mols

    # Create an EMLE calculator.
    calculator = EMLECalculator(device="cpu", log=0, save_settings=False)

    # Create a dynamics object.
    d = mols.dynamics(timestep="1fs", constraint="none", platform="cpu")

    # Get the pure MM energy.
    nrg_mm = d.current_potential_energy()

    # Create an EMLE engine bound to the calculator.
    mols, engine = emle(mols, mols[0], calculator)

    # Create a QM/MM capable dynamics object.
    d = mols.dynamics(
        timestep="1fs", constraint="none", qm_engine=engine, platform="cpu"
    )

    # Get the pure EMLE energy.
    nrg_emle = d.current_potential_energy()

    # Get interpolated MM energy.
    d.set_lambda(0.0)
    nrg_mm_interp = d.current_potential_energy()

    # Make sure this agrees with the standard MM energy.
    assert math.isclose(nrg_mm_interp.value(), nrg_mm.value(), rel_tol=1e-4)

    # Now get the interpolated energy at lambda = 0.5.
    d.set_lambda(0.5)
    nrg_interp = d.current_potential_energy()

    # Make sure the interpolated energy is correct.
    assert math.isclose(
        nrg_interp.value(), 0.5 * (nrg_mm + nrg_emle).value(), rel_tol=1e-4
    )
