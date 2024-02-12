import math
import pytest
import tempfile

from sire.legacy.Convert import EMLECallback, EMLEEngine

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


@pytest.mark.parametrize(
    "selection, expected",
    [
        (
            "residx 0",
            (
                {6: 4},
                {6: [7, 8]},
                {6: 0.8164794007490638},
            ),
        ),
        (
            "residx 1",
            (
                {4: 6, 16: 14},
                {4: [1, 5], 16: [17, 18]},
                {4: 0.7565543071161049, 16: 0.8164794007490638},
            ),
        ),
        (
            "residx 2",
            (
                {14: 16},
                {14: [8, 15]},
                {14: 0.7565543071161049},
            ),
        ),
    ],
)
def test_link_atoms(ala_mols, selection, expected):
    """
    Make sure that the link atoms are correctly identified.
    """

    from sire.base import create_map as _create_map
    from sire.qm._utils import _create_qm_mol_to_atoms, _get_link_atoms

    # Create a local copy of the test system.
    mols = ala_mols

    # Extract the QM atom selection.
    qm_atoms = mols[0][selection].atoms()

    # Create the mapping between molecule numbers and QM atoms.
    qm_mol_to_atoms = _create_qm_mol_to_atoms(qm_atoms)

    # Get link atom information.
    mm1_to_qm, mm1_to_mm2, bond_scale_factors, mm1_indices = _get_link_atoms(
        mols, qm_mol_to_atoms, _create_map({})
    )

    assert mm1_to_qm == expected[0]
    assert mm1_to_mm2 == expected[1]
    assert bond_scale_factors == expected[2]


@pytest.mark.skipif(not has_emle, reason="emle-engine is not installed")
@pytest.mark.parametrize("selection", ["molidx 0", "resname ALA"])
def test_interpolate(ala_mols, selection):
    """
    Make sure that lambda interpolation between pure MM and EMLE potentials works.
    """

    # Create a local copy of the test system.
    mols = ala_mols.__copy__()

    # Create an EMLE calculator.
    calculator = EMLECalculator(device="cpu")

    # Create a dynamics object.
    d = mols.dynamics(timestep="1fs", constraint="none", platform="cpu")

    # Get the pure MM energy.
    nrg_mm = d.current_potential_energy()

    # Create an EMLE engine bound to the calculator.
    mols, engine = emle(mols, selection, calculator)

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
