import numpy as np
import pytest
import tempfile

from sire.legacy.Convert import PyQMCallback

import sire as sr

try:
    from emle.calculator import EMLECalculator

    has_emle = True
except:
    has_emle = False

try:
    from openmmml import MLPotential

    has_openmm_ml = True
except:
    has_openmm_ml = False


def test_callback_method():
    """Makes sure that a callback method works correctly"""

    class Test:
        def callback(self, a, b, c, d):
            return (42, d, c)

    # Instantiate the class.
    test = Test()

    # Create a callback object.
    cb = PyQMCallback(test, "callback")

    # Create some lists to hold test data.
    a = [1, 2]
    b = [3, 4]
    c = [a, b]
    d = [b, a]

    # Call the callback.
    result = cb.call(a, b, c, d)

    # Make sure the result is correct.
    assert result == (42, d, c) == test.callback(a, b, c, d)


def test_callback_function():
    """Makes sure that a callback function works correctly"""

    def callback(a, b, c, d):
        return (42, d, c)

    # Create a callback object.
    cb = PyQMCallback(callback, "")

    # Create some lists to hold test data.
    a = [1, 2]
    b = [3, 4]
    c = [a, b]
    d = [b, a]

    # Call the callback.
    result = cb.call(a, b, c, d)

    # Make sure the result is correct.
    assert result == (42, d, c) == callback(a, b, c, d)


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


def test_charge_redistribution():
    """
    Make sure that charge redistribution works correctly.
    """

    import sire as sr

    from sire.base import create_map
    from sire.qm._utils import _check_charge
    from sire.mol import selection_to_atoms

    # A selection for the QM region. This is a subset of a TRP residue.
    selection = "residx 118 and not atomname C, CA, H, HA, N, O"

    # The inverse selection.
    not_selection = f"not ({selection})"

    # Load the AbyU test system.
    mols = sr.load_test_files("abyu.prm7", "abyu.rst7")

    # Selet the QM atoms.
    qm_atoms = selection_to_atoms(mols, selection)

    # Get the charges of both regions.
    charge0 = mols[selection].charge()
    charge1 = mols[not_selection].charge()

    # Check the charges, redistributing to the nearest integer.
    _check_charge(mols, qm_atoms, create_map({}), redistribute_charge=True)

    # Get the new charges.
    new_charge0 = mols[selection].charge()
    new_charge1 = mols[not_selection].charge()

    # Make sure the QM charge has been redistributed to the nearest integer.
    assert np.isclose(round(charge0.value()), new_charge0.value(), rtol=1e-4)

    # Make sure the remainder has beeen redistributed to the other atoms.
    assert np.isclose((charge0 + charge1).value(), new_charge1.value(), rtol=1e-4)

    # Make sure the check fails if we don't redistribute the charge.
    with pytest.raises(Exception):
        _check_charge(mols, qm_atoms, create_map({}), redistribute_charge=False)


@pytest.mark.skipif(not has_emle, reason="emle-engine is not installed")
@pytest.mark.parametrize("selection", ["molidx 0", "resname ALA"])
def test_emle_interpolate(ala_mols, selection):
    """
    Make sure that lambda interpolation between pure MM and EMLE potentials works.
    """

    # Create a local copy of the test system.
    mols = ala_mols.clone()

    # Create an EMLE calculator.
    calculator = EMLECalculator(device="cpu")

    # Create a dynamics object.
    d = mols.dynamics(timestep="1fs", constraint="none", platform="cpu")

    # Get the pure MM energy.
    nrg_mm = d.current_potential_energy()

    # Create an EMLE engine bound to the calculator.
    mols, engine = sr.qm.emle(mols, selection, calculator)

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
    assert np.isclose(nrg_mm_interp.value(), nrg_mm.value(), rtol=1e-4)

    # Now get the interpolated energy at lambda = 0.5.
    d.set_lambda(0.5)
    nrg_interp = d.current_potential_energy()

    # Make sure the interpolated energy is correct. Note that the interpolation
    # is actually non-linear so the energies are not exactly the average of the
    # two states.
    assert np.isclose(nrg_interp.value(), 0.5 * (nrg_mm + nrg_emle).value(), rtol=1e-4)


@pytest.mark.skipif(not has_emle, reason="emle-engine is not installed")
@pytest.mark.skipif(not has_openmm_ml, reason="openmm-ml is not installed")
def test_emle_openmm_ml(ala_mols):
    """
    Make sure that the EMLE engine can be used with OpenMM-ML.
    """

    import openmm
    import sire as sr

    # Create a local copy of the test system.
    mols = ala_mols.clone()

    # Create an EMLE calculator.
    calculator = EMLECalculator(backend="torchani", device="cpu")

    # Create an EMLE engine bound to the calculator.
    emle_mols, engine = sr.qm.emle(mols, mols[0], calculator)

    # Create a QM/MM capable dynamics object.
    d = emle_mols.dynamics(
        timestep="1fs",
        constraint="none",
        qm_engine=engine,
        cutoff_type="pme",
        cutoff="7.5 A",
        platform="cpu",
    )

    # Get the energy.
    nrg_sire = d.current_potential_energy().to("kJ_per_mol")

    with tempfile.TemporaryDirectory() as tmpdir:
        # Create a new EMLECalculator without a backend.
        calculator = EMLECalculator(backend=None, device="cpu")

        # Create an EMLE engine bound to the calculator.
        emle_mols, engine = sr.qm.emle(mols, mols[0], calculator)

        # Write the sytem to an AMBER coordinate and topology file.
        files = sr.expand(tmpdir, ["ala.rst7", "ala.prm7"])
        for file in files:
            sr.save(mols, file)

        # Load back the files and create an OpenMM topology.
        inpcrd = openmm.app.AmberInpcrdFile(f"{tmpdir}/ala.rst7")
        prmtop = openmm.app.AmberPrmtopFile(f"{tmpdir}/ala.prm7")
        topology = prmtop.topology

        # Create the MM system.
        mm_system = prmtop.createSystem(
            nonbondedMethod=openmm.app.PME,
            nonbondedCutoff=7.5 * openmm.unit.angstrom,
            constraints=openmm.app.HBonds,
        )

        # Define the ML region.
        ml_atoms = list(range(mols[0].num_atoms()))

        # Create the mixed ML/MM system.
        potential = MLPotential("ani2x")
        ml_system = potential.createMixedSystem(
            topology, mm_system, ml_atoms, interpolate=False
        )

        # Get the OpenMM forces from the engine.
        emle_force, interpolation_force = engine.get_forces()

        # Add the EMLE force to the system.
        ml_system.addForce(emle_force)

        # Create the integrator.
        integrator = openmm.LangevinMiddleIntegrator(
            300 * openmm.unit.kelvin,
            1.0 / openmm.unit.picosecond,
            0.002 * openmm.unit.picosecond,
        )

        # Create the context.
        context = openmm.Context(ml_system, integrator)

        # Set the positions.
        context.setPositions(inpcrd.positions)

        # Get the energy.
        state = context.getState(getEnergy=True)
        nrg_openmm = state.getPotentialEnergy().value_in_unit(
            openmm.unit.kilojoules_per_mole
        )

        # Make sure the energies are close.
        assert np.isclose(nrg_openmm, nrg_sire, rtol=1e-3)


@pytest.mark.skipif(not has_emle, reason="emle-engine is not installed")
def test_emle_indirect(ala_mols):
    """
    Make sure that a QM/MM dynamics object can be created using the indirect
    setup for EMLE engines.
    """

    # Create a local copy of the test system.
    mols = ala_mols.clone()

    # Create an EMLE calculator.
    calculator = EMLECalculator(backend="torchani", device="cpu")

    # Create an EMLE engine bound to the calculator.
    emle_mols, engine = sr.qm.create_engine(
        mols, mols[0], calculator, callback="_sire_callback"
    )

    # Create a QM/MM capable dynamics object.
    d = emle_mols.dynamics(
        timestep="1fs",
        constraint="none",
        qm_engine=engine,
        cutoff_type="pme",
        cutoff="7.5 A",
        platform="cpu",
    )

    # Get the potential energy. This will fail if the callback can't be found.
    d.current_potential_energy()


def test_create_engine(ala_mols):
    """
    Make sure that a QM/MM engine can be created and used via a simple callback
    function.
    """

    # A test callback function. Returns a known energy and dummy forces.
    def callback(numbers_qm, charges_mm, xyz_qm, xyz_mm):
        return (42, xyz_qm, xyz_mm)

    # Create a local copy of the test system.
    mols = ala_mols.clone()

    # Create a QM engine bound to the callback.
    qm_mols, engine = sr.qm.create_engine(
        mols,
        mols[0],
        callback,
        callback=None,
    )

    # Create a QM/MM capable dynamics object for the QM molecule only.
    d = qm_mols[0].dynamics(
        timestep="1fs",
        constraint="none",
        qm_engine=engine,
        cutoff_type="pme",
        cutoff="7.5 A",
        platform="cpu",
    )

    # Get the potential energy. This should equal the value returned by the
    # callback, i.e. 42.
    nrg = d.current_potential_energy().to("kJ_per_mol")

    # Make sure the energy is correct.
    assert nrg == 42
