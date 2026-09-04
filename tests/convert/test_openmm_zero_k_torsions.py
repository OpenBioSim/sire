import sire as sr
import pytest


@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="openmm support is not available",
)
def test_openmm_zero_k_torsions(openmm_platform):
    """
    GROMACS topologies can contain torsions with a zero force constant, which
    the GroTop reader drops. In AMBER-derived topologies those torsions are what
    carry the 1-4 scaling, so AmberParams::validateAndFix() has to rebuild a null
    dihedral for every 1-4 pair that no longer has one. Check that every 1-4 pair
    in the file survives that repair, and that repeating the conversion is stable
    (the repair loop is run in parallel and used to race).
    """
    import openmm

    mols = sr.load_test_files("zero_k_torsions.gro", "zero_k_torsions.top")

    # the number of entries in the [ pairs ] section of the topology
    num_pairs = 8575

    for _ in range(5):
        omm = sr.convert.to(mols, "openmm", map={"platform": openmm_platform})

        nonbonded = None

        for force in omm.getSystem().getForces():
            if isinstance(force, openmm.NonbondedForce):
                nonbonded = force
                break

        assert nonbonded is not None

        # count the exceptions that are real 1-4 pairs, i.e. those that are
        # scaled rather than fully excluded
        num_14 = 0

        for i in range(nonbonded.getNumExceptions()):
            _, _, chg, _, lj = nonbonded.getExceptionParameters(i)

            if chg._value != 0.0 or lj._value != 0.0:
                num_14 += 1

        assert num_14 == num_pairs
