"""
Unit tests for boresch_search().

Test data (boresch_restraints.prm7 + boresch_restraints.dcd) is a lambda=0
trajectory of the 1jr5 protein with decoupled ligand01.
"""

import numpy as np
import pytest
import sire as sr
import sire.legacy.MM as _SireMM

PROTOCOLS = ["rxrx", "aldeghi"]


@pytest.fixture(scope="module")
def abfe_system():
    """
    Load the protein-ligand test system with embedded trajectory frames.
    The topology is loaded from prm7; the trajectory from DCD. Molecule 0
    is the protein, molecule 1 is the ligand, the rest is water and ions.
    The ligand is then decoupled so that is_perturbable is set, matching
    what SOMD2's runner provides to boresch_search().
    """
    mols = sr.load_test_files("boresch_restraints.prm7", "boresch_restraints.dcd")
    lig_decoupled = sr.morph.decouple(mols.molecule(1), as_new_molecule=False)
    mols.update(lig_decoupled)
    return mols


@pytest.fixture(scope="module", params=PROTOCOLS)
def boresch_result(request, abfe_system):
    from sire.restraints import boresch_search

    return boresch_search(abfe_system, protocol=request.param, temperature="298 K")


class TestGenerateBoreschRestraint:
    def test_returns_tuple(self, boresch_result):
        assert isinstance(boresch_result, tuple) and len(boresch_result) == 2

    def test_restraints_type(self, boresch_result):
        restraints, _ = boresch_result
        assert isinstance(restraints, _SireMM.BoreschRestraints)

    def test_single_restraint(self, boresch_result):
        restraints, _ = boresch_result
        assert restraints.at(0) is not None

    def test_angle_potential_defaults_to_restricted_bending(self, boresch_result):
        """
        boresch_search() defaults to "restricted_bending" (unlike
        sire.restraints.boresch's own default of "harmonic"), since it
        directly avoids the collinearity singularity the restraint search
        protocols are designed to steer away from in the first place.
        """
        restraints, _ = boresch_result
        assert restraints.angle_potential() == "restricted_bending"

    def test_correction_is_negative(self, boresch_result):
        """Standard state correction is always negative (costs free energy to restrain)."""
        _, correction = boresch_result
        assert float(correction.to(sr.units.kcal_per_mol)) < 0

    def test_distance_positive(self, boresch_result):
        restraints, _ = boresch_result
        r = restraints.at(0)
        assert float(r.r0().value()) > 0

    def test_angles_not_collinear(self, boresch_result):
        """Both anchor angles must be away from 0° and 180° (our stability filter)."""
        restraints, _ = boresch_result
        r = restraints.at(0)
        for theta in r.theta0():
            deg = float(theta.to(sr.units.degrees))
            assert 10.0 < deg < 170.0

    def test_force_constants_positive(self, boresch_result):
        restraints, _ = boresch_result
        r = restraints.at(0)
        assert float(r.kr().value()) > 0
        for k in r.ktheta():
            assert float(k.value()) > 0
        for k in r.kphi():
            assert float(k.value()) > 0

    def test_receptor_atoms_count(self, boresch_result):
        """Receptor selection must yield exactly 3 anchor atoms."""
        restraints, _ = boresch_result
        assert len(list(restraints.at(0).receptor_atoms())) == 3

    def test_ligand_atoms_count(self, boresch_result):
        """Ligand selection must yield exactly 3 anchor atoms."""
        restraints, _ = boresch_result
        assert len(list(restraints.at(0).ligand_atoms())) == 3

    @pytest.mark.parametrize("protocol", PROTOCOLS)
    def test_restraint_idx_selects_different_candidate(self, abfe_system, protocol):
        from sire.restraints import boresch_search

        r0, _ = boresch_search(abfe_system, protocol=protocol, restraint_idx=0)
        r1, _ = boresch_search(abfe_system, protocol=protocol, restraint_idx=1)
        # Different candidates must differ in the receptor and/or ligand anchor
        # atoms: two candidates can share the same receptor anchor triplet
        # while using a different ligand triplet branching off the same l1,
        # or vice versa.
        r0_sextuplet = (
            list(r0.at(0).receptor_atoms()),
            list(r0.at(0).ligand_atoms()),
        )
        r1_sextuplet = (
            list(r1.at(0).receptor_atoms()),
            list(r1.at(0).ligand_atoms()),
        )
        assert r0_sextuplet != r1_sextuplet

    @pytest.mark.parametrize("protocol", PROTOCOLS)
    def test_too_few_frames_raises(self, abfe_system, protocol):
        from sire.restraints import boresch_search

        with pytest.raises(ValueError, match="frame"):
            boresch_search(abfe_system, protocol=protocol, min_frames=10_000)

    @pytest.mark.parametrize("protocol", PROTOCOLS)
    def test_angle_potential_harmonic_override(self, abfe_system, protocol):
        from sire.restraints import boresch_search

        restraints, _ = boresch_search(
            abfe_system, protocol=protocol, angle_potential="harmonic"
        )
        assert restraints.angle_potential() == "harmonic"

    @pytest.mark.parametrize("protocol", PROTOCOLS)
    def test_angle_potential_invalid_raises(self, abfe_system, protocol):
        from sire.restraints import boresch_search

        with pytest.raises(ValueError, match="angle_potential"):
            boresch_search(abfe_system, protocol=protocol, angle_potential="nonsense")

    def test_force_constant_override(self, abfe_system):
        """'force_constant' is an Aldeghi-only kwarg (RXRX uses two separate
        force constants, matching the fixed protocol values from the paper)."""
        from sire.restraints import boresch_search

        kval = 10.0
        restraints, _ = boresch_search(
            abfe_system,
            protocol="aldeghi",
            force_constant=f"{kval} kcal mol-1 A-2",
        )
        r = restraints.at(0)
        assert np.isclose(float(r.kr().value()), kval, atol=1e-3)

    def test_tight_cutoff_raises(self, abfe_system):
        """'cutoff' is an Aldeghi-only kwarg (RXRX has no equivalent distance
        cutoff; candidates are seeded from hydrogen-bond occupancy instead)."""
        from sire.restraints import boresch_search

        with pytest.raises(ValueError, match="cutoff"):
            boresch_search(abfe_system, protocol="aldeghi", cutoff="0.1 A")
