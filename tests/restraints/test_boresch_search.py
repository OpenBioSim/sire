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
        assert isinstance(boresch_result, tuple) and len(boresch_result) == 3

    def test_starting_structure_is_system(self, boresch_result):
        """
        The third return value is a single trajectory frame (the least-strained
        starting structure), returned as a sire System with the same molecules
        as the search system.
        """
        _, _, starting_structure = boresch_result
        assert isinstance(starting_structure, sr.system.System)
        assert starting_structure.num_molecules() > 0

    def test_equilibria_match_anchor_geometry(self, boresch_result):
        """
        The restraint's stored equilibrium values must match the geometry of its
        own anchor atoms, measured on the returned starting structure. This guards
        against the anchor atoms being reordered relative to the sampled DOFs: an
        "atomidx a, b, c" selection returns the atoms sorted by index rather than
        in the required [r1, r2, r3] / [l1, l2, l3] order, which would scramble the
        equilibria (r0/theta0/phi0) onto the wrong atoms and badly strain the
        restraint. The starting structure is the least-strained frame, so on it the
        deviation from the equilibria is small; a scrambled restraint deviates by
        up to ~180 degrees.
        """
        restraints, _, start = boresch_result
        r = restraints.at(0)
        rec = list(r.receptor_atoms())
        lig = list(r.ligand_atoms())
        a = start.atoms()
        R1, R2, R3 = a[rec[0]], a[rec[1]], a[rec[2]]
        L1, L2, L3 = a[lig[0]], a[lig[1]], a[lig[2]]

        # Distance.
        assert abs(float(sr.measure(L1, R1).value()) - float(r.r0().value())) < 1.5

        # Angles and dihedrals (wrap the dihedral deviation into [0, 180]).
        measured_angles = [
            float(sr.measure(R2, R1, L1).to(sr.units.degrees)),
            float(sr.measure(R1, L1, L2).to(sr.units.degrees)),
        ]
        target_angles = [float(t.to(sr.units.degrees)) for t in r.theta0()]
        for got, tgt in zip(measured_angles, target_angles):
            assert abs(got - tgt) < 45.0

        measured_dih = [
            float(sr.measure(R3, R2, R1, L1).to(sr.units.degrees)),
            float(sr.measure(R2, R1, L1, L2).to(sr.units.degrees)),
            float(sr.measure(R1, L1, L2, L3).to(sr.units.degrees)),
        ]
        target_dih = [float(p.to(sr.units.degrees)) for p in r.phi0()]
        for got, tgt in zip(measured_dih, target_dih):
            delta = abs(got - tgt) % 360.0
            assert min(delta, 360.0 - delta) < 45.0

    def test_restraints_type(self, boresch_result):
        restraints, _, _ = boresch_result
        assert isinstance(restraints, _SireMM.BoreschRestraints)

    def test_single_restraint(self, boresch_result):
        restraints, _, _ = boresch_result
        assert restraints.at(0) is not None

    def test_angle_potential_defaults_to_restricted_bending(self, boresch_result):
        """
        boresch_search() defaults to "restricted_bending" (unlike
        sire.restraints.boresch's own default of "harmonic"), since it
        directly avoids the collinearity singularity the restraint search
        protocols are designed to steer away from in the first place.
        """
        restraints, _, _ = boresch_result
        assert restraints.angle_potential() == "restricted_bending"

    def test_restraint_lever_matches_protocol(self, boresch_result, request):
        """
        boresch_search() defaults restraint_lever to match the protocol:
        "split" for rxrx (matching the RXRX protocol's staged restraint
        turn-on), "combined" for aldeghi (matching the Aldeghi protocol,
        where all six restraint terms are turned on together).
        """
        restraints, _, _ = boresch_result
        protocol = request.node.callspec.params["boresch_result"]
        expected = "split" if protocol == "rxrx" else "combined"
        assert restraints.restraint_lever() == expected

    def test_correction_is_negative(self, boresch_result):
        """Standard state correction is always negative (costs free energy to restrain)."""
        _, correction, _ = boresch_result
        assert float(correction.to(sr.units.kcal_per_mol)) < 0

    def test_distance_positive(self, boresch_result):
        restraints, _, _ = boresch_result
        r = restraints.at(0)
        assert float(r.r0().value()) > 0

    def test_angles_not_collinear(self, boresch_result):
        """Both anchor angles must be away from 0° and 180° (our stability filter)."""
        restraints, _, _ = boresch_result
        r = restraints.at(0)
        for theta in r.theta0():
            deg = float(theta.to(sr.units.degrees))
            assert 10.0 < deg < 170.0

    def test_force_constants_positive(self, boresch_result):
        restraints, _, _ = boresch_result
        r = restraints.at(0)
        assert float(r.kr().value()) > 0
        for k in r.ktheta():
            assert float(k.value()) > 0
        for k in r.kphi():
            assert float(k.value()) > 0

    def test_receptor_atoms_count(self, boresch_result):
        """Receptor selection must yield exactly 3 anchor atoms."""
        restraints, _, _ = boresch_result
        assert len(list(restraints.at(0).receptor_atoms())) == 3

    def test_ligand_atoms_count(self, boresch_result):
        """Ligand selection must yield exactly 3 anchor atoms."""
        restraints, _, _ = boresch_result
        assert len(list(restraints.at(0).ligand_atoms())) == 3

    @pytest.mark.parametrize("protocol", PROTOCOLS)
    def test_restraint_idx_selects_different_candidate(self, abfe_system, protocol):
        from sire.restraints import boresch_search

        r0, _, _ = boresch_search(abfe_system, protocol=protocol, restraint_idx=0)
        r1, _, _ = boresch_search(abfe_system, protocol=protocol, restraint_idx=1)
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

        restraints, _, _ = boresch_search(
            abfe_system, protocol=protocol, angle_potential="harmonic"
        )
        assert restraints.angle_potential() == "harmonic"

    @pytest.mark.parametrize("protocol", PROTOCOLS)
    def test_angle_potential_invalid_raises(self, abfe_system, protocol):
        from sire.restraints import boresch_search

        with pytest.raises(ValueError, match="angle_potential"):
            boresch_search(abfe_system, protocol=protocol, angle_potential="nonsense")

    @pytest.mark.parametrize("protocol", PROTOCOLS)
    def test_restraint_lever_override(self, abfe_system, protocol):
        from sire.restraints import boresch_search

        # Explicitly pick whichever value isn't the protocol's own default,
        # to confirm the override actually takes effect.
        override = "combined" if protocol == "rxrx" else "split"
        restraints, _, _ = boresch_search(
            abfe_system, protocol=protocol, restraint_lever=override
        )
        assert restraints.restraint_lever() == override

    @pytest.mark.parametrize("protocol", PROTOCOLS)
    def test_restraint_lever_invalid_raises(self, abfe_system, protocol):
        from sire.restraints import boresch_search

        with pytest.raises(ValueError, match="restraint_lever"):
            boresch_search(abfe_system, protocol=protocol, restraint_lever="nonsense")

    def test_force_constant_override(self, abfe_system):
        """'force_constant' is an Aldeghi-only kwarg (RXRX uses two separate
        force constants, matching the fixed protocol values from the paper)."""
        from sire.restraints import boresch_search

        kval = 10.0
        restraints, _, _ = boresch_search(
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
