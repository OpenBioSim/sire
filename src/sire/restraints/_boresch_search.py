"""
Automatic Boresch restraint generation for ABFE simulations.
"""

__all__ = ["boresch_search"]

from collections import deque as _deque

import numpy as _np

import sire as _sr

# ---------------------------------------------------------------------------
# Shared helpers, used by both the "rxrx" and "aldeghi" search protocols.
# ---------------------------------------------------------------------------


def _nonH_bonded(connectivity, at_idx, mol, is_lig, ghost_elem, h_elem):
    """Non-H atoms bonded to at_idx (lambda=0 non-ghost for ligand atoms)."""
    bonded = connectivity.connections_to(at_idx)
    result = []
    for b_idx in bonded:
        if is_lig:
            elem = mol.atom(b_idx).property("element0")
            if elem == ghost_elem or elem == h_elem:
                continue
        else:
            if mol.atom(b_idx).property("element") == h_elem:
                continue
        result.append(b_idx)
    return result


def _build_triplets(connectivity, anchor_idx, mol, is_lig, ghost_elem, h_elem):
    """
    All bonded chains anchor-a2-a3 of non-H (non-ghost) atoms.

    A given anchor can have several bonded heavy-atom neighbours (e.g. ring
    branch points), and picking an arbitrary one can point a2/a3 back into
    the parent molecule, making the corresponding Boresch angle collinear
    with the restraint distance vector. All chains are returned so that the
    caller's angle filter/scoring can reject bad geometries in favour of a
    better one, rather than committing to a single (possibly degenerate)
    choice up front.
    """
    triplets = []
    for a2 in _nonH_bonded(connectivity, anchor_idx, mol, is_lig, ghost_elem, h_elem):
        for a3 in _nonH_bonded(connectivity, a2, mol, is_lig, ghost_elem, h_elem):
            if a3 != anchor_idx:
                triplets.append((anchor_idx, a2, a3))
    return triplets


def _circular_mean_std(values):
    """
    Circular mean/std (radians) for periodic dihedral samples.

    Dihedral angles wrap at +-pi, so a plain arithmetic mean/std is wrong
    whenever the true distribution straddles the wrap point (e.g. a mean
    dihedral near +-180 degrees): samples flip between +179 and -179
    degrees, giving a mean near 0 and a hugely inflated std. This uses the
    standard circular statistics (mean of unit vectors, then wrapped
    deviation) instead.
    """
    mean = float(_np.arctan2(_np.sin(values).mean(), _np.cos(values).mean()))
    dtheta = _np.abs(values - mean)
    dtheta = _np.minimum(dtheta, 2.0 * _np.pi - dtheta)
    std = float(_np.sqrt(_np.mean(dtheta**2)))
    return mean, std


def _sample_dofs(system, pert_mol_num, sextuplets, n_frames):
    """
    Second trajectory pass: sample the six Boresch DOFs for every candidate
    sextuplet. Returns arrays of shape (n_sextuplets, n_frames).
    """
    n_sext = len(sextuplets)

    dof_r = _np.zeros((n_sext, n_frames))
    dof_tA = _np.zeros((n_sext, n_frames))  # thetaA: r2-r1-l1
    dof_tB = _np.zeros((n_sext, n_frames))  # thetaB: r1-l1-l2
    dof_pA = _np.zeros((n_sext, n_frames))  # phiA: r3-r2-r1-l1
    dof_pB = _np.zeros((n_sext, n_frames))  # phiB: r2-r1-l1-l2
    dof_pC = _np.zeros((n_sext, n_frames))  # phiC: r1-l1-l2-l3

    for f_idx, frame in enumerate(system.trajectory()):
        space = frame.space()
        lig_frame = frame.molecule(pert_mol_num)
        rec_frame_cache = {}

        for s_idx, s in enumerate(sextuplets):
            r1_mol_num = s["r1_mol_num"]
            if r1_mol_num not in rec_frame_cache:
                rec_frame_cache[r1_mol_num] = frame.molecule(r1_mol_num)
            rec_frame = rec_frame_cache[r1_mol_num]

            l1 = lig_frame.atom(s["l1_idx"])
            l2 = lig_frame.atom(s["l2_idx"])
            l3 = lig_frame.atom(s["l3_idx"])
            r1 = rec_frame.atom(s["r1_idx"])
            r2 = rec_frame.atom(s["r2_idx"])
            r3 = rec_frame.atom(s["r3_idx"])

            dof_r[s_idx, f_idx] = float(
                space.calc_dist(l1.coordinates(), r1.coordinates())
            )
            dof_tA[s_idx, f_idx] = float(_sr.measure(r2, r1, l1).value())
            dof_tB[s_idx, f_idx] = float(_sr.measure(r1, l1, l2).value())
            dof_pA[s_idx, f_idx] = float(_sr.measure(r3, r2, r1, l1).value())
            dof_pB[s_idx, f_idx] = float(_sr.measure(r2, r1, l1, l2).value())
            dof_pC[s_idx, f_idx] = float(_sr.measure(r1, l1, l2, l3).value())

    return dof_r, dof_tA, dof_tB, dof_pA, dof_pB, dof_pC


def _best_starting_frame(
    system,
    best_frame_idx,
):
    """
    Extract a single trajectory frame as a standalone starting structure.

    The equilibrium values of the generated restraint are trajectory
    *averages*, so the structure the search was seeded from (e.g. the output
    of a separate equilibration pipeline) is generally not consistent with
    them - restarting production from it can leave the restraint badly
    strained at t=0 and blow the simulation up as the restraint is switched
    on. Returning the frame closest to the restraint equilibrium (see
    ``_least_strained_frame``) lets the caller seed production from a
    structure at which the restraint is essentially relaxed.
    """
    return system.trajectory()[best_frame_idx].current()


def _least_strained_frame(
    dof_r_b,
    dof_tA_b,
    dof_tB_b,
    dof_pA_b,
    dof_pB_b,
    dof_pC_b,
    r0,
    tA0,
    tB0,
    pA0,
    pB0,
    pC0,
    kr,
    ktheta_A,
    ktheta_B,
    kphi_A,
    kphi_B,
    kphi_C,
):
    """
    Index of the trajectory frame at which the generated restraint is least
    strained, i.e. the frame whose six Boresch DOFs give the lowest total
    restraint bias energy against the chosen equilibrium values (and force
    constants). Dihedral deviations are wrapped into [-pi, pi]. All arguments
    are per-frame arrays (the ``[best]`` slices) or scalars for the chosen
    sextuplet; angles/dihedrals are in radians and r in Angstrom, matching
    the units the equilibrium values and force constants are expressed in, so
    the weighted sum is a genuine (kcal mol-1) restraint energy.
    """

    def _wrap(d):
        return (d + _np.pi) % (2.0 * _np.pi) - _np.pi

    strain = (
        kr * (dof_r_b - r0) ** 2
        + ktheta_A * (dof_tA_b - tA0) ** 2
        + ktheta_B * (dof_tB_b - tB0) ** 2
        + kphi_A * _wrap(dof_pA_b - pA0) ** 2
        + kphi_B * _wrap(dof_pB_b - pB0) ** 2
        + kphi_C * _wrap(dof_pC_b - pC0) ** 2
    )
    return int(_np.argmin(strain))


def _assemble_restraints(
    system,
    pert_mol_num,
    s,
    r0,
    tA0,
    tB0,
    pA0,
    pB0,
    pC0,
    kr_str,
    kt_str,
    kp_str,
    temperature,
    angle_potential,
    restraint_lever,
):
    """Build the sire.mm.BoreschRestraints object and standard state correction."""
    from ._restraints import boresch as _boresch
    from ._standard_state_correction import get_standard_state_correction as _get_ssc

    # Select the anchor atoms in the exact [r1, r2, r3] / [l1, l2, l3] order.
    # NB: an "atomidx a, b, c" search string returns the atoms sorted by index,
    # NOT in listed order, so it cannot be used here - it would scramble the
    # anchor ordering relative to the sampled DOFs (r0/theta0/phi0), producing a
    # restraint whose equilibrium values do not match its own anchor geometry.
    # List-indexing a molecule's atoms preserves the given order.
    receptor_mol_atoms = system.molecule(s["r1_mol_num"]).atoms()
    receptor_atoms = receptor_mol_atoms[
        [s["r1_idx"].value(), s["r2_idx"].value(), s["r3_idx"].value()]
    ]
    ligand_mol_atoms = system.molecule(pert_mol_num).atoms()
    ligand_atoms = ligand_mol_atoms[
        [s["l1_idx"].value(), s["l2_idx"].value(), s["l3_idx"].value()]
    ]

    restraints = _boresch(
        system,
        receptor=receptor_atoms,
        ligand=ligand_atoms,
        kr=kr_str,
        ktheta=kt_str,
        kphi=kp_str,
        r0=f"{r0:.6f} A",
        theta0=[f"{tA0:.6f} rad", f"{tB0:.6f} rad"],
        phi0=[f"{pA0:.6f} rad", f"{pB0:.6f} rad", f"{pC0:.6f} rad"],
        temperature=temperature,
        angle_potential=angle_potential,
        restraint_lever=restraint_lever,
    )

    correction = _get_ssc(restraints[0], temperature)

    return restraints, correction


def _parse_unit(value, name, unit, unit_desc):
    """
    Parse 'value' (a unit str or a GeneralUnit) into a GeneralUnit with the
    dimensions of 'unit', raising a clear error otherwise. Every str-or-
    GeneralUnit parameter in this module (temperature, distance cutoffs,
    angle cutoffs, force constants, ...) is validated the same way, so that
    they can all be supplied consistently as unit strings, e.g. "298 K",
    "3.5 A", "120 degrees", "1 kcal mol-1 A-2".
    """
    if isinstance(value, str):
        try:
            value = _sr.u(value)
        except Exception as e:
            raise ValueError(f"Could not parse '{name}' as a unit: {e}") from e
    else:
        try:
            float(value.to(unit))
        except Exception as e:
            raise TypeError(
                f"'{name}' must be a str or a GeneralUnit with {unit_desc} "
                f"dimensions, got {type(value)}: {e}"
            ) from e
    return value


def _validate_angle_potential(angle_potential):
    """
    Validate 'angle_potential', shared by both search protocols. Must be
    either "harmonic" or "restricted_bending" (see
    sire.restraints.boresch's angle_potential parameter). Unlike
    sire.restraints.boresch's own default of "harmonic" (kept for backwards
    compatibility with existing callers), boresch_search defaults to
    "restricted_bending" for the restraints it generates itself, since the
    collinearity singularity this avoids is exactly the failure mode the
    restraint search protocols are designed to steer away from in the first
    place.
    """
    if angle_potential not in ("harmonic", "restricted_bending"):
        raise ValueError(
            "'angle_potential' must be either 'harmonic' or 'restricted_bending', "
            f"got {angle_potential!r}"
        )


def _validate_restraint_lever(restraint_lever):
    """
    Validate 'restraint_lever', shared by both search protocols. Must be
    either "combined" or "split" (see sire.restraints.boresch's
    restraint_lever parameter).
    """
    if restraint_lever not in ("combined", "split"):
        raise ValueError(
            "'restraint_lever' must be either 'combined' or 'split', "
            f"got {restraint_lever!r}"
        )


# ---------------------------------------------------------------------------
# RXRX protocol (default): H-bond-driven restraint search algorithm named
# and described in:
# "Optimizing Absolute Binding Free Energy Calculations for Production Usage"
# https://doi.org/10.1021/acs.jctc.5c00861
#
# Atoms from protein-ligand hydrogen bonds rather than from bulk distance
# variance, and scores candidates with a formula that explicitly penalises
# theta angles near 0/180 degrees, in order to avoid the Boresch angle
# singularities that plague automated restraint search.
# ---------------------------------------------------------------------------


def _boresch_search_rxrx(
    system,
    temperature="298 K",
    protein_selection="not water",
    hbond_distance_cutoff="3.5 A",
    hbond_angle_cutoff="120 degrees",
    occupancy_cutoff=0.5,
    restraint_idx=0,
    force_constant_r=None,
    force_constant_angle=None,
    min_frames=50,
    angle_potential="restricted_bending",
    restraint_lever="split",
):
    """
    Generate a Boresch restraint using the RXRX restraint search algorithm.

    The system must contain trajectory frames (i.e. be the result of
    ``dynamics.commit()`` after a run with ``frame_frequency > 0``). Ligand
    anchor atoms are restricted to non-terminal heavy atoms, and candidate
    protein-ligand atom sextuplets are seeded from hydrogen bonds with
    occupancy above ``occupancy_cutoff``, with the protein anchor mapped to
    the Cα of the hydrogen-bonded residue. Among candidates that survive an
    angle-range filter, the ligand anchor closest to the ligand's centre of
    mass is chosen, and the lowest-scoring sextuplet using that anchor
    (Equation 1 of the RXRX paper) is returned.

    Parameters
    ----------

    system : sire.system.System
        A Sire system with embedded trajectory frames. Must contain exactly
        one perturbable molecule.

    temperature : str or GeneralUnit, optional
        Simulation temperature. Defaults to 298 K.

    protein_selection : str
        Sire selection string used to identify protein atoms considered as
        hydrogen-bond partners. Defaults to all non-water atoms.

    hbond_distance_cutoff : str or GeneralUnit
        Maximum donor-acceptor heavy atom distance for a hydrogen bond.

    hbond_angle_cutoff : str or GeneralUnit
        Minimum donor-H...acceptor angle for a hydrogen bond.

    occupancy_cutoff : float
        Minimum fraction of frames (0-1) in which a hydrogen bond must be
        present for it to be considered as a restraint anchor.

    restraint_idx : int
        Index into the candidate list (sharing the chosen ligand anchor
        atom) sorted by ascending RXRX score (0 = best).

    force_constant_r : str or GeneralUnit, optional
        Distance restraint force constant. Defaults to the RXRX protocol
        value of 1 kcal mol-1 A-2.

    force_constant_angle : str or GeneralUnit, optional
        Angle and dihedral restraint force constant. Defaults to the RXRX
        protocol value of 80 kcal mol-1 rad-2.

    min_frames : int
        Minimum number of trajectory frames required.

    angle_potential : str
        The functional form used for the two Boresch angle restraint terms,
        either "harmonic" or "restricted_bending" (see
        sire.restraints.boresch). Defaults to "restricted_bending", which
        diverges as the angle approaches 0 or 180 degrees, preventing it
        from ever reaching the Boresch collinearity singularity - the exact
        failure mode this restraint search is designed to avoid.

    restraint_lever : str
        How the restraint's six degrees of freedom are grouped into
        lambda-addressable OpenMM Forces (see sire.restraints.boresch).
        Defaults to "split", matching the RXRX protocol's staged restraint
        turn-on, where the dihedral terms and the distance/angle terms are
        turned on according to different lambda schedule equations.

    Returns
    -------

    restraints : sire.mm.BoreschRestraints
        The generated Boresch restraints object, ready to pass to
        ``Config.restraints``.

    correction : sire.units.GeneralUnit
        Standard state correction in kcal mol-1.

    starting_structure : sire.system.System
        The trajectory frame at which the generated restraint is least
        strained (its six Boresch degrees of freedom are closest to the
        restraint equilibrium values). Because the equilibrium values are
        trajectory averages, the structure the search was seeded from is
        generally not consistent with them; seeding production from this
        frame instead avoids a large restraint force at t=0 that can
        otherwise destabilise the simulation as the restraint is switched on.
    """
    from ..legacy import Mol as _SireMol
    from ..system import System as _System

    if not isinstance(system, _System):
        raise TypeError(
            f"'system' must be of type 'sire.system.System', got {type(system)}"
        )

    temperature = _parse_unit(
        temperature, "temperature", _sr.units.kelvin, "temperature"
    )

    if not isinstance(protein_selection, str):
        raise TypeError(
            f"'protein_selection' must be a str, got {type(protein_selection)}"
        )

    hbond_distance_cutoff = _parse_unit(
        hbond_distance_cutoff, "hbond_distance_cutoff", _sr.units.angstrom, "length"
    )
    hbond_cutoff_A = float(hbond_distance_cutoff.to(_sr.units.angstrom))

    hbond_angle_cutoff = _parse_unit(
        hbond_angle_cutoff, "hbond_angle_cutoff", _sr.units.degrees, "angle"
    )
    hbond_angle_cutoff_deg = float(hbond_angle_cutoff.to(_sr.units.degrees))

    if not (0.0 < occupancy_cutoff <= 1.0):
        raise ValueError(
            f"'occupancy_cutoff' must be in (0, 1], got {occupancy_cutoff!r}"
        )

    if not isinstance(restraint_idx, int) or restraint_idx < 0:
        raise ValueError(
            f"'restraint_idx' must be a non-negative int, got {restraint_idx!r}"
        )

    r_fc_unit = _sr.units.kcal_per_mol / _sr.units.angstrom**2
    angle_fc_unit = _sr.units.kcal_per_mol / _sr.units.radians**2

    if force_constant_r is not None:
        force_constant_r = _parse_unit(
            force_constant_r, "force_constant_r", r_fc_unit, "energy length-2"
        )
    if force_constant_angle is not None:
        force_constant_angle = _parse_unit(
            force_constant_angle,
            "force_constant_angle",
            angle_fc_unit,
            "energy angle-2",
        )

    if not isinstance(min_frames, int) or min_frames < 1:
        raise ValueError(f"'min_frames' must be a positive int, got {min_frames!r}")

    _validate_angle_potential(angle_potential)
    _validate_restraint_lever(restraint_lever)

    n_frames = system.num_frames()
    if n_frames < min_frames:
        raise ValueError(
            f"Trajectory has only {n_frames} frame(s); at least {min_frames} are required "
            "for reliable Boresch restraint generation."
        )

    ghost_elem = _SireMol.Element(0)
    h_elem = _SireMol.Element("H")
    n_elem = _SireMol.Element("N")
    o_elem = _SireMol.Element("O")

    # -------------------------------------------------------------------------
    # 1. Locate the perturbable molecule (ligand) and its candidate atoms.
    #
    # The full heavy-atom set is used as the hydrogen-bond donor/acceptor
    # pool; the candidate anchor set is restricted to non-terminal atoms
    # (bonded to >= 2 non-H heavy atoms), matching the RXRX paper's Figure 2B.
    # -------------------------------------------------------------------------
    pert_mols = system.molecules("property is_perturbable")
    if pert_mols.num_molecules() != 1:
        raise ValueError(
            "System must contain exactly one perturbable molecule for Boresch "
            f"restraint generation; found {pert_mols.num_molecules()}."
        )
    pert_mol = pert_mols.molecule(0)
    pert_mol_num = pert_mol.number()
    lig_connectivity = pert_mol.connectivity()

    lig_heavy_idxs = []
    for atom in pert_mol.atoms():
        elem0 = atom.property("element0")
        if elem0 != ghost_elem and elem0 != h_elem:
            lig_heavy_idxs.append(atom.index())

    if len(lig_heavy_idxs) < 3:
        raise ValueError(
            f"Ligand has only {len(lig_heavy_idxs)} non-ghost, non-hydrogen atom(s); "
            "need at least 3 for Boresch restraints."
        )

    # Lambda=0 atomic masses, used for the mass-weighted centre of mass below
    # (topology doesn't change between frames, so this is computed once).
    lig_masses = _np.array(
        [float(pert_mol.atom(idx).property("mass0").value()) for idx in lig_heavy_idxs]
    )

    candidate_set = set()
    for idx in lig_heavy_idxs:
        n_heavy = len(
            _nonH_bonded(lig_connectivity, idx, pert_mol, True, ghost_elem, h_elem)
        )
        if n_heavy >= 2:
            candidate_set.add(idx)

    if not candidate_set:
        raise ValueError(
            "No non-terminal ligand heavy atoms (bonded to >= 2 heavy atoms) "
            "were found for RXRX restraint search."
        )

    def _attached_H(connectivity, at_idx, mol, is_lig):
        bonded = connectivity.connections_to(at_idx)
        result = []
        for b_idx in bonded:
            elem = mol.atom(b_idx).property("element0" if is_lig else "element")
            if elem == h_elem:
                result.append(b_idx)
        return result

    lig_hbond_idxs = [
        idx
        for idx in lig_heavy_idxs
        if pert_mol.atom(idx).property("element0") in (n_elem, o_elem)
    ]
    if not lig_hbond_idxs:
        raise ValueError(
            "Ligand has no N/O heavy atoms to use as hydrogen-bond partners."
        )
    lig_donor_H = {
        idx: _attached_H(lig_connectivity, idx, pert_mol, True)
        for idx in lig_hbond_idxs
    }

    # -------------------------------------------------------------------------
    # 2. Locate protein hydrogen-bond partner atoms (N/O) and cache each
    #    residue's C-alpha atom index for the anchor mapping.
    #
    # 'protein_selection' defaults to the whole protein, which can be
    # thousands of atoms; a hydrogen bond to the ligand can only involve
    # atoms within a few Angstrom of it, so the N/O search is narrowed with
    # a generous spatial "within" clause up front, using Sire's own (native,
    # spatially-indexed) selection engine rather than a manual per-atom
    # Python distance loop. This is evaluated once, using the system's
    # default (first-frame) coordinates.
    # -------------------------------------------------------------------------
    _prefilter_margin_A = 15.0
    try:
        prot_hbond_sel = system[
            f"(element N,O) and ({protein_selection}) and not (molnum {pert_mol_num.value()}) and "
            f"(atoms within {hbond_cutoff_A + _prefilter_margin_A} of (molnum {pert_mol_num.value()}))"
        ]
    except Exception as e:
        raise ValueError(
            f"Could not apply protein selection '{protein_selection}': {e}"
        ) from e

    prot_mol_cache = {}
    prot_hbond_atoms = []  # dicts: mol_num, idx, res_key, attached_H
    res_objs = {}  # res_key -> Residue view, for the lazy Cα lookup below
    for atom in prot_hbond_sel.atoms():
        mn = atom.molecule().number()
        if mn not in prot_mol_cache:
            prot_mol_cache[mn] = system.molecule(mn)
        mol = prot_mol_cache[mn]
        residue = atom.residue()
        res_key = (mn, residue.number())
        res_objs.setdefault(res_key, residue)
        prot_hbond_atoms.append(
            {
                "mol_num": mn,
                "idx": atom.index(),
                "res_key": res_key,
                "attached_H": _attached_H(mol.connectivity(), atom.index(), mol, False),
            }
        )

    if not prot_hbond_atoms:
        raise ValueError(
            f"No protein N/O atoms found within {hbond_cutoff_A + _prefilter_margin_A:.1f} A "
            f"of the ligand for protein selection '{protein_selection}'. Try adjusting "
            "'protein_selection' or the hydrogen-bond cutoffs."
        )

    # Cα lookups are only needed for the (typically small) set of residues
    # that turn out to be hydrogen-bonded below, so resolve them lazily by
    # scanning just that residue's own (small) atom list, rather than every
    # atom of every protein molecule up front.
    res_ca_idx = {}

    def _get_ca_idx(res_key):
        if res_key not in res_ca_idx:
            res_ca_idx[res_key] = None
            for a in res_objs[res_key].atoms():
                if a.name().value() == "CA":
                    res_ca_idx[res_key] = a.index()
                    break
        return res_ca_idx[res_key]

    # -------------------------------------------------------------------------
    # 3. First trajectory pass: hydrogen-bond occupancy per (ligand atom,
    #    protein residue) pair, plus per-frame ligand centre of mass and
    #    candidate anchor positions (for the CoM tie-break in step 5).
    # -------------------------------------------------------------------------
    def _vec3(v):
        return _np.array(
            [
                float(v.x().to(_sr.units.angstrom)),
                float(v.y().to(_sr.units.angstrom)),
                float(v.z().to(_sr.units.angstrom)),
            ]
        )

    def _min_image_vec3(space, from_coord, to_coord):
        """
        Minimum-image displacement vector (Angstrom, numpy) from from_coord
        to to_coord. Sire wraps at the atom level (not per-molecule), so even
        directly-bonded atoms cannot be assumed to share the same periodic
        image; raw coordinate subtraction is not safe. Mirrors the approach
        used in BioSimSpace's SireWrapper._getCenterOfMass.
        """
        from ..legacy.Maths import Vector as _SireVector

        return _vec3(_SireVector(space.calc_dist_vector(from_coord, to_coord)))

    def _angle_deg(space, d_coord, h_coord, a_coord):
        v1 = _min_image_vec3(space, h_coord, d_coord)
        v2 = _min_image_vec3(space, h_coord, a_coord)
        cosang = _np.dot(v1, v2) / (_np.linalg.norm(v1) * _np.linalg.norm(v2))
        return _np.degrees(_np.arccos(_np.clip(cosang, -1.0, 1.0)))

    occ_counts = {}  # (lig_idx, res_key) -> frame count with a satisfied H-bond
    candidate_com_dists = {idx: [] for idx in candidate_set}

    for frame in system.trajectory():
        space = frame.space()
        lig_frame = frame.molecule(pert_mol_num)

        lig_coords = {idx: lig_frame.atom(idx).coordinates() for idx in lig_heavy_idxs}
        ref_coord = lig_coords[lig_heavy_idxs[0]]
        ref_vec3 = _vec3(ref_coord)
        # Unwrap every ligand atom relative to a single reference atom before
        # averaging, so the mass-weighted centre of mass is computed from a
        # mutually consistent (non-split-across-the-box) set of positions.
        lig_unwrapped = {
            idx: ref_vec3 + _min_image_vec3(space, ref_coord, lig_coords[idx])
            for idx in lig_heavy_idxs
        }
        lig_positions = _np.array([lig_unwrapped[idx] for idx in lig_heavy_idxs])
        com = _np.average(lig_positions, axis=0, weights=lig_masses)
        for idx in candidate_set:
            candidate_com_dists[idx].append(
                float(_np.linalg.norm(lig_unwrapped[idx] - com))
            )

        prot_frame_cache = {}

        def _prot_frame(mn):
            if mn not in prot_frame_cache:
                prot_frame_cache[mn] = frame.molecule(mn)
            return prot_frame_cache[mn]

        prot_coords = [
            _prot_frame(p["mol_num"]).atom(p["idx"]).coordinates()
            for p in prot_hbond_atoms
        ]

        pairs_this_frame = set()

        for li in lig_hbond_idxs:
            l_coord = lig_coords[li]
            l_H_coords = [lig_frame.atom(h).coordinates() for h in lig_donor_H[li]]

            for p_idx, p in enumerate(prot_hbond_atoms):
                p_coord = prot_coords[p_idx]

                dist = float(space.calc_dist(l_coord, p_coord))
                if dist > hbond_cutoff_A:
                    continue

                bonded = False
                # Ligand donor -> protein acceptor.
                for h_coord in l_H_coords:
                    if (
                        _angle_deg(space, l_coord, h_coord, p_coord)
                        >= hbond_angle_cutoff_deg
                    ):
                        bonded = True
                        break
                if not bonded:
                    # Protein donor -> ligand acceptor.
                    p_frame = _prot_frame(p["mol_num"])
                    for h_idx in p["attached_H"]:
                        h_coord = p_frame.atom(h_idx).coordinates()
                        if (
                            _angle_deg(space, p_coord, h_coord, l_coord)
                            >= hbond_angle_cutoff_deg
                        ):
                            bonded = True
                            break

                if bonded:
                    pairs_this_frame.add((li, p["res_key"]))

        for pair in pairs_this_frame:
            occ_counts[pair] = occ_counts.get(pair, 0) + 1

    occupied_pairs = [
        pair
        for pair, count in occ_counts.items()
        if count / n_frames > occupancy_cutoff
    ]

    if not occupied_pairs:
        raise ValueError(
            "No protein-ligand hydrogen bonds exceeded the occupancy cutoff "
            f"({occupancy_cutoff}). Try lowering 'occupancy_cutoff' or widening "
            "'hbond_distance_cutoff'/'hbond_angle_cutoff'."
        )

    # -------------------------------------------------------------------------
    # 4. Map each occupied (ligand atom, protein residue) pair to a candidate
    #    six-atom combination: ligand atom -> nearest candidate anchor (BFS
    #    over connectivity), protein residue -> C-alpha, then extend both
    #    into bonded triplets.
    # -------------------------------------------------------------------------
    def _nearest_candidate(atom_idx, connectivity):
        if atom_idx in candidate_set:
            return atom_idx
        visited = {atom_idx}
        queue = _deque([atom_idx])
        while queue:
            cur = queue.popleft()
            for nb in connectivity.connections_to(cur):
                if nb in visited:
                    continue
                if nb in candidate_set:
                    return nb
                visited.add(nb)
                queue.append(nb)
        return None

    seen_sextuplets = set()
    valid_sextuplets = []

    for li, res_key in occupied_pairs:
        l1_idx = _nearest_candidate(li, lig_connectivity)
        if l1_idx is None:
            continue

        r1_mol_num = res_key[0]
        r1_idx = _get_ca_idx(res_key)
        if r1_idx is None:
            continue
        r1_mol = prot_mol_cache[r1_mol_num]
        r1_conn = r1_mol.connectivity()

        ligand_triplets = _build_triplets(
            lig_connectivity, l1_idx, pert_mol, True, ghost_elem, h_elem
        )
        receptor_triplets = _build_triplets(
            r1_conn, r1_idx, r1_mol, False, ghost_elem, h_elem
        )

        for l1_idx, l2_idx, l3_idx in ligand_triplets:
            for r1_idx, r2_idx, r3_idx in receptor_triplets:
                key = (l1_idx, l2_idx, l3_idx, r1_mol_num, r1_idx, r2_idx, r3_idx)
                if key in seen_sextuplets:
                    continue
                seen_sextuplets.add(key)

                valid_sextuplets.append(
                    {
                        "l1_idx": l1_idx,
                        "l2_idx": l2_idx,
                        "l3_idx": l3_idx,
                        "r1_mol_num": r1_mol_num,
                        "r1_idx": r1_idx,
                        "r2_idx": r2_idx,
                        "r3_idx": r3_idx,
                    }
                )

    if not valid_sextuplets:
        raise ValueError(
            "Could not construct a valid Boresch sextuplet from any hydrogen-bonded "
            "candidate. Try adjusting 'protein_selection' or the hydrogen-bond cutoffs."
        )

    # -------------------------------------------------------------------------
    # 5. Second trajectory pass: sample all 6 DOFs for every valid sextuplet.
    # -------------------------------------------------------------------------
    dof_r, dof_tA, dof_tB, dof_pA, dof_pB, dof_pC = _sample_dofs(
        system, pert_mol_num, valid_sextuplets, n_frames
    )

    # -------------------------------------------------------------------------
    # 6. Filter to the RXRX angle window (45-135 degrees), then reduce to the
    #    single ligand anchor atom closest to the ligand's centre of mass,
    #    and score the survivors sharing that anchor with Equation 1.
    # -------------------------------------------------------------------------
    _min_angle = _np.deg2rad(45.0)
    _max_angle = _np.deg2rad(135.0)

    angle_ok_idxs = [
        s_idx
        for s_idx in range(len(valid_sextuplets))
        if _min_angle < dof_tA[s_idx].mean() < _max_angle
        and _min_angle < dof_tB[s_idx].mean() < _max_angle
    ]

    if not angle_ok_idxs:
        raise ValueError(
            "All hydrogen-bonded Boresch sextuplets were rejected because one or "
            "both anchor angles have a mean outside the 45-135 degree RXRX window."
        )

    candidate_l1_idxs = {valid_sextuplets[s_idx]["l1_idx"] for s_idx in angle_ok_idxs}
    mean_com_dist = {
        idx: float(_np.mean(candidate_com_dists[idx])) for idx in candidate_l1_idxs
    }
    best_l1 = min(mean_com_dist, key=mean_com_dist.get)

    scored = []
    for s_idx in angle_ok_idxs:
        if valid_sextuplets[s_idx]["l1_idx"] != best_l1:
            continue

        mean_tA = dof_tA[s_idx].mean()
        mean_tB = dof_tB[s_idx].mean()
        _, std_pA = _circular_mean_std(dof_pA[s_idx])
        _, std_pB = _circular_mean_std(dof_pB[s_idx])
        _, std_pC = _circular_mean_std(dof_pC[s_idx])

        score = (
            dof_r[s_idx].std()
            * dof_tA[s_idx].std()
            * dof_tB[s_idx].std()
            * std_pA
            * std_pB
            * std_pC
            * dof_r[s_idx].mean() ** 2
            * (1.0 - abs(_np.sin(mean_tA)))
            * (1.0 - abs(_np.sin(mean_tB)))
        )
        scored.append((score, s_idx))

    scored.sort()

    if restraint_idx >= len(scored):
        raise ValueError(
            f"restraint_idx={restraint_idx} exceeds the number of valid "
            f"candidates ({len(scored)})."
        )

    _, best = scored[restraint_idx]
    s = valid_sextuplets[best]

    # -------------------------------------------------------------------------
    # 7. Compute equilibrium values (trajectory means, circular for dihedrals)
    #    and the RXRX protocol force constants.
    # -------------------------------------------------------------------------
    r0 = float(dof_r[best].mean())
    tA0 = float(dof_tA[best].mean())
    tB0 = float(dof_tB[best].mean())
    pA0, _ = _circular_mean_std(dof_pA[best])
    pB0, _ = _circular_mean_std(dof_pB[best])
    pC0, _ = _circular_mean_std(dof_pC[best])

    if force_constant_r is None:
        kr_val = 1.0
    else:
        kr_val = float(force_constant_r.to(r_fc_unit))

    if force_constant_angle is None:
        ka_val = 80.0
    else:
        ka_val = float(force_constant_angle.to(angle_fc_unit))

    kr_str = f"{kr_val:.6f} kcal mol-1 A-2"
    kt_str = [f"{ka_val:.6f} kcal mol-1 rad-2"] * 2
    kp_str = [f"{ka_val:.6f} kcal mol-1 rad-2"] * 3

    restraints, correction = _assemble_restraints(
        system,
        pert_mol_num,
        s,
        r0,
        tA0,
        tB0,
        pA0,
        pB0,
        pC0,
        kr_str,
        kt_str,
        kp_str,
        temperature,
        angle_potential,
        restraint_lever,
    )

    starting_structure = _best_starting_frame(
        system,
        _least_strained_frame(
            dof_r[best],
            dof_tA[best],
            dof_tB[best],
            dof_pA[best],
            dof_pB[best],
            dof_pC[best],
            r0,
            tA0,
            tB0,
            pA0,
            pB0,
            pC0,
            kr_val,
            ka_val,
            ka_val,
            ka_val,
            ka_val,
            ka_val,
        ),
    )

    return restraints, correction, starting_structure


# ---------------------------------------------------------------------------
# Aldeghi protocol (reference implementation): candidate sextuplets are
# generated from bulk ligand-receptor distance variance rather than
# hydrogen-bond occupancy, matching the MDRestraintsGenerator/BioSimSpace
# approach. Kept for comparison; RXRX is the default protocol.
# ---------------------------------------------------------------------------


def _boresch_search_aldeghi(
    system,
    temperature="298 K",
    receptor_selection="(not water) and (atomidx > 1) and (atomname CA, C, N)",
    cutoff="10 A",
    restraint_idx=0,
    force_constant=None,
    max_candidates=100,
    min_frames=50,
    angle_potential="restricted_bending",
    restraint_lever="combined",
):
    """
    Generate a Boresch restraint for an ABFE simulation by analysing a
    trajectory of the protein-ligand complex, using the Aldeghi-style
    (MDRestraintsGenerator/BioSimSpace) restraint search protocol.

    The system must contain trajectory frames (i.e. be the result of
    ``dynamics.commit()`` after a run with ``frame_frequency > 0``). The
    perturbable molecule is used as the ligand; the receptor anchor atoms are
    selected via ``receptor_selection``.

    The six Boresch degrees of freedom are sampled over all frames. Candidate
    sextuplets are scored by configurational volume (lower = tighter restraint)
    and the winner is used to construct a ``sire.mm.BoreschRestraints`` object.
    Force constants are derived from per-DOF trajectory variance via the
    equipartition theorem unless ``force_constant`` is given.

    Parameters
    ----------

    system : sire.system.System
        A Sire system with embedded trajectory frames. Must contain exactly
        one perturbable molecule.

    temperature : str or GeneralUnit, optional
        Simulation temperature. Defaults to 298 K.

    receptor_selection : str
        Sire selection string for receptor anchor atom candidates.
        Defaults to backbone heavy atoms (CA, C, N) in non-water molecules,
        which is appropriate for AMBER-format input.

    cutoff : str or GeneralUnit
        Maximum mean ligand-receptor anchor distance. Pairs whose mean
        distance exceeds this value are excluded.

    restraint_idx : int
        Index into the candidate list sorted by ascending configurational
        volume (0 = tightest).

    force_constant : str or GeneralUnit, optional
        Override for all force constants. If None (default), force constants
        are fitted from trajectory variance.

    max_candidates : int
        Maximum number of (l1, r1) pairs to evaluate for full DOF sampling
        in the second trajectory pass. The pairs with the lowest distance
        variance are evaluated first.

    min_frames : int
        Minimum number of trajectory frames required. Raise ``ValueError``
        if the trajectory has fewer frames. Default 50.

    angle_potential : str
        The functional form used for the two Boresch angle restraint terms,
        either "harmonic" or "restricted_bending" (see
        sire.restraints.boresch). Defaults to "restricted_bending".

    restraint_lever : str
        How the restraint's six degrees of freedom are grouped into
        lambda-addressable OpenMM Forces (see sire.restraints.boresch).
        Defaults to "combined", matching the Aldeghi protocol, where all
        six restraint terms are turned on together.

    Returns
    -------

    restraints : sire.mm.BoreschRestraints
        The generated Boresch restraints object, ready to pass to
        ``Config.restraints``.

    correction : sire.units.GeneralUnit
        Standard state correction in kcal mol-1.

    starting_structure : sire.system.System
        The trajectory frame at which the generated restraint is least
        strained (its six Boresch degrees of freedom are closest to the
        restraint equilibrium values). Because the equilibrium values are
        trajectory averages, the structure the search was seeded from is
        generally not consistent with them; seeding production from this
        frame instead avoids a large restraint force at t=0 that can
        otherwise destabilise the simulation as the restraint is switched on.
    """

    from ..legacy import Mol as _SireMol
    from ..system import System as _System
    from ..legacy.Units import k_boltz as _k_boltz

    # -------------------------------------------------------------------------
    # Parameter validation.
    # -------------------------------------------------------------------------
    if not isinstance(system, _System):
        raise TypeError(
            f"'system' must be of type 'sire.system.System', got {type(system)}"
        )

    temperature = _parse_unit(
        temperature, "temperature", _sr.units.kelvin, "temperature"
    )

    if not isinstance(receptor_selection, str):
        raise TypeError(
            f"'receptor_selection' must be a str, got {type(receptor_selection)}"
        )

    cutoff = _parse_unit(cutoff, "cutoff", _sr.units.angstrom, "length")

    if not isinstance(restraint_idx, int) or restraint_idx < 0:
        raise ValueError(
            f"'restraint_idx' must be a non-negative int, got {restraint_idx!r}"
        )

    if force_constant is not None:
        force_constant = _parse_unit(
            force_constant,
            "force_constant",
            _sr.units.kcal_per_mol / _sr.units.angstrom**2,
            "energy length-2",
        )

    if not isinstance(max_candidates, int) or max_candidates < 1:
        raise ValueError(
            f"'max_candidates' must be a positive int, got {max_candidates!r}"
        )

    if not isinstance(min_frames, int) or min_frames < 1:
        raise ValueError(f"'min_frames' must be a positive int, got {min_frames!r}")

    _validate_angle_potential(angle_potential)
    _validate_restraint_lever(restraint_lever)

    n_frames = system.num_frames()
    if n_frames < min_frames:
        raise ValueError(
            f"Trajectory has only {n_frames} frame(s); at least {min_frames} are required "
            "for reliable Boresch restraint generation."
        )

    # kBT as a plain float in kcal mol-1 (for equipartition).
    kBT = float((_k_boltz * temperature).to(_sr.units.kcal_per_mol))

    # Cutoff as a plain float in Angstroms.
    cutoff_A = float(cutoff.to(_sr.units.angstrom))

    # -------------------------------------------------------------------------
    # 1. Locate the perturbable molecule (ligand).
    # -------------------------------------------------------------------------
    pert_mols = system.molecules("property is_perturbable")
    if pert_mols.num_molecules() != 1:
        raise ValueError(
            "System must contain exactly one perturbable molecule for Boresch "
            f"restraint generation; found {pert_mols.num_molecules()}."
        )
    pert_mol = pert_mols.molecule(0)
    pert_mol_num = pert_mol.number()

    ghost_elem = _SireMol.Element(0)
    h_elem = _SireMol.Element("H")

    # Collect non-ghost, non-H ligand AtomIdx values (lambda=0 state).
    lig_atom_idxs = []
    for atom in pert_mol.atoms():
        elem0 = atom.property("element0")
        if elem0 != ghost_elem and elem0 != h_elem:
            lig_atom_idxs.append(atom.index())

    if len(lig_atom_idxs) < 3:
        raise ValueError(
            f"Ligand has only {len(lig_atom_idxs)} non-ghost, non-hydrogen atom(s); "
            "need at least 3 for Boresch restraints."
        )

    lig_connectivity = pert_mol.connectivity()

    # -------------------------------------------------------------------------
    # 2. Locate receptor atoms.
    # -------------------------------------------------------------------------
    try:
        rec_sel = system[receptor_selection]
    except Exception as e:
        raise ValueError(
            f"Could not apply receptor selection '{receptor_selection}': {e}"
        ) from e

    # (mol_num, AtomIdx-within-mol) for each receptor atom.
    rec_atom_info = []
    for atom in rec_sel.atoms():
        rec_atom_info.append((atom.molecule().number(), atom.index()))

    if len(rec_atom_info) < 3:
        raise ValueError(
            f"Receptor selection '{receptor_selection}' matched only "
            f"{len(rec_atom_info)} atom(s); need at least 3."
        )

    n_lig = len(lig_atom_idxs)
    n_rec = len(rec_atom_info)

    # -------------------------------------------------------------------------
    # 3. First trajectory pass: accumulate all ligand-receptor distances.
    # -------------------------------------------------------------------------
    frame_dists = []

    for frame in system.trajectory():
        space = frame.space()
        lig_frame = frame.molecule(pert_mol_num)
        frame_rec_mol_cache = {}
        d = _np.empty((n_lig, n_rec))
        for i, l_idx in enumerate(lig_atom_idxs):
            l_coord = lig_frame.atom(l_idx).coordinates()
            for j, (r_mol_num, r_idx) in enumerate(rec_atom_info):
                if r_mol_num not in frame_rec_mol_cache:
                    frame_rec_mol_cache[r_mol_num] = frame.molecule(r_mol_num)
                r_coord = frame_rec_mol_cache[r_mol_num].atom(r_idx).coordinates()
                d[i, j] = float(space.calc_dist(l_coord, r_coord))
        frame_dists.append(d)

    # Shape: (n_frames, n_lig, n_rec)
    dists = _np.array(frame_dists)
    mean_dists = dists.mean(axis=0)
    var_dists = dists.var(axis=0)

    # Candidate (l1, r1) pairs within cutoff, sorted by ascending variance.
    candidate_pairs = sorted(
        (var_dists[i, j], i, j)
        for i in range(n_lig)
        for j in range(n_rec)
        if mean_dists[i, j] <= cutoff_A
    )

    if not candidate_pairs:
        raise ValueError(
            f"No ligand-receptor atom pairs found within cutoff {cutoff}. "
            "Try increasing the cutoff or adjusting the receptor selection."
        )

    candidate_pairs = candidate_pairs[:max_candidates]

    # -------------------------------------------------------------------------
    # 4. Build valid sextuplets (l1,l2,l3, r1,r2,r3).
    #
    # Convention (matching sire.restraints.boresch ordering):
    #   r = dist(l1, r1)
    #   θA = angle(r2, r1, l1)
    #   θB = angle(r1, l1, l2)
    #   φA = dihedral(r3, r2, r1, l1)   -- r3 bonded to r2
    #   φB = dihedral(r2, r1, l1, l2)
    #   φC = dihedral(r1, l1, l2, l3)   -- l3 bonded to l2
    # -------------------------------------------------------------------------
    rec_mol_cache = {}

    def _rec_mol(mol_num):
        if mol_num not in rec_mol_cache:
            rec_mol_cache[mol_num] = system.molecule(mol_num)
        return rec_mol_cache[mol_num]

    valid_sextuplets = []

    for _, l1_pos, r1_pos in candidate_pairs:
        l1_idx = lig_atom_idxs[l1_pos]
        r1_mol_num, r1_idx = rec_atom_info[r1_pos]
        r1_mol = _rec_mol(r1_mol_num)
        r1_conn = r1_mol.connectivity()

        ligand_triplets = _build_triplets(
            lig_connectivity, l1_idx, pert_mol, True, ghost_elem, h_elem
        )
        receptor_triplets = _build_triplets(
            r1_conn, r1_idx, r1_mol, False, ghost_elem, h_elem
        )

        for l1_idx, l2_idx, l3_idx in ligand_triplets:
            for r1_idx, r2_idx, r3_idx in receptor_triplets:
                valid_sextuplets.append(
                    {
                        "l1_idx": l1_idx,
                        "l2_idx": l2_idx,
                        "l3_idx": l3_idx,
                        "r1_mol_num": r1_mol_num,
                        "r1_idx": r1_idx,
                        "r2_idx": r2_idx,
                        "r3_idx": r3_idx,
                    }
                )

    if not valid_sextuplets:
        raise ValueError(
            "Could not construct a valid Boresch sextuplet from the candidate "
            "pairs. Try increasing the cutoff or adjusting the receptor selection."
        )

    # -------------------------------------------------------------------------
    # 5. Second trajectory pass: sample all 6 DOFs for every valid sextuplet.
    # -------------------------------------------------------------------------
    dof_r, dof_tA, dof_tB, dof_pA, dof_pB, dof_pC = _sample_dofs(
        system, pert_mol_num, valid_sextuplets, n_frames
    )

    # -------------------------------------------------------------------------
    # 6. Score sextuplets by configurational volume; apply stability filter.
    #
    # score ∝ r₀² · |sin θA₀| · |sin θB₀| · σr · σθA · σθB · σφA · σφB · σφC
    # (from Boresch et al. 2003, J. Phys. Chem. B, Eqn 32, lower = tighter)
    #
    # Dihedral (φ) statistics use circular mean/std since these DOFs wrap at
    # +-180 degrees; a plain arithmetic std would be wildly inflated whenever
    # the true distribution straddles the wrap point.
    # -------------------------------------------------------------------------
    _min_angle = _np.deg2rad(10.0)
    _max_angle = _np.deg2rad(170.0)

    scored = []
    for s_idx in range(len(valid_sextuplets)):
        mean_tA = dof_tA[s_idx].mean()
        mean_tB = dof_tB[s_idx].mean()

        if not (
            _min_angle < mean_tA < _max_angle and _min_angle < mean_tB < _max_angle
        ):
            continue

        _, std_pA = _circular_mean_std(dof_pA[s_idx])
        _, std_pB = _circular_mean_std(dof_pB[s_idx])
        _, std_pC = _circular_mean_std(dof_pC[s_idx])

        score = (
            dof_r[s_idx].mean() ** 2
            * abs(_np.sin(mean_tA))
            * abs(_np.sin(mean_tB))
            * dof_r[s_idx].std()
            * dof_tA[s_idx].std()
            * dof_tB[s_idx].std()
            * std_pA
            * std_pB
            * std_pC
        )
        scored.append((score, s_idx))

    if not scored:
        raise ValueError(
            "All candidate Boresch sextuplets were rejected because one or both "
            "anchor angles have a mean near 0° or 180° (collinearity instability). "
            "Try increasing the cutoff or adjusting the receptor selection."
        )

    scored.sort()

    if restraint_idx >= len(scored):
        raise ValueError(
            f"restraint_idx={restraint_idx} exceeds the number of valid "
            f"candidates ({len(scored)})."
        )

    _, best = scored[restraint_idx]
    s = valid_sextuplets[best]

    # -------------------------------------------------------------------------
    # 7. Compute equilibrium values (trajectory means) and force constants.
    # -------------------------------------------------------------------------
    r0 = float(dof_r[best].mean())
    tA0 = float(dof_tA[best].mean())
    tB0 = float(dof_tB[best].mean())
    pA0, std_pA_best = _circular_mean_std(dof_pA[best])
    pB0, std_pB_best = _circular_mean_std(dof_pB[best])
    pC0, std_pC_best = _circular_mean_std(dof_pC[best])

    if force_constant is not None:
        # Use the numeric magnitude for all DOFs with appropriate units per DOF
        # type (A-2 for distance, rad-2 for angles). This matches BSS behaviour.
        fc_val = float(
            force_constant.to(_sr.units.kcal_per_mol / _sr.units.angstrom**2)
        )
        kr_str = f"{fc_val:.6f} kcal mol-1 A-2"
        kt_str = [f"{fc_val:.6f} kcal mol-1 rad-2"] * 2
        kp_str = [f"{fc_val:.6f} kcal mol-1 rad-2"] * 3
        # Per-DOF weights (kcal mol-1 units) for the least-strained-frame pick.
        kr_w = ktA_w = ktB_w = kpA_w = kpB_w = kpC_w = fc_val
    else:
        # Equipartition: k = kBT / (2σ²).
        # Sire's boresch() uses E = k·x² (half-spring-constant convention).
        def _k(sigma):
            return kBT / (2.0 * float(sigma) ** 2)

        kr_w = _k(dof_r[best].std())
        ktA_w = _k(dof_tA[best].std())
        ktB_w = _k(dof_tB[best].std())
        kpA_w = _k(std_pA_best)
        kpB_w = _k(std_pB_best)
        kpC_w = _k(std_pC_best)

        kr_str = f"{kr_w:.6f} kcal mol-1 A-2"
        kt_str = [
            f"{ktA_w:.6f} kcal mol-1 rad-2",
            f"{ktB_w:.6f} kcal mol-1 rad-2",
        ]
        kp_str = [
            f"{kpA_w:.6f} kcal mol-1 rad-2",
            f"{kpB_w:.6f} kcal mol-1 rad-2",
            f"{kpC_w:.6f} kcal mol-1 rad-2",
        ]

    # -------------------------------------------------------------------------
    # 8. Build the BoreschRestraints object and pick the least-strained frame.
    # -------------------------------------------------------------------------
    restraints, correction = _assemble_restraints(
        system,
        pert_mol_num,
        s,
        r0,
        tA0,
        tB0,
        pA0,
        pB0,
        pC0,
        kr_str,
        kt_str,
        kp_str,
        temperature,
        angle_potential,
        restraint_lever,
    )

    starting_structure = _best_starting_frame(
        system,
        _least_strained_frame(
            dof_r[best],
            dof_tA[best],
            dof_tB[best],
            dof_pA[best],
            dof_pB[best],
            dof_pC[best],
            r0,
            tA0,
            tB0,
            pA0,
            pB0,
            pC0,
            kr_w,
            ktA_w,
            ktB_w,
            kpA_w,
            kpB_w,
            kpC_w,
        ),
    )

    return restraints, correction, starting_structure


# ---------------------------------------------------------------------------
# Public entry point.
# ---------------------------------------------------------------------------


def boresch_search(
    system,
    protocol="rxrx",
    temperature="298 K",
    restraint_idx=0,
    min_frames=50,
    protein_selection="not water",
    hbond_distance_cutoff="3.5 A",
    hbond_angle_cutoff="120 degrees",
    occupancy_cutoff=0.5,
    force_constant_r=None,
    force_constant_angle=None,
    receptor_selection="(not water) and (atomidx > 1) and (atomname CA, C, N)",
    cutoff="10 A",
    force_constant=None,
    max_candidates=100,
    angle_potential="restricted_bending",
    restraint_lever=None,
):
    """
    Generate a Boresch restraint for an ABFE simulation by analysing a
    trajectory of the protein-ligand complex.

    Two restraint search protocols are available, selected via ``protocol``;
    the parameters below are grouped by which protocol(s) use them.

    - ``"rxrx"`` (default): the RXRX restraint search algorithm, which seeds
      candidate anchor atoms from protein-ligand hydrogen bonds and scores
      them with a formula that explicitly penalises angles near the Boresch
      0/180 degree singularities.

    - ``"aldeghi"``: a reference implementation of the MDRestraintsGenerator/
      BioSimSpace restraint search, kept for comparison.

    Parameters
    ----------

    system : sire.system.System
        A Sire system with embedded trajectory frames. Must contain exactly
        one perturbable molecule. Used by both protocols.

    protocol : str
        The restraint search protocol to use: ``"rxrx"`` (default) or
        ``"aldeghi"``.

    temperature : str or GeneralUnit, optional
        Simulation temperature. Defaults to 298 K. Used by both protocols.

    restraint_idx : int
        Index into the candidate list sorted by ascending score (0 = best).
        Used by both protocols.

    min_frames : int
        Minimum number of trajectory frames required. Used by both protocols.

    angle_potential : str
        The functional form used for the two Boresch angle restraint terms,
        either "harmonic" or "restricted_bending" (see
        sire.restraints.boresch). Defaults to "restricted_bending", which
        diverges as the angle approaches 0 or 180 degrees, preventing it
        from ever reaching the Boresch collinearity singularity. Used by
        both protocols.

    restraint_lever : str, optional
        How the restraint's six degrees of freedom are grouped into
        lambda-addressable OpenMM Forces (see sire.restraints.boresch).
        Defaults to None, which is matched to ``protocol``: "split" for
        "rxrx" (matching the RXRX protocol's staged restraint turn-on) and
        "combined" for "aldeghi" (matching the Aldeghi protocol, where all
        six restraint terms are turned on together). Used by both
        protocols.

    protein_selection : str
        [rxrx only] Sire selection string used to identify protein atoms
        considered as hydrogen-bond partners. Defaults to all non-water
        atoms.

    hbond_distance_cutoff : str or GeneralUnit
        [rxrx only] Maximum donor-acceptor heavy atom distance for a
        hydrogen bond.

    hbond_angle_cutoff : str or GeneralUnit
        [rxrx only] Minimum donor-H...acceptor angle for a hydrogen bond.

    occupancy_cutoff : float
        [rxrx only] Minimum fraction of frames (0-1) in which a hydrogen
        bond must be present for it to be considered as a restraint anchor.

    force_constant_r : str or GeneralUnit, optional
        [rxrx only] Distance restraint force constant. Defaults to the RXRX
        protocol value of 1 kcal mol-1 A-2.

    force_constant_angle : str or GeneralUnit, optional
        [rxrx only] Angle and dihedral restraint force constant. Defaults to
        the RXRX protocol value of 80 kcal mol-1 rad-2.

    receptor_selection : str
        [aldeghi only] Sire selection string for receptor anchor atom
        candidates. Defaults to backbone heavy atoms (CA, C, N) in non-water
        molecules, which is appropriate for AMBER-format input.

    cutoff : str or GeneralUnit
        [aldeghi only] Maximum mean ligand-receptor anchor distance. Pairs
        whose mean distance exceeds this value are excluded.

    force_constant : str or GeneralUnit, optional
        [aldeghi only] Override for all force constants. If None (default),
        force constants are fitted from trajectory variance.

    max_candidates : int
        [aldeghi only] Maximum number of (l1, r1) pairs to evaluate for full
        DOF sampling in the second trajectory pass. The pairs with the
        lowest distance variance are evaluated first.

    Returns
    -------

    restraints : sire.mm.BoreschRestraints
        The generated Boresch restraints object, ready to pass to
        ``Config.restraints``.

    correction : sire.units.GeneralUnit
        Standard state correction in kcal mol-1.

    starting_structure : sire.system.System
        The trajectory frame at which the generated restraint is least
        strained (its six Boresch degrees of freedom are closest to the
        restraint equilibrium values). Because the equilibrium values are
        trajectory averages, the structure the search was seeded from is
        generally not consistent with them; seeding production from this
        frame instead avoids a large restraint force at t=0 that can
        otherwise destabilise the simulation as the restraint is switched on.
    """
    if not isinstance(protocol, str):
        raise TypeError(f"'protocol' must be a str, got {type(protocol)}")

    if restraint_lever is None:
        # Matches the paper: the RXRX protocol turns the dihedral and
        # distance/angle restraint terms on according to different lambda
        # schedule equations, whereas the Aldeghi protocol turns all six
        # terms on together.
        restraint_lever = "split" if protocol == "rxrx" else "combined"

    if protocol == "rxrx":
        return _boresch_search_rxrx(
            system,
            temperature=temperature,
            protein_selection=protein_selection,
            hbond_distance_cutoff=hbond_distance_cutoff,
            hbond_angle_cutoff=hbond_angle_cutoff,
            occupancy_cutoff=occupancy_cutoff,
            restraint_idx=restraint_idx,
            force_constant_r=force_constant_r,
            force_constant_angle=force_constant_angle,
            min_frames=min_frames,
            angle_potential=angle_potential,
            restraint_lever=restraint_lever,
        )
    elif protocol == "aldeghi":
        return _boresch_search_aldeghi(
            system,
            temperature=temperature,
            receptor_selection=receptor_selection,
            cutoff=cutoff,
            restraint_idx=restraint_idx,
            force_constant=force_constant,
            max_candidates=max_candidates,
            min_frames=min_frames,
            angle_potential=angle_potential,
            restraint_lever=restraint_lever,
        )
    else:
        raise ValueError(
            f"Unknown 'protocol'={protocol!r}; must be 'rxrx' or 'aldeghi'."
        )
