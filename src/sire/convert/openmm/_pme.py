__all__ = ["tune_pme"]


def tune_pme(mols, target_error=None, max_spacing=0.16, map=None, **kwargs):
    """
    Find the smallest PME grid, and the splitting parameter for it, whose
    forces are at least as accurate as the target. Since the real-space cost
    doesn't depend on the splitting parameter, the smallest grid is the fastest.
    This is intended for the CUDA and OpenCL platforms.

    Parameters
    ----------

    mols: sire.system.System, sire.mol.SelectorMol
        The molecules to tune for. These should have sensible (e.g.
        equilibrated) coordinates, since clashes dominate the force error.

    target_error: float
        The target RMS force error, relative to the RMS reference force. By
        default this is the error of the PME parameters that OpenMM would
        choose from the Ewald error tolerance.

    max_spacing: float
        The coarsest grid spacing to consider, in nanometers.

    map: dict
        The property map, as passed to `dynamics`.

    kwargs:
        Any other arguments to pass to `dynamics`, e.g. the cutoff, lambda
        value and platform, so that the system matches the simulation.

    Returns
    -------

    dict
        The `pme_alpha` and `pme_grid` map options to use. This is empty
        if no grid smaller than OpenMM's choice meets the target.
    """
    import math

    import numpy as np
    from openmm import NonbondedForce, unit

    d = mols.dynamics(map=map, **kwargs)
    context = d._d._omm_mols

    nbff = None

    for force in context.getSystem().getForces():
        if isinstance(force, NonbondedForce):
            nbff = force

    if nbff is None or nbff.getNonbondedMethod() != NonbondedForce.PME:
        raise ValueError("PME tuning requires a periodic system that uses PME.")

    cutoff = nbff.getCutoffDistance().value_in_unit(unit.nanometer)
    tolerance = nbff.getEwaldErrorTolerance()

    box = [
        np.linalg.norm(v.value_in_unit(unit.nanometer))
        for v in context.getSystem().getDefaultPeriodicBoxVectors()
    ]

    force_unit = unit.kilojoule_per_mole / unit.nanometer

    def get_forces():
        state = context.getState(getForces=True)
        return state.getForces(asNumpy=True).value_in_unit(force_unit)

    def set_pme(alpha, grid):
        nbff.setPMEParameters(alpha, *grid)
        context.reinitialize(preserveState=True)

    _, *default_grid = nbff.getPMEParametersInContext(context)
    default_forces = get_forces()

    # reference forces from a tolerance an order of magnitude below the target
    ref_tolerance = min(tolerance, target_error or tolerance) / 10
    nbff.setEwaldErrorTolerance(ref_tolerance)
    set_pme(0.0, [0, 0, 0])
    ref_forces = get_forces()
    nbff.setEwaldErrorTolerance(tolerance)

    ref_norm = np.mean(np.sum(ref_forces**2, axis=1))

    def error(forces):
        return math.sqrt(np.mean(np.sum((forces - ref_forces) ** 2, axis=1)) / ref_norm)

    if target_error is None:
        target_error = error(default_forces)

    def grid_for(n):
        # grid with n points along the longest box vector, at the same spacing in the others
        spacing = max(box) / n
        return [max(6, math.ceil(length / spacing)) for length in box]

    def is_fft_friendly(n):
        for f in (2, 3, 5, 7):
            while n % f == 0:
                n //= f
        return n == 1

    n_lo = max(6, math.ceil(max(box) / max_spacing))
    n_hi = default_grid[box.index(max(box))]
    candidates = [n for n in range(n_lo, n_hi + 1) if is_fft_friendly(n)]

    # bracket the splitting parameter between OpenMM's choice for the target
    # and that for a real-space error twenty times smaller
    alpha_lo = math.sqrt(-math.log(2.0 * target_error)) / cutoff
    alpha_hi = math.sqrt(-math.log(0.1 * target_error)) / cutoff

    def best_alpha(grid, iterations=5):
        ratio = (math.sqrt(5) - 1) / 2

        def evaluate(alpha):
            set_pme(alpha, grid)
            return error(get_forces()), alpha

        a, b = alpha_lo, alpha_hi
        c = b - ratio * (b - a)
        e = a + ratio * (b - a)
        fc, fe = evaluate(c), evaluate(e)
        best = min(fc, fe)

        for _ in range(iterations):
            if fc[0] < fe[0]:
                b, e, fe = e, c, fc
                c = b - ratio * (b - a)
                fc = evaluate(c)
                best = min(best, fc)
            else:
                a, c, fc = c, e, fe
                e = a + ratio * (b - a)
                fe = evaluate(e)
                best = min(best, fe)

        return best

    best = None
    lo, hi = 0, len(candidates) - 1

    while lo <= hi:
        mid = (lo + hi) // 2
        grid = grid_for(candidates[mid])
        err, alpha = best_alpha(grid)

        if err <= target_error:
            best = (alpha, grid)
            hi = mid - 1
        else:
            lo = mid + 1

    del d

    if best is None:
        return {}

    alpha, grid = best

    return {"pme_alpha": alpha, "pme_grid": grid}
