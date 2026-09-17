"""
Regression tests for SireAnalysis::TIPMF (ti.cpp), which fits TI
gradient-vs-lambda data to a polynomial and integrates it. This used to be
implemented on top of the third_party/regress.{h,cpp} library (GPLv3,
Conrad Shyu), now replaced with an Eigen-based least-squares fit; these
tests give confidence that rewrite didn't change behaviour.

Only 3 of the old Regress class's methods were ever actually used by
TIPMF::recalculate():
  - the constructor/fit itself (classic normal-equations polynomial least
    squares - a Vandermonde/moment matrix solved via Gaussian elimination)
  - GetPolynomial() - feeds TIPMF::integral() (the last cumulative value of
    an analytic power-rule integration of the fit, see ti.cpp ~L294-361) and
    TIPMF::smoothedGradients() (the fit evaluated at ~101 points across the
    range)
  - DoQuadrature() - feeds TIPMF::quadrature(), a **plain trapezoidal rule
    on the raw (possibly endpoint-extended) gradient points**, independent
    of the polynomial fit entirely

So there are two independent things to check against independent oracles:
  - integral()/smoothedGradients() against numpy.polynomial.polynomial's own
    least-squares polyfit (SVD-based, not the same algorithm as Regress's
    Gaussian elimination on the normal equations, so compared with a loose-
    ish tolerance rather than expecting bit-identical agreement - the new
    Eigen-based fit uses a similarly more-stable method, not a reproduction
    of Regress's exact algorithm)
  - quadrature() against a direct, hand-written trapezoidal rule, replicating
    ti.cpp's endpoint-extension logic (if the integration range extends
    beyond the first/last raw gradient point, a synthetic point holding the
    nearest gradient's y value is added at that end before integrating)
"""

import numpy as np
import numpy.polynomial.polynomial as npoly
import pytest

import sire.legacy.Analysis as A


@pytest.fixture(autouse=True, scope="module")
def _ensure_pythonized():
    # pythonize TIPMF directly, rather than relying on some other test
    # module triggering the global use_new_api() sweep first (unreliable)
    # or doing it here at collection time (too early - would clash with
    # test files that must be the first to pick an API mode, e.g.
    # tests/biosimspace/test_select.py's use_mixed_api())
    from sire._pythonize import _pythonize

    _pythonize(A.TIPMF, delete_old=True)


def _fit_coeffs(x, y, degree):
    return npoly.polyfit(x, y, degree)


def _analytic_integral(coeffs, xmin, xmax):
    total = 0.0
    for i, c in enumerate(coeffs):
        power = i + 1
        total += c * (xmax**power - xmin**power) / power
    return total


def _trapezoidal(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    return float(np.sum((y[1:] + y[:-1]) * 0.5 * (x[1:] - x[:-1])))


def _extend_with_endpoints(x, y, range_min, range_max):
    x = list(x)
    y = list(y)

    if x[0] > range_min:
        x.insert(0, range_min)
        y.insert(0, y[0])

    if x[-1] < range_max:
        x.append(range_max)
        y.append(y[-1])

    return x, y


def _make_pmf(x, y, range_min, range_max, order):
    dps = [A.DataPoint(float(xi), float(yi)) for xi, yi in zip(x, y)]
    pmf = A.TIPMF(range_min, range_max, order)
    pmf.set_gradients(dps)
    return pmf


# (x values, y = f(x), range_min, range_max, polynomial order)
DATASETS = [
    (np.linspace(0.0, 1.0, 6), lambda x: x**2 + 0.1 * x, 0.0, 1.0, 2),
    (np.linspace(0.0, 1.0, 11), lambda x: np.sin(3 * x) + 0.2, 0.0, 1.0, 4),
    # gradients don't span the full range - exercises the endpoint-extension
    # logic in both quadrature() and (indirectly) integral()
    (np.linspace(0.1, 0.9, 9), lambda x: x**3 - 0.5 * x, 0.0, 1.0, 3),
    (np.linspace(0.0, 1.0, 21), lambda x: 2.0 - x, 0.0, 1.0, 1),
    (np.array([0.0, 0.25, 0.5, 0.75, 1.0]), lambda x: np.exp(x) - 1, 0.0, 1.0, 3),
    (np.linspace(0.2, 0.8, 4), lambda x: np.cos(x), 0.0, 1.0, 2),
]


@pytest.mark.parametrize("x_raw,y_fn,range_min,range_max,order", DATASETS)
def test_ti_pmf_integral_matches_polyfit_oracle(
    x_raw, y_fn, range_min, range_max, order
):
    x = np.asarray(x_raw, dtype=float)
    y = y_fn(x)

    pmf = _make_pmf(x, y, range_min, range_max, order)

    coeffs = _fit_coeffs(x, y, order)
    expected_integral = _analytic_integral(coeffs, range_min, range_max)

    assert pmf.integral() == pytest.approx(expected_integral, rel=2e-3, abs=1e-6)


@pytest.mark.parametrize("x_raw,y_fn,range_min,range_max,order", DATASETS)
def test_ti_pmf_quadrature_matches_trapezoidal_oracle(
    x_raw, y_fn, range_min, range_max, order
):
    x = np.asarray(x_raw, dtype=float)
    y = y_fn(x)

    pmf = _make_pmf(x, y, range_min, range_max, order)

    ext_x, ext_y = _extend_with_endpoints(x, y, range_min, range_max)
    expected_quad = _trapezoidal(ext_x, ext_y)

    assert pmf.quadrature() == pytest.approx(expected_quad, rel=1e-8, abs=1e-10)


@pytest.mark.parametrize("x_raw,y_fn,range_min,range_max,order", DATASETS)
def test_ti_pmf_smoothed_gradients_match_polyfit_oracle(
    x_raw, y_fn, range_min, range_max, order
):
    x = np.asarray(x_raw, dtype=float)
    y = y_fn(x)

    pmf = _make_pmf(x, y, range_min, range_max, order)
    coeffs = _fit_coeffs(x, y, order)

    smoothed = pmf.smoothed_gradients()
    assert len(smoothed) > 0

    # the range endpoints are the most numerically sensitive (extrapolation-
    # adjacent), so explicitly check those are present and correct as well
    # as everything in between
    xs = [dp.x() for dp in smoothed]
    assert xs[0] == pytest.approx(range_min)
    assert xs[-1] == pytest.approx(range_max)

    for dp in smoothed:
        expected_y = npoly.polyval(dp.x(), coeffs)
        assert dp.y() == pytest.approx(expected_y, rel=2e-3, abs=1e-6)


def test_ti_pmf_order_and_range_accessors():
    x = np.linspace(0.0, 1.0, 6)
    y = x**2

    pmf = _make_pmf(x, y, 0.0, 1.0, 2)

    assert pmf.order() == 2
    assert pmf.range_min() == pytest.approx(0.0)
    assert pmf.range_max() == pytest.approx(1.0)


def test_ti_pmf_empty_gives_zero():
    pmf = A.TIPMF(0.0, 1.0, 2)

    assert pmf.integral() == 0
    assert pmf.quadrature() == 0
