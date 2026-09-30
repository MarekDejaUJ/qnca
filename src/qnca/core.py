"""Core estimator for Quantile Necessary Condition Analysis (QNCA).

QNCA generalises Dul's Necessary Condition Analysis (NCA). NCA draws a single
deterministic ceiling along the upper-left edge of an (X, Y) scatter and reads a
bottleneck off it. That ceiling is fixed by a handful of extreme points, so one
unusually efficient case can drag the whole envelope toward the floor. QNCA
replaces the single deterministic envelope by a *family* of frontiers indexed by
a tolerance parameter ``pi`` in (0, 1], where each frontier is a conditional
quantile of the condition among cases that reached an outcome level. The
deterministic NCA ceiling re-emerges as the ``pi == 1`` limit, which is the
method's validation anchor (see :func:`nca_ce_fdh_d`).

The estimator has four moving parts, kept deliberately small:

1. ``qnca_frontier``    -- the quantile necessity frontier phi_pi(y), isotonic-fit.
2. ``qnca_d``           -- the effect size d_pi: empty-zone area / scope.
3. ``permutation_test`` -- the identification p-value, weight-light, Dul-style.
4. ``qnca``             -- the user-facing wrapper returning all of the above.

Identification (does a necessity relation exist) is a property of the raw
scatter and is answered by the permutation test on uniform weights. Any
downstream causal/weighted claims belong to a different inferential layer and are
out of scope for this package by design.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from fractions import Fraction
from statistics import NormalDist

import numpy as np

__all__ = [
    "QNCAResult",
    "SpuriousnessBandResult",
    "ConsistencyProbeResult",
    "quantile_type1",
    "quantile_type1_pi",
    "qnca_rank",
    "isotonic_increasing",
    "qnca_frontier",
    "qnca_d",
    "permutation_test",
    "qnca",
    "nca_ce_fdh_d",
    "spuriousness_band",
    "consistency_probe",
]

_NORMAL = NormalDist()


@dataclass(frozen=True)
class QNCAResult:
    """Container for a single QNCA fit at one tolerance ``pi``."""

    pi: float
    d_pi: float
    p_pi: float | None
    y_grid: np.ndarray = field(repr=False)
    x_required: np.ndarray = field(repr=False)
    scope: float = field(repr=False)
    n_perm: int | None = None

    @property
    def bottleneck(self) -> np.ndarray:
        """Two-column (outcome_level, X_required) array: the bottleneck table."""
        return np.column_stack([self.y_grid, self.x_required])


@dataclass(frozen=True)
class SpuriousnessBandResult:
    """Result from the proposed pi-resolved spuriousness band."""

    pi_grid: np.ndarray = field(repr=False)
    d_obs: np.ndarray = field(repr=False)
    lower: np.ndarray = field(repr=False)
    median: np.ndarray = field(repr=False)
    upper: np.ndarray = field(repr=False)
    excess: np.ndarray = field(repr=False)
    p_null: np.ndarray = field(repr=False)
    d_null: np.ndarray = field(repr=False)
    rank_correlation: float

    @property
    def table(self) -> np.ndarray:
        """Seven-column table with one row per tolerance."""
        return np.column_stack([
            self.pi_grid,
            self.d_obs,
            self.lower,
            self.median,
            self.upper,
            self.excess,
            self.p_null,
        ])


@dataclass(frozen=True)
class ConsistencyProbeResult:
    """Result from the proposed extreme-vs-consistent consistency probe."""

    sizes: np.ndarray = field(repr=False)
    pi_values: np.ndarray = field(repr=False)
    mean_d: np.ndarray = field(repr=False)
    sd_d: np.ndarray = field(repr=False)
    beta_drift: np.ndarray = field(repr=False)
    divergence: float
    effects: np.ndarray = field(repr=False)


def quantile_type1(values: np.ndarray, prob: float) -> float:
    """Inverse-empirical-CDF quantile, matching R ``quantile(type = 1)``.

    Type 1 is the discontinuous inverse of the empirical CDF: the smallest order
    statistic x such that F(x) >= prob. We use it (rather than the default type 7
    linear interpolation) so that ``pi == 1`` returns the exact sample minimum and
    therefore reproduces NCA's CE-FDH boundary point-for-point.
    """
    x = np.asarray(values, dtype=float)
    if x.size == 0:
        return math.nan
    x = np.sort(x)
    if prob <= 0.0:
        return float(x[0])
    if prob >= 1.0:
        return float(x[-1])
    idx = int(math.ceil(x.size * prob)) - 1
    idx = min(max(idx, 0), x.size - 1)
    return float(x[idx])


def _tail_fraction(pi) -> Fraction:
    """Return ``1 - pi`` as an exact fraction of the configured decimal tolerance.

    A float tolerance is read through its shortest round-trip decimal, so
    ``0.95`` means nineteen twentieths exactly. Fractions are accepted as given.
    """
    if isinstance(pi, Fraction):
        p = pi
    else:
        value = float(pi)
        if not math.isfinite(value):
            raise ValueError("pi must be finite")
        p = Fraction(repr(value))
    if not (0 < p <= 1):
        raise ValueError("pi must be in (0, 1]")
    return 1 - p


def _rank_from_tail(k: int, q: Fraction) -> int:
    rank = -((-k * q.numerator) // q.denominator)
    return min(max(rank, 1), k)


def qnca_rank(k: int, pi) -> int:
    """Selected type-1 order-statistic rank ``max(1, ceil(k * (1 - pi)))``.

    The ceiling is evaluated exactly from the decimal tolerance, so ``pi = 0.95``
    with ``k = 100`` selects rank 5.
    """
    k = int(k)
    if k < 1:
        raise ValueError("conditioning-set size must be positive")
    return _rank_from_tail(k, _tail_fraction(pi))


def quantile_type1_pi(values: np.ndarray, pi) -> float:
    """Type-1 quantile parameterised by the QNCA tolerance ``pi``."""
    x = np.sort(np.asarray(values, dtype=float))
    if x.size == 0:
        return math.nan
    return float(x[qnca_rank(x.size, pi) - 1])


def _validate_xy(X: np.ndarray, Y: np.ndarray) -> None:
    if X.shape != Y.shape or X.ndim != 1:
        raise ValueError("X and Y must be one-dimensional with equal lengths")
    if X.size < 2:
        raise ValueError("X and Y must contain at least two observations")
    if not np.all(np.isfinite(X)):
        raise ValueError("X must contain only finite values")
    if not np.all(np.isfinite(Y)):
        raise ValueError("Y must contain only finite values")
    if not np.min(X) < np.max(X):
        raise ValueError("X must be non-degenerate")
    if not np.min(Y) < np.max(Y):
        raise ValueError("Y must be non-degenerate")


def _validated_scope(
    X: np.ndarray,
    Y: np.ndarray,
    scope: tuple[float, float, float, float] | None,
) -> tuple[float, float, float, float]:
    if scope is None:
        bounds = (float(np.min(X)), float(np.max(X)), float(np.min(Y)), float(np.max(Y)))
    else:
        bounds = tuple(float(v) for v in scope)
    if len(bounds) != 4:
        raise ValueError("scope must contain x_min, x_max, y_min, y_max")
    if not all(math.isfinite(v) for v in bounds):
        raise ValueError("scope bounds must be finite")
    x_min, x_max, y_min, y_max = bounds
    if not x_min < x_max:
        raise ValueError("scope requires x_min < x_max")
    if not y_min < y_max:
        raise ValueError("scope requires y_min < y_max")
    if np.min(X) < x_min or np.max(X) > x_max:
        raise ValueError("X observations must lie inside scope")
    if np.min(Y) < y_min or np.max(Y) > y_max:
        raise ValueError("Y observations must lie inside scope")
    return bounds  # type: ignore[return-value]


def _validate_grid(y_grid: np.ndarray) -> None:
    if y_grid.size == 0:
        raise ValueError("y_grid must not be empty")
    if not np.all(np.isfinite(y_grid)):
        raise ValueError("y_grid must contain only finite values")
    if np.any(np.diff(y_grid) < 0):
        raise ValueError("y_grid must be non-decreasing")


def isotonic_increasing(y: np.ndarray) -> np.ndarray:
    """Pool-adjacent-violators isotonic regression (non-decreasing, unit weights).

    This mirrors ``stats::isoreg`` and the Julia PAVA used in the sibling packages,
    so the three implementations agree to floating point. A genuine necessity
    relation makes the frontier climb with the outcome level; enforcing
    monotonicity here removes finite-sample wobble and matches the non-decreasing
    shape of Dul's ceiling.
    """
    y = np.asarray(y, dtype=float)
    vals: list[float] = []
    wts: list[int] = []
    for yi in y:
        vals.append(float(yi))
        wts.append(1)
        while len(vals) > 1 and vals[-2] > vals[-1]:
            v2, w2 = vals.pop(), wts.pop()
            v1, w1 = vals.pop(), wts.pop()
            vals.append((v1 * w1 + v2 * w2) / (w1 + w2))
            wts.append(w1 + w2)
    out = np.empty_like(y, dtype=float)
    pos = 0
    for v, w in zip(vals, wts):
        out[pos : pos + w] = v
        pos += w
    return out


def qnca_frontier(
    X: np.ndarray,
    Y: np.ndarray,
    pi: float,
    y_grid: np.ndarray,
    x_max: float | None = None,
) -> np.ndarray:
    """Quantile necessity frontier phi_pi(y), isotonic-fit to be non-decreasing.

    For each outcome level y on the grid, look only at the cases that reached it
    ({i : Y_i >= y}) and take the (1 - pi)-quantile of their condition values.
    At pi == 1 this is the per-level minimum -- Dul's CE-FDH leftmost boundary.
    Levels with no qualifying case are carried to x_max when a fixed x ceiling is
    supplied; otherwise they keep the last finite frontier value.
    """
    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float)
    y_grid = np.asarray(y_grid, dtype=float)
    _validate_xy(X, Y)
    _validate_grid(y_grid)
    return _frontier_validated(X, Y, _tail_fraction(pi), y_grid, x_max=x_max)


def _frontier_validated(
    X: np.ndarray,
    Y: np.ndarray,
    q: Fraction,
    y_grid: np.ndarray,
    x_max: float | None = None,
) -> np.ndarray:
    phi = np.full(y_grid.shape, np.nan, dtype=float)
    for j, y in enumerate(y_grid):
        hit = Y >= y
        if np.any(hit):
            x = np.sort(X[hit])
            phi[j] = x[_rank_from_tail(x.size, q) - 1]
    finite = np.isfinite(phi)
    if np.count_nonzero(finite) >= 2:
        phi[finite] = isotonic_increasing(phi[finite])
    # Scope reading: above the highest observed outcome the conditioning set is
    # empty. Under a fixed scope Dul's CE-FDH counts the whole top band as empty,
    # so the inverse frontier carries to the scope's x ceiling. Frontier-only
    # calls without an x ceiling keep the last finite value.
    idx_finite = np.flatnonzero(finite)
    if idx_finite.size:
        last = idx_finite[-1]
        if last + 1 < phi.size:
            phi[last + 1 :] = phi[last] if x_max is None else float(x_max)
    return phi


def qnca_d(phi: np.ndarray, y_grid: np.ndarray, x_min: float, scope: float) -> float:
    """Effect size d_pi: trapezoidal empty-zone area left of the frontier / scope.

    The empty zone at outcome level y has width ``phi(y) - x_min``; integrating
    that width over y and dividing by the scope rectangle gives a number in [0, 1]
    with the same reading as Dul's d.
    """
    w = np.maximum(0.0, phi - x_min)
    ok = np.isfinite(w)
    if np.count_nonzero(ok) < 2 or scope <= 0.0:
        return 0.0
    yg = y_grid[ok]
    wv = w[ok]
    area = float(np.sum(np.diff(yg) * (wv[:-1] + wv[1:]) / 2.0))
    return max(0.0, min(1.0, area / scope))


def permutation_test(
    X: np.ndarray,
    Y: np.ndarray,
    pi: float,
    d_obs: float,
    y_grid: np.ndarray,
    x_min: float,
    scope: float,
    B: int,
    rng: np.random.Generator,
    permutations: np.ndarray | None = None,
    x_max: float | None = None,
) -> float:
    """Finite-sample permutation p-value for identification.

    Null hypothesis: X and Y are unrelated. Shuffling Y against X breaks the
    dependence; the empty space that survives is what randomness alone produces.
    The (1 + .) / (B + 1) correction makes the observed sample count as one draw
    from the null, so the p-value is exact in the permutation sense and never
    exactly zero. This is Dul, van der Laan & Kuik (2020), transcribed to the
    quantile statistic, and it runs on uniform weights on purpose.
    """
    if permutations is not None:
        perm = _validate_permutation_matrix(permutations, Y.size)
        iterator = (Y[row - 1] for row in perm)
        B = int(perm.shape[0])
    else:
        iterator = (rng.permutation(Y) for _ in range(B))

    q = _tail_fraction(pi)
    count = 0
    y_range = y_grid[-1] - y_grid[0]
    frontier_x_max = x_max if x_max is not None else (
        x_min + scope / y_range if y_range > 0 else None
    )
    for Yp in iterator:
        d_p = qnca_d(
            _frontier_validated(X, Yp, q, y_grid, x_max=frontier_x_max),
            y_grid,
            x_min,
            scope,
        )
        if np.isfinite(d_p) and d_p >= d_obs:
            count += 1
    return (1.0 + count) / (B + 1.0)


def _validate_permutation_matrix(permutations: np.ndarray, n: int) -> np.ndarray:
    """Return a one-based permutation matrix after validating its shape."""
    perm = np.asarray(permutations, dtype=int)
    if perm.ndim == 1:
        perm = perm.reshape(1, -1)
    if perm.ndim != 2 or perm.shape[1] != n:
        raise ValueError("permutations must be a B x n matrix")
    expected = np.arange(1, n + 1)
    for row in perm:
        if not np.array_equal(np.sort(row), expected):
            raise ValueError("each permutation row must contain 1:n exactly once")
    return perm


def _scope_tuple(
    X: np.ndarray,
    Y: np.ndarray,
    scope: tuple[float, float, float, float] | None,
) -> tuple[float, float, float, float]:
    return _validated_scope(X, Y, scope)


def _d_for_pi(
    X: np.ndarray,
    Y: np.ndarray,
    pi: float,
    scope: tuple[float, float, float, float],
    n_grid: int,
) -> float:
    return qnca(X, Y, pi=pi, n_grid=n_grid, B=None, scope=scope).d_pi


def _average_ranks(values: np.ndarray) -> np.ndarray:
    x = np.asarray(values, dtype=float)
    order = np.argsort(x, kind="mergesort")
    ranks = np.empty(x.size, dtype=float)
    i = 0
    while i < x.size:
        j = i + 1
        while j < x.size and x[order[j]] == x[order[i]]:
            j += 1
        ranks[order[i:j]] = (i + 1 + j) / 2.0
        i = j
    return ranks


def _normal_ppf_array(p: np.ndarray) -> np.ndarray:
    p = np.clip(np.asarray(p, dtype=float), 1e-12, 1.0 - 1e-12)
    return np.fromiter((_NORMAL.inv_cdf(float(v)) for v in p), dtype=float, count=p.size)


def _normal_cdf_array(z: np.ndarray) -> np.ndarray:
    z = np.asarray(z, dtype=float)
    return np.fromiter((_NORMAL.cdf(float(v)) for v in z), dtype=float, count=z.size)


def _rank_normal_correlation(X: np.ndarray, Y: np.ndarray) -> float:
    n = X.size
    ax = _normal_ppf_array(_average_ranks(X) / (n + 1.0))
    ay = _normal_ppf_array(_average_ranks(Y) / (n + 1.0))
    sx = float(np.std(ax))
    sy = float(np.std(ay))
    if sx == 0.0 or sy == 0.0:
        return 0.0
    r = float(np.corrcoef(ax, ay)[0, 1])
    if not np.isfinite(r):
        return 0.0
    return max(-0.999999, min(0.999999, r))


def _empirical_quantile_values(values: np.ndarray, probs: np.ndarray) -> np.ndarray:
    x = np.sort(np.asarray(values, dtype=float))
    p = np.clip(np.asarray(probs, dtype=float), 0.0, 1.0)
    idx = np.ceil(x.size * p).astype(int) - 1
    idx = np.clip(idx, 0, x.size - 1)
    return x[idx]


def _quantile_type7(values: np.ndarray, prob: float) -> float:
    x = np.sort(np.asarray(values, dtype=float))
    if x.size == 0:
        return math.nan
    if x.size == 1 or prob <= 0.0:
        return float(x[0])
    if prob >= 1.0:
        return float(x[-1])
    h = 1.0 + (x.size - 1.0) * prob
    j = int(math.floor(h))
    gamma = h - j
    if j >= x.size:
        return float(x[-1])
    return float((1.0 - gamma) * x[j - 1] + gamma * x[j])


def qnca(
    X: np.ndarray,
    Y: np.ndarray,
    pi: float = 1.0,
    n_grid: int = 50,
    B: int | None = 1999,
    scope: tuple[float, float, float, float] | None = None,
    seed: int | None = None,
    permutations: np.ndarray | None = None,
) -> QNCAResult:
    """Fit QNCA at a single tolerance ``pi``.

    Parameters
    ----------
    X, Y : array-like
        Condition and outcome, on a comparable scale.
    pi : float
        Tolerance in (0, 1]. pi == 1 reproduces Dul's deterministic CE-FDH ceiling.
    n_grid : int
        Number of outcome levels on the frontier grid.
    B : int or None
        Permutation repetitions. If None, the p-value is skipped (fast path for a
        bottleneck-only call).
    scope : (x_min, x_max, y_min, y_max) or None
        Fixed scope rectangle. If None, the observed data range is used.
    seed : int or None
        Seed for the permutation RNG.
    permutations : array-like or None
        Optional one-based B x n permutation matrix. When supplied, these
        permutations are used instead of drawing from the RNG.
    """
    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float)
    _validate_xy(X, Y)
    if int(n_grid) < 2:
        raise ValueError("n_grid must be at least 2")
    if B is not None and B < 0:
        raise ValueError("B must be non-negative or None")
    q = _tail_fraction(pi)
    x_min, x_max, y_min, y_max = _validated_scope(X, Y, scope)
    scope_area = (x_max - x_min) * (y_max - y_min)
    y_grid = np.linspace(y_min, y_max, int(n_grid))

    phi_obs = _frontier_validated(X, Y, q, y_grid, x_max=x_max)
    d_obs = qnca_d(phi_obs, y_grid, x_min, scope_area)

    p_pi: float | None = None
    n_perm: int | None = B
    if permutations is not None:
        perm = _validate_permutation_matrix(permutations, Y.size)
        if perm.shape[0] > 0:
            rng = np.random.default_rng(seed)
            p_pi = permutation_test(
                X, Y, pi, d_obs, y_grid, x_min, scope_area, perm.shape[0], rng, perm,
                x_max=x_max
            )
            n_perm = int(perm.shape[0])
    elif B is not None and B > 0:
        rng = np.random.default_rng(seed)
        p_pi = permutation_test(
            X, Y, pi, d_obs, y_grid, x_min, scope_area, B, rng, x_max=x_max
        )

    return QNCAResult(
        pi=pi,
        d_pi=d_obs,
        p_pi=p_pi,
        y_grid=y_grid,
        x_required=phi_obs,
        scope=scope_area,
        n_perm=n_perm,
    )


def nca_ce_fdh_d(
    X: np.ndarray,
    Y: np.ndarray,
    scope: tuple[float, float, float, float] | None = None,
) -> float:
    """Native CE-FDH effect size, computed the way Dul's NCA defines it.

    NCA's CE-FDH ceiling is the upper-left free-disposal hull of the scatter; the
    effect size d is the area of the empty zone above that ceiling, divided by the
    scope. We build the hull directly from the data (no quantiles, no grid) and
    integrate the empty area exactly as a sum of rectangles. This lets the QNCA
    validation invariant -- d_pi(pi = 1) == this number, in the continuous reading
    -- be checked without calling the R NCA package.

    The CE-FDH ceiling, expressed as a non-decreasing required-X function of the
    outcome y, is r(y) = min{ X_i : Y_i >= y }. The empty zone is everything to
    the left of that step function within scope, so
        C_E = integral over y of ( r(y) - x_min ) dy
    evaluated exactly on the breakpoints of the step function.
    """
    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float)
    _validate_xy(X, Y)
    x_min, x_max, y_min, y_max = _validated_scope(X, Y, scope)
    scope_area = (x_max - x_min) * (y_max - y_min)
    if scope_area <= 0.0:
        return 0.0

    # r(y) = min{ X_i : Y_i >= y } is a right-continuous step function that only
    # changes at the observed Y values. Sort points by Y descending; run_min[j] is
    # the smallest X among the j+1 highest-outcome cases. As y falls we admit more
    # cases, so r(y) is non-decreasing in y (matching the isotonic frontier).
    order = np.argsort(-Y)  # descending Y
    Ys = Y[order]
    Xs = X[order]
    run_min = np.minimum.accumulate(Xs)  # run_min[j] = min(Xs[0..j])

    # We integrate ( r(y) - x_min ) dy from y_min to y_max as a sum of horizontal
    # bands. On the band y in ( Ys[k], Ys[k-1] ], the cases reaching that outcome
    # are exactly the k highest ones (indices 0..k-1), so r = run_min[k-1]. The
    # topmost band ( Ys[0], y_max ] uses run_min[0]; below the lowest outcome all
    # cases qualify, so r = run_min[-1] (the global min X) down to y_min.
    def clamp_y(v: float) -> float:
        return min(max(v, y_min), y_max)

    area = 0.0
    # top band: above the highest observed outcome Dul's CE-FDH ceiling leaves
    # the whole fixed-scope band empty.
    top_lo = clamp_y(Ys[0])
    if y_max > top_lo:
        area += (y_max - top_lo) * max(0.0, x_max - x_min)
    # interior bands between consecutive distinct outcomes
    for k in range(1, len(Ys)):
        hi = clamp_y(Ys[k - 1])
        lo = clamp_y(Ys[k])
        if hi > lo:
            area += (hi - lo) * max(0.0, run_min[k - 1] - x_min)
    # bottom band: from the lowest outcome down to the scope floor (all cases qualify)
    bot_hi = clamp_y(Ys[-1])
    if bot_hi > y_min:
        area += (bot_hi - y_min) * max(0.0, run_min[-1] - x_min)

    d = area / scope_area
    return max(0.0, min(1.0, d))


def spuriousness_band(
    X: np.ndarray,
    Y: np.ndarray,
    pi_grid: tuple[float, ...] | list[float] | np.ndarray = (1.0, 0.95, 0.90),
    n_grid: int = 50,
    M: int = 999,
    scope: tuple[float, float, float, float] | None = None,
    seed: int | None = None,
    normal_draws: np.ndarray | None = None,
    uniform_draws: np.ndarray | None = None,
) -> SpuriousnessBandResult:
    """Proposed pi-resolved spuriousness band.

    The diagnostic compares the observed ``d_pi`` curve with datasets generated
    from a Gaussian copula at the observed normal-score rank correlation and the
    empirical marginals of ``X`` and ``Y``. Supplying ``uniform_draws`` uses those
    empirical-copula uniforms directly, which is useful for cross-language parity.
    It is a descriptive diagnostic with unestablished operating characteristics;
    it does not change the QNCA estimator.
    """
    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float)
    _validate_xy(X, Y)
    if int(n_grid) < 2:
        raise ValueError("n_grid must be at least 2")
    pi_vals = np.asarray(pi_grid, dtype=float)
    if pi_vals.size == 0:
        raise ValueError("pi_grid must not be empty")
    for pi in pi_vals:
        _tail_fraction(float(pi))
    scope_tuple = _scope_tuple(X, Y, scope)

    if normal_draws is not None and uniform_draws is not None:
        raise ValueError("supply at most one of normal_draws and uniform_draws")
    if uniform_draws is not None:
        uniforms = np.asarray(uniform_draws, dtype=float)
        if uniforms.ndim != 3 or uniforms.shape[1] != X.size or uniforms.shape[2] != 2:
            raise ValueError("uniform_draws must have shape M x n x 2")
        M = int(uniforms.shape[0])
        draws = None
    elif normal_draws is None:
        rng = np.random.default_rng(seed)
        draws = rng.standard_normal((int(M), X.size, 2))
        uniforms = None
    else:
        draws = np.asarray(normal_draws, dtype=float)
        if draws.ndim != 3 or draws.shape[1] != X.size or draws.shape[2] != 2:
            raise ValueError("normal_draws must have shape M x n x 2")
        M = int(draws.shape[0])
        uniforms = None
    if M <= 0:
        raise ValueError("M must be positive")

    d_obs = np.array([_d_for_pi(X, Y, pi, scope_tuple, n_grid) for pi in pi_vals])
    d_null = np.empty((int(M), pi_vals.size), dtype=float)
    r = _rank_normal_correlation(X, Y)
    scale = math.sqrt(max(0.0, 1.0 - r * r))

    for j in range(int(M)):
        if uniforms is not None:
            Xj = _empirical_quantile_values(X, uniforms[j, :, 0])
            Yj = _empirical_quantile_values(Y, uniforms[j, :, 1])
        else:
            zx = draws[j, :, 0]
            zy = r * zx + scale * draws[j, :, 1]
            Xj = _empirical_quantile_values(X, _normal_cdf_array(zx))
            Yj = _empirical_quantile_values(Y, _normal_cdf_array(zy))
        for k, pi in enumerate(pi_vals):
            d_null[j, k] = _d_for_pi(Xj, Yj, float(pi), scope_tuple, n_grid)

    lower = np.array([_quantile_type7(d_null[:, k], 0.025) for k in range(pi_vals.size)])
    median = np.array([_quantile_type7(d_null[:, k], 0.500) for k in range(pi_vals.size)])
    upper = np.array([_quantile_type7(d_null[:, k], 0.975) for k in range(pi_vals.size)])
    excess = d_obs - upper
    p_null = (1.0 + np.sum(d_null >= d_obs.reshape(1, -1), axis=0)) / (int(M) + 1.0)

    return SpuriousnessBandResult(
        pi_grid=pi_vals,
        d_obs=d_obs,
        lower=lower,
        median=median,
        upper=upper,
        excess=excess,
        p_null=p_null,
        d_null=d_null,
        rank_correlation=r,
    )


def consistency_probe(
    X: np.ndarray,
    Y: np.ndarray,
    pi_pair: tuple[float, float] = (1.0, 0.90),
    sizes: tuple[int, ...] | list[int] | np.ndarray | None = None,
    reps: int = 100,
    n_grid: int = 50,
    scope: tuple[float, float, float, float] | None = None,
    seed: int | None = None,
    subsamples: list[np.ndarray] | None = None,
) -> ConsistencyProbeResult:
    """Proposed extreme-vs-consistent consistency probe.

    The diagnostic tracks how ``d_pi`` changes across increasing subsample sizes.
    It is descriptive and can be confounded by hard support bounds and marginal
    skewness, so it should be read beside the spuriousness band rather than alone.
    """
    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float)
    n = X.size
    pi_vals = np.asarray(pi_pair, dtype=float)
    if pi_vals.size != 2:
        raise ValueError("pi_pair must contain exactly two tolerances")
    if sizes is None:
        raw_sizes = np.array([max(10, n // 4), max(10, n // 2), max(10, 3 * n // 4), n])
        size_vals = np.unique(np.clip(raw_sizes, 2, n)).astype(int)
    else:
        size_vals = np.asarray(sizes, dtype=int)
    if size_vals.size < 2:
        raise ValueError("at least two subsample sizes are required")
    if np.any(size_vals < 2) or np.any(size_vals > n):
        raise ValueError("subsample sizes must be between 2 and n")
    if reps <= 0:
        raise ValueError("reps must be positive")

    scope_tuple = _scope_tuple(X, Y, scope)
    rng = np.random.default_rng(seed)
    effects = np.empty((size_vals.size, int(reps), pi_vals.size), dtype=float)

    if subsamples is not None and len(subsamples) != size_vals.size:
        raise ValueError("subsamples must contain one matrix per subsample size")

    for k, m in enumerate(size_vals):
        if subsamples is None:
            rows = np.vstack([
                rng.choice(n, size=int(m), replace=False) + 1 for _ in range(int(reps))
            ])
        else:
            rows = np.asarray(subsamples[k], dtype=int)
            if rows.ndim == 1:
                rows = rows.reshape(1, -1)
            if rows.shape != (int(reps), int(m)):
                raise ValueError("each subsample matrix must be reps x size")
            if np.any(rows < 1) or np.any(rows > n):
                raise ValueError("subsample indices must be one-based and within 1:n")
        for r_idx in range(int(reps)):
            idx = rows[r_idx, :] - 1
            for p_idx, pi in enumerate(pi_vals):
                effects[k, r_idx, p_idx] = _d_for_pi(X[idx], Y[idx], float(pi), scope_tuple, n_grid)

    mean_d = np.mean(effects, axis=1)
    sd_d = np.std(effects, axis=1, ddof=1) if reps > 1 else np.zeros_like(mean_d)
    g = 1.0 / size_vals.astype(float)
    denom = float(np.sum((g - np.mean(g)) ** 2))
    beta = np.full(pi_vals.size, math.nan, dtype=float)
    if denom > 0.0:
        for p_idx in range(pi_vals.size):
            ybar = mean_d[:, p_idx]
            beta[p_idx] = float(np.sum((g - np.mean(g)) * (ybar - np.mean(ybar))) / denom)
    divergence = float(abs(beta[0]) - abs(beta[1])) if np.all(np.isfinite(beta)) else math.nan

    return ConsistencyProbeResult(
        sizes=size_vals,
        pi_values=pi_vals,
        mean_d=mean_d,
        sd_d=sd_d,
        beta_drift=beta,
        divergence=divergence,
        effects=effects,
    )
