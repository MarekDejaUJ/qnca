"""Test suite for the QNCA Python package.

The headline test is the validation invariant: d_pi at pi = 1 equals the native
CE-FDH effect size (the continuous reading of Dul's d). The remaining tests pin
down the structural properties that make the method behave like a necessity
analysis: a non-decreasing frontier, an effect size in [0, 1], monotone growth of
d_pi as pi falls (outlier discarding raises the floor), and a permutation p-value
bounded below by 1 / (B + 1).
"""

from __future__ import annotations

import numpy as np
import pytest

from qnca import (
    consistency_probe,
    generate_reverse_l,
    isotonic_increasing,
    nca_ce_fdh_d,
    qnca,
    qnca_frontier,
    spuriousness_band,
)

SCOPE = (0.0, 100.0, 0.0, 100.0)


@pytest.fixture(params=[100, 200, 500, 1000])
def data(request):
    rng = np.random.default_rng(42)
    return generate_reverse_l(request.param, rng=rng)


def test_validation_invariant_pi1_equals_ce_fdh(data):
    """d_pi(pi = 1) must match the native CE-FDH d on a fine grid.

    Checked under both the fixed [0,100]^2 scope (with frontier carry-forward to
    the scope ceiling) and the data-range scope. The residual is grid
    discretization only, so a tight tolerance is appropriate.
    """
    X, Y = data
    for scope in (SCOPE, None):
        res = qnca(X, Y, pi=1.0, n_grid=400, B=None, scope=scope)
        d_ref = nca_ce_fdh_d(X, Y, scope=scope)
        assert abs(res.d_pi - d_ref) < 2e-3


def test_frontier_is_non_decreasing(data):
    X, Y = data
    y_grid = np.linspace(0.0, 100.0, 50)
    phi = qnca_frontier(X, Y, pi=0.95, y_grid=y_grid)
    finite = phi[np.isfinite(phi)]
    assert np.all(np.diff(finite) >= -1e-9)


def test_effect_size_in_unit_interval(data):
    X, Y = data
    for pi in (1.0, 0.95, 0.90):
        res = qnca(X, Y, pi=pi, B=None, scope=SCOPE)
        assert 0.0 <= res.d_pi <= 1.0


def test_d_pi_increases_as_pi_falls(data):
    """Discarding efficient outliers (lower pi) raises the floor, enlarging d."""
    X, Y = data
    d1 = qnca(X, Y, pi=1.0, n_grid=200, B=None, scope=SCOPE).d_pi
    d95 = qnca(X, Y, pi=0.95, n_grid=200, B=None, scope=SCOPE).d_pi
    d90 = qnca(X, Y, pi=0.90, n_grid=200, B=None, scope=SCOPE).d_pi
    assert d1 <= d95 + 1e-9 <= d90 + 1e-9


def test_permutation_p_lower_bound():
    X, Y = generate_reverse_l(200, rng=np.random.default_rng(7))
    B = 199
    res = qnca(X, Y, pi=0.95, B=B, scope=SCOPE, seed=1)
    assert res.p_pi is not None
    assert res.p_pi >= 1.0 / (B + 1.0) - 1e-12
    assert res.p_pi <= 1.0


def test_supplied_permutations_control_p_value():
    X, Y = generate_reverse_l(50, rng=np.random.default_rng(11))
    permutations = np.tile(np.arange(1, len(Y) + 1), (5, 1))
    res = qnca(X, Y, pi=0.95, B=None, scope=SCOPE, permutations=permutations)
    assert res.p_pi == 1.0
    assert res.n_perm == 5


def test_isotonic_known_case():
    x = np.array([3.0, 1.0, 2.0, 5.0, 4.0])
    out = isotonic_increasing(x)
    assert np.all(np.diff(out) >= -1e-12)
    assert abs(out.sum() - x.sum()) < 1e-9  # PAVA preserves the total


def test_no_signal_gives_small_effect():
    """Random X, Y (no necessity) should give a small effect size."""
    rng = np.random.default_rng(0)
    X = rng.uniform(0, 100, 500)
    Y = rng.uniform(0, 100, 500)
    res = qnca(X, Y, pi=1.0, n_grid=200, B=None, scope=SCOPE)
    assert res.d_pi < 0.10


def test_spuriousness_band_shape_and_bounds():
    X, Y = generate_reverse_l(80, rng=np.random.default_rng(21))
    draws = np.random.default_rng(3).standard_normal((12, len(X), 2))
    res = spuriousness_band(
        X, Y, pi_grid=(1.0, 0.95), n_grid=40, M=12, scope=SCOPE, normal_draws=draws
    )
    assert res.table.shape == (2, 7)
    assert res.d_null.shape == (12, 2)
    assert np.all(np.isfinite(res.table))
    assert np.all(res.lower <= res.median + 1e-12)
    assert np.all(res.median <= res.upper + 1e-12)
    assert np.all((res.p_null >= 1.0 / 13.0 - 1e-12) & (res.p_null <= 1.0))


def test_spuriousness_band_accepts_uniform_draws():
    X, Y = generate_reverse_l(60, rng=np.random.default_rng(31))
    uniforms = np.random.default_rng(5).uniform(size=(9, len(X), 2))
    res = spuriousness_band(
        X, Y, pi_grid=(1.0, 0.95), n_grid=35, M=9, scope=SCOPE, uniform_draws=uniforms
    )
    assert res.d_null.shape == (9, 2)
    assert np.all(np.isfinite(res.table))


def test_spuriousness_band_rejects_two_supplied_streams():
    X, Y = generate_reverse_l(30, rng=np.random.default_rng(32))
    uniforms = np.random.default_rng(6).uniform(size=(4, len(X), 2))
    normals = np.random.default_rng(7).standard_normal((4, len(X), 2))
    with pytest.raises(ValueError, match="at most one"):
        spuriousness_band(
            X, Y, M=4, scope=SCOPE, uniform_draws=uniforms, normal_draws=normals
        )


def test_consistency_probe_shape_and_finite_slopes():
    X, Y = generate_reverse_l(90, rng=np.random.default_rng(22))
    res = consistency_probe(
        X, Y, pi_pair=(1.0, 0.90), sizes=(30, 60, 90), reps=8,
        n_grid=40, scope=SCOPE, seed=4
    )
    assert res.mean_d.shape == (3, 2)
    assert res.sd_d.shape == (3, 2)
    assert res.effects.shape == (3, 8, 2)
    assert np.all(np.isfinite(res.mean_d))
    assert np.all(np.isfinite(res.sd_d))
    assert np.all(np.isfinite(res.beta_drift))
    assert np.isfinite(res.divergence)


def test_decimal_tolerance_selects_stated_rank():
    from fractions import Fraction

    from qnca import qnca_rank, quantile_type1_pi

    assert qnca_rank(100, 0.95) == 5
    assert qnca_rank(20, 0.95) == 1
    assert qnca_rank(21, 0.95) == 2
    assert qnca_rank(10, 0.90) == 1
    assert qnca_rank(11, 0.90) == 2
    assert qnca_rank(100, 1.0) == 1
    assert qnca_rank(100, Fraction(19, 20)) == 5
    assert quantile_type1_pi(np.arange(1.0, 101.0), 0.95) == 5.0
    for k in range(1, 501):
        assert qnca_rank(k, 0.95) == max(1, -(-k // 20))
        assert qnca_rank(k, 0.90) == max(1, -(-k // 10))
    with pytest.raises(ValueError):
        qnca_rank(100, 0.0)
    with pytest.raises(ValueError):
        qnca_rank(100, 1.01)


def test_orientation_sentinel_distinguishes_x_from_y():
    X = np.array([0.10, 0.20, 0.40, 0.80, 0.90])
    Y = np.array([0.10, 0.20, 0.60, 0.70, 0.95])
    assert qnca_frontier(X, Y, 1.0, np.array([0.60]))[0] == 0.40
    assert qnca_frontier(Y, X, 1.0, np.array([0.60]))[0] == 0.70


def test_frontier_uses_corrected_rank_at_decimal_boundary():
    v = np.arange(1.0, 101.0)
    assert qnca_frontier(v, v, 0.95, np.array([1.0]))[0] == 5.0


def test_invalid_inputs_and_scopes_fail_explicitly():
    with pytest.raises(ValueError):
        qnca([1.0, 2.0], [1.0], B=None)
    with pytest.raises(ValueError):
        qnca([1.0, np.nan], [1.0, 2.0], B=None)
    with pytest.raises(ValueError):
        qnca([1.0, 1.0], [1.0, 2.0], B=None)
    with pytest.raises(ValueError):
        qnca([1.0, 2.0], [1.0, 2.0], B=None, scope=(0.0, 1.5, 0.0, 2.0))
    with pytest.raises(ValueError):
        qnca([1.0, 2.0], [1.0, 2.0], B=None, n_grid=1)


def test_zero_valued_cohorts_remain_valid():
    res = qnca([0.0, 0.2, 0.4, 0.8], [0.0, 0.25, 0.5, 0.75], pi=0.95,
               n_grid=101, B=None, scope=(0.0, 1.0, 0.0, 1.0))
    assert np.isfinite(res.d_pi)
    assert 0.0 <= res.d_pi <= 1.0


def test_envelope_leaves_strict_member_unchanged():
    rng = np.random.default_rng(401)
    for n in (60, 300):
        X, Y = generate_reverse_l(n, rng=rng)
        a = qnca(X, Y, pi=1.0, n_grid=120, B=None, scope=SCOPE)
        b = qnca(X, Y, pi=1.0, n_grid=120, B=None, scope=SCOPE, monotone="isotonic")
        assert np.array_equal(a.x_required, b.x_required)
        assert a.d_pi == b.d_pi
    from qnca import monotone_envelope

    assert monotone_envelope([3.0, 1.0, 4.0, 2.0, 5.0]).tolist() == [3.0, 3.0, 4.0, 4.0, 5.0]
    with pytest.raises(ValueError):
        qnca([1.0, 2.0], [1.0, 2.0], B=None, monotone="median")


def test_tolerance_path_is_monotone_in_pi():
    rng = np.random.default_rng(402)
    grid = [1.0, 0.99, 0.975, 0.95, 0.925, 0.90, 0.85, 0.80, 0.50]
    for _ in range(10):
        X, Y = generate_reverse_l(150, rng=rng)
        for mode in ("envelope", "isotonic"):
            d = [qnca(X, Y, pi=p, n_grid=60, B=None, scope=SCOPE, monotone=mode).d_pi for p in grid]
            assert np.all(np.diff(d) >= -1e-12)


def test_resolution_table_agrees_with_fitted_frontier():
    from qnca import qnca_resolution

    rng = np.random.default_rng(404)
    X, Y = generate_reverse_l(250, rng=rng)
    Xa = np.append(X, 2.0)
    Ya = np.append(Y, 95.0)
    for mode in ("envelope", "isotonic"):
        res = qnca_resolution(Xa, Ya, pi=0.95, n_grid=80, scope=SCOPE, monotone=mode)
        fit = qnca(Xa, Ya, pi=0.95, n_grid=80, B=None, scope=SCOPE, monotone=mode)
        assert np.array_equal(res.fitted, fit.x_required)
        ok = res.k > 0
        assert np.all(res.resistance[ok] == res.rank[ok] - 1)
        if mode == "envelope":
            keep = ok & ~res.inherited
            assert np.all(res.below[keep] <= res.rank[keep] - 1)
            assert np.any(res.inherited)


def test_deletion_screen_across_tolerance_grid():
    from qnca import qnca_outliers

    rng = np.random.default_rng(405)
    X, Y = generate_reverse_l(600, rng=rng)
    Xa = np.append(X, 1.0)
    Ya = np.append(Y, 50.0)
    out = qnca_outliers(Xa, Ya, pi_grid=(1.0, 0.95), scope=SCOPE, n_grid=60)
    assert len(out.cases) == 601
    assert out.cases[0] == [601]
    assert out.flagged[0]
    assert out.dif_abs[0, 0] > 10 * abs(out.dif_abs[0, 1])
    joint = qnca_outliers(Xa[:200], Ya[:200], pi_grid=(1.0, 0.95), k=2, scope=SCOPE,
                          n_grid=60, max_candidates=6)
    assert len(joint.cases) == 15
    assert all(len(c) == 2 for c in joint.cases)


def test_breakdown_curve_under_added_efficient_cases():
    from qnca import qnca_sensitivity

    rng = np.random.default_rng(406)
    X, Y = generate_reverse_l(300, rng=rng)
    res = qnca_sensitivity(X, Y, points=(3.0, 60.0), counts=range(7), scope=SCOPE, n_grid=80)
    assert np.all(res.retention[0] == 1.0)
    assert np.all(np.diff(res.d, axis=0) <= 1e-12)
    assert res.retention[1, 0] < res.retention[1, 1]
    with pytest.raises(ValueError):
        qnca_sensitivity(X, Y, points=(150.0, 60.0), scope=SCOPE)
