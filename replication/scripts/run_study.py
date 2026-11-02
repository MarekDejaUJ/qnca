#!/usr/bin/env python3
"""Reproduce every QNCA table the paper reports, from the installed package.

Outputs (written to --out, default ../data):
  qnca_dataset_n{N}.csv ............ the shipped reverse-L datasets
  validation_invariant.csv ......... QNCA d_pi(1) vs native CE-FDH vs (optional) Dul NCA
  pi_grid_results.csv .............. d_pi and p_pi across the pi grid and dataset sizes
  bottleneck_n1000_pi095.csv ....... a practitioner bottleneck table
  outlier_diagnostic.csv ........... d_pi(1) vs d_pi(0.95) with/without a planted efficient case
  calibration_null.csv ............. permutation false-positive rate under no signal
  power_curve.csv .................. replicated power as the empty corner sharpens
  nca_dul_bridge.R + dul_nca_results.csv (only if Rscript + NCA are present)

All seeds are fixed so the paper reports exact, reproducible numbers.
Run:  python scripts/run_study.py --B 1999
"""

from __future__ import annotations

import argparse
import csv
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

from qnca import generate_reverse_l, nca_ce_fdh_d, qnca
from qnca.nca_bridge import run_r_nca

SCOPE = (0.0, 100.0, 0.0, 100.0)
SIZES = [100, 200, 500, 1000]
PI_GRID = [1.00, 0.95, 0.90]


def _write(path: Path, header, rows) -> None:
    with path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)


def make_datasets(out: Path, seed: int) -> dict[int, tuple[np.ndarray, np.ndarray]]:
    rng = np.random.default_rng(seed)
    data = {}
    for n in SIZES:
        X, Y = generate_reverse_l(n, rng=rng)
        data[n] = (X, Y)
        _write(out / f"qnca_dataset_n{n}.csv", ["X", "Y"], zip(X, Y))
    return data


def validation_invariant(out: Path, data) -> None:
    rows = []
    for n in SIZES:
        X, Y = data[n]
        d1 = qnca(X, Y, pi=1.0, n_grid=400, B=None, scope=SCOPE).d_pi
        native = nca_ce_fdh_d(X, Y, scope=SCOPE)
        rows.append([f"n{n}", n, f"{d1:.6f}", f"{native:.6f}", f"{abs(d1 - native):.2e}"])
    _write(out / "validation_invariant.csv",
           ["dataset", "n", "qnca_d_pi1", "native_ce_fdh_d", "abs_diff"], rows)


def pi_grid(out: Path, data, B: int) -> None:
    rows = []
    for n in SIZES:
        X, Y = data[n]
        for pi in PI_GRID:
            res = qnca(X, Y, pi=pi, n_grid=50, B=B, scope=SCOPE, seed=1000 + n)
            rows.append([f"n{n}", n, f"{pi:.2f}", f"{res.d_pi:.6f}", f"{res.p_pi:.6f}"])
    _write(out / "pi_grid_results.csv",
           ["dataset", "n", "pi", "d_pi", "p_pi"], rows)


def bottleneck(out: Path, data) -> None:
    X, Y = data[1000]
    res = qnca(X, Y, pi=0.95, n_grid=50, B=None, scope=SCOPE)
    sel = np.floor(np.linspace(1, len(res.y_grid), 11)).astype(int) - 1
    rows = [[f"{res.y_grid[j]:.4f}", f"{res.x_required[j]:.4f}"] for j in sel]
    _write(out / "bottleneck_n1000_pi095.csv", ["outcome_level", "X_required"], rows)


def outlier_diagnostic(out: Path, seed: int) -> None:
    """Plant one very efficient case (high Y, tiny X) and show the gap react.

    Deterministic NCA / pi=1 is dragged down by the planted case; pi<1 discards it
    and the gap d_pi(0.95) - d_pi(1) widens -- the diagnostic the paper features.
    """
    rng = np.random.default_rng(seed)
    X, Y = generate_reverse_l(300, rng=rng)
    rows = []
    for label, (Xx, Yy) in {
        "clean": (X, Y),
        "with_efficient_outlier": (np.append(X, 2.0), np.append(Y, 95.0)),
    }.items():
        d1 = qnca(Xx, Yy, pi=1.00, n_grid=200, B=None, scope=SCOPE).d_pi
        d95 = qnca(Xx, Yy, pi=0.95, n_grid=200, B=None, scope=SCOPE).d_pi
        rows.append([label, f"{d1:.6f}", f"{d95:.6f}", f"{d95 - d1:.6f}"])
    _write(out / "outlier_diagnostic.csv",
           ["case", "d_pi_1.00", "d_pi_0.95", "gap"], rows)


def _calibration_trial(args: tuple[float, int, int]) -> tuple[float, bool]:
    pi, seed, B = args
    rng = np.random.default_rng(seed)
    X = rng.uniform(0, 100, 150)
    Y = rng.uniform(0, 100, 150)
    res = qnca(X, Y, pi=pi, n_grid=40, B=B, scope=SCOPE, seed=seed + 1_000_000)
    return pi, bool(res.p_pi is not None and res.p_pi < 0.05)


def _power_trial(args: tuple[float, int, int]) -> tuple[float, bool, float]:
    strength, seed, B = args
    rng = np.random.default_rng(seed)
    X, Y = generate_ceiling_strength_sample(120, strength, rng)
    res = qnca(X, Y, pi=0.95, n_grid=40, B=B, scope=SCOPE, seed=seed + 1_000_000)
    return strength, bool(res.p_pi is not None and res.p_pi < 0.05), float(res.d_pi)


def _run_tasks(worker, tasks, workers: int):
    if workers <= 1:
        return [worker(task) for task in tasks]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        return list(pool.map(worker, tasks))


def calibration_null(out: Path, seed: int, B: int, reps: int = 300, workers: int = 1) -> None:
    """False-positive rate of the permutation test when X, Y are independent."""
    rng = np.random.default_rng(seed)
    rows = []
    for pi in PI_GRID:
        seeds = rng.integers(1, 1_000_000_000, size=reps, dtype=np.int64)
        tasks = [(pi, int(seed_i), B) for seed_i in seeds]
        results = _run_tasks(_calibration_trial, tasks, workers)
        rejects = sum(reject for _, reject in results)
        rows.append([f"{pi:.2f}", reps, rejects, f"{rejects / reps:.4f}"])
    _write(out / "calibration_null.csv",
           ["pi", "reps", "rejections_at_0.05", "fpr"], rows)


def generate_ceiling_strength_sample(
    n: int,
    strength: float,
    rng: np.random.Generator,
    noise: float = 10.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Generate data whose attainable ceiling tightens as strength increases."""
    X = rng.uniform(0, 100, n)
    ceiling = (1.0 - strength) * 100.0 + strength * X
    Y = rng.uniform(0, ceiling) + rng.normal(0, noise, n)
    return X, np.clip(Y, 0, 100)


def power_curve(out: Path, seed: int, B: int, reps: int = 150, workers: int = 1) -> None:
    """Replicated power as the attainable empty corner sharpens."""
    rng = np.random.default_rng(seed)
    rows = []
    alpha = 0.05
    for strength in [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]:
        seeds = rng.integers(1, 1_000_000_000, size=reps, dtype=np.int64)
        tasks = [(strength, int(seed_i), B) for seed_i in seeds]
        results = _run_tasks(_power_trial, tasks, workers)
        rejects = sum(reject for _, reject, _ in results)
        effects = [effect for _, _, effect in results]
        rows.append([
            f"{strength:.1f}", reps, B, f"{alpha:.2f}", rejects,
            f"{rejects / reps:.4f}", f"{float(np.mean(effects)):.6f}",
        ])
    _write(
        out / "power_curve.csv",
        ["ceiling_strength", "reps", "B", "alpha", "rejections", "power", "mean_d_pi_0.95"],
        rows,
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--B", type=int, default=1999)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--calibration-B", type=int, default=999)
    ap.add_argument("--calibration-reps", type=int, default=300)
    ap.add_argument("--power-B", type=int, default=499)
    ap.add_argument("--power-reps", type=int, default=150)
    ap.add_argument("--workers", type=int, default=max(1, min(4, os.cpu_count() or 1)))
    ap.add_argument("--out", type=Path, default=Path(__file__).resolve().parent.parent / "data")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    data = make_datasets(args.out, args.seed)
    validation_invariant(args.out, data)
    pi_grid(args.out, data, args.B)
    bottleneck(args.out, data)
    outlier_diagnostic(args.out, args.seed + 5)
    calibration_null(
        args.out, args.seed + 6, args.calibration_B,
        reps=args.calibration_reps, workers=args.workers,
    )
    power_curve(
        args.out, args.seed + 7, args.power_B,
        reps=args.power_reps, workers=args.workers,
    )

    ok, msg = run_r_nca(args.out, SIZES, test_rep=args.B)
    bridge_status = (
        "OK - dul_nca_results.csv written."
        if ok else
        f"SKIPPED - {msg}"
    )
    (args.out / "run_status.txt").write_text(
        f"QNCA study complete. B={args.B}, seed={args.seed}.\n"
        f"Calibration: reps={args.calibration_reps}, B={args.calibration_B}.\n"
        f"Power: reps={args.power_reps}, B={args.power_B}.\n"
        f"Workers: {args.workers}.\n"
        f"Dul NCA bridge: {bridge_status}\n"
    )
    print("Study complete. Outputs in", args.out)
    print("Dul NCA bridge:", msg)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
