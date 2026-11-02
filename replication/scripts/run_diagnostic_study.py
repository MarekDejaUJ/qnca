#!/usr/bin/env python3
"""Run the synthetic diagnostic validation workflow for QNCA.

Outputs (written to --out, default ../data):
  diagnostic_band_illustration.csv
  diagnostic_band_sample_size.csv
  diagnostic_probe_sample_size.csv
  diagnostic_cross_language_agreement.csv
  calibration_reconciliation.csv
  qnca_diagnostics.log

The study keeps the core QNCA outputs untouched and writes only new diagnostic
artefacts. The spuriousness band is screened against a reverse-L floor design
and a positive-dependence Gaussian-copula design with no necessity floor. The
consistency probe is logged alongside the same generators and treated as a
diagnostic screen rather than a calibrated test.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from statistics import NormalDist

import numpy as np

from qnca import consistency_probe, generate_reverse_l, qnca, spuriousness_band

ROOT = Path(__file__).resolve().parent.parent.parent
SCOPE = (0.0, 100.0, 0.0, 100.0)
PI_GRID = (1.00, 0.95, 0.90)
PI_PAIR = (1.00, 0.90)
SIZES = (100, 200, 500, 1000)
_NORMAL = NormalDist()


def _write_csv(path: Path, header: list[str], rows: list[list[object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)


def _run_tasks(worker, tasks, workers: int):
    if workers <= 1:
        return [worker(task) for task in tasks]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        return list(pool.map(worker, tasks))


def _normal_cdf_array(z: np.ndarray) -> np.ndarray:
    z = np.asarray(z, dtype=float)
    return np.fromiter((_NORMAL.cdf(float(v)) for v in z), dtype=float, count=z.size)


def gaussian_copula_no_floor(
    n: int,
    rng: np.random.Generator,
    rho: float = 0.70,
) -> tuple[np.ndarray, np.ndarray]:
    z1 = rng.standard_normal(n)
    z2 = rho * z1 + math.sqrt(max(0.0, 1.0 - rho * rho)) * rng.standard_normal(n)
    X = 100.0 * _normal_cdf_array(z1)
    Y = 100.0 * _normal_cdf_array(z2)
    return X, Y


def generate_strength_zero(
    n: int,
    rng: np.random.Generator,
    noise: float = 10.0,
) -> tuple[np.ndarray, np.ndarray]:
    X = rng.uniform(0.0, 100.0, n)
    Y = rng.uniform(0.0, 100.0, n) + rng.normal(0.0, noise, n)
    return X, np.clip(Y, 0.0, 100.0)


def independent_uniform(
    n: int,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    return rng.uniform(0.0, 100.0, n), rng.uniform(0.0, 100.0, n)


def make_dataset(
    generator: str,
    n: int,
    seed: int,
) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    if generator == "floor":
        return generate_reverse_l(n, rng=rng)
    if generator == "no_floor":
        return gaussian_copula_no_floor(n, rng=rng)
    raise ValueError(f"unknown generator {generator}")


def diagnostic_seed(base_seed: int, generator: str, n: int, rep: int = 0) -> int:
    offset = 0 if generator == "floor" else 10_000_000
    return int(base_seed + offset + n * 1000 + rep)


def band_screen_task(args: tuple[str, int, int, int, int, int]) -> tuple[str, int, int, np.ndarray, np.ndarray, np.ndarray]:
    generator, n, rep, M, n_grid, base_seed = args
    X, Y = make_dataset(generator, n, diagnostic_seed(base_seed, generator, n, rep))
    result = spuriousness_band(
        X, Y, pi_grid=PI_GRID, n_grid=n_grid, M=M, scope=SCOPE,
        seed=diagnostic_seed(base_seed + 500_000, generator, n, rep),
    )
    return generator, n, rep, result.d_obs, result.excess, result.p_null


def probe_screen_task(args: tuple[str, int, int, int, int, int]) -> tuple[str, int, int, np.ndarray, np.ndarray, np.ndarray, float]:
    generator, n, rep, subsample_reps, n_grid, base_seed = args
    X, Y = make_dataset(generator, n, diagnostic_seed(base_seed, generator, n, rep))
    sizes = np.unique(np.clip(np.array([n // 4, n // 2, 3 * n // 4, n]), 2, n)).astype(int)
    result = consistency_probe(
        X, Y, pi_pair=PI_PAIR, sizes=sizes, reps=subsample_reps,
        n_grid=n_grid, scope=SCOPE, seed=diagnostic_seed(base_seed + 900_000, generator, n, rep),
    )
    return generator, n, rep, sizes, result.mean_d, result.beta_drift, result.divergence


def calibration_task(args: tuple[str, int, float, int, int, int]) -> tuple[str, int, float, bool]:
    design, n, pi, B, seed, n_grid = args
    rng = np.random.default_rng(seed)
    if design == "independent_uniform":
        X, Y = independent_uniform(n, rng)
    elif design == "strength_zero":
        X, Y = generate_strength_zero(n, rng)
    else:
        raise ValueError(f"unknown design {design}")
    res = qnca(X, Y, pi=pi, n_grid=n_grid, B=B, scope=SCOPE, seed=seed + 1_000_000)
    return design, n, pi, bool(res.p_pi is not None and res.p_pi < 0.05)


def run_band_illustration(
    out: Path,
    seed: int,
    M: int,
    n_grid: int,
) -> dict[str, dict[float, dict[str, float]]]:
    rows: list[list[object]] = []
    summary: dict[str, dict[float, dict[str, float]]] = {}
    for generator in ("floor", "no_floor"):
        X, Y = make_dataset(generator, 500, diagnostic_seed(seed, generator, 500, 0))
        result = spuriousness_band(
            X, Y, pi_grid=PI_GRID, n_grid=n_grid, M=M, scope=SCOPE,
            seed=diagnostic_seed(seed + 500_000, generator, 500, 0),
        )
        summary[generator] = {}
        for idx, pi in enumerate(result.pi_grid):
            summary[generator][float(pi)] = {
                "d_obs": float(result.d_obs[idx]),
                "upper": float(result.upper[idx]),
                "delta": float(result.excess[idx]),
                "p_null": float(result.p_null[idx]),
            }
            rows.append([
                generator,
                500,
                f"{float(pi):.2f}",
                M,
                f"{float(result.d_obs[idx]):.6f}",
                f"{float(result.lower[idx]):.6f}",
                f"{float(result.median[idx]):.6f}",
                f"{float(result.upper[idx]):.6f}",
                f"{float(result.excess[idx]):.6f}",
                f"{float(result.p_null[idx]):.6f}",
                f"{float(result.rank_correlation):.6f}",
            ])
    _write_csv(
        out / "diagnostic_band_illustration.csv",
        [
            "generator", "n", "pi", "M", "d_obs", "lower", "median", "upper",
            "delta", "p_null", "rank_correlation",
        ],
        rows,
    )
    return summary


def run_band_sample_size(
    out: Path,
    seed: int,
    dataset_reps: int,
    M: int,
    n_grid: int,
    workers: int,
) -> dict[tuple[str, int], dict[float, dict[str, float]]]:
    tasks = [
        (generator, n, rep, M, n_grid, seed)
        for generator in ("floor", "no_floor")
        for n in SIZES
        for rep in range(dataset_reps)
    ]
    results = _run_tasks(band_screen_task, tasks, workers)
    grouped: dict[tuple[str, int], list[tuple[np.ndarray, np.ndarray, np.ndarray]]] = {}
    for generator, n, _, d_obs, excess, p_null in results:
        grouped.setdefault((generator, n), []).append((d_obs, excess, p_null))

    rows: list[list[object]] = []
    summary: dict[tuple[str, int], dict[float, dict[str, float]]] = {}
    for generator in ("floor", "no_floor"):
        for n in SIZES:
            triples = grouped[(generator, n)]
            d_mat = np.vstack([triple[0] for triple in triples])
            ex_mat = np.vstack([triple[1] for triple in triples])
            p_mat = np.vstack([triple[2] for triple in triples])
            summary[(generator, n)] = {}
            for idx, pi in enumerate(PI_GRID):
                positive_delta = ex_mat[:, idx] > 0.0
                small_p = p_mat[:, idx] <= 0.05
                detection = positive_delta & small_p
                summary[(generator, n)][pi] = {
                    "mean_d_obs": float(np.mean(d_mat[:, idx])),
                    "mean_delta": float(np.mean(ex_mat[:, idx])),
                    "mean_p_null": float(np.mean(p_mat[:, idx])),
                    "detection_rate": float(np.mean(detection)),
                }
                rows.append([
                    generator,
                    n,
                    f"{pi:.2f}",
                    dataset_reps,
                    M,
                    f"{float(np.mean(d_mat[:, idx])):.6f}",
                    f"{float(np.mean(ex_mat[:, idx])):.6f}",
                    f"{float(np.mean(p_mat[:, idx])):.6f}",
                    f"{float(np.mean(positive_delta)):.4f}",
                    f"{float(np.mean(small_p)):.4f}",
                    f"{float(np.mean(detection)):.4f}",
                ])
    _write_csv(
        out / "diagnostic_band_sample_size.csv",
        [
            "generator", "n", "pi", "dataset_reps", "M", "mean_d_obs",
            "mean_delta", "mean_p_null", "positive_delta_rate", "tail_rate",
            "detection_rate",
        ],
        rows,
    )
    return summary


def run_probe_sample_size(
    out: Path,
    seed: int,
    dataset_reps: int,
    subsample_reps: int,
    n_grid: int,
    workers: int,
) -> dict[tuple[str, int], dict[str, float]]:
    tasks = [
        (generator, n, rep, subsample_reps, n_grid, seed)
        for generator in ("floor", "no_floor")
        for n in SIZES
        for rep in range(dataset_reps)
    ]
    results = _run_tasks(probe_screen_task, tasks, workers)
    grouped: dict[tuple[str, int], list[tuple[np.ndarray, np.ndarray, float]]] = {}
    for generator, n, _, sizes, mean_d, beta_drift, divergence in results:
        grouped.setdefault((generator, n), []).append((sizes, mean_d, beta_drift, divergence))

    rows: list[list[object]] = []
    summary: dict[tuple[str, int], dict[str, float]] = {}
    for generator in ("floor", "no_floor"):
        for n in SIZES:
            items = grouped[(generator, n)]
            sizes = items[0][0]
            beta_mat = np.vstack([item[2] for item in items])
            divergence = np.array([item[3] for item in items], dtype=float)
            mean_d_small = np.array([item[1][0, :] for item in items], dtype=float)
            mean_d_large = np.array([item[1][-1, :] for item in items], dtype=float)
            summary[(generator, n)] = {
                "mean_beta_pi1": float(np.mean(beta_mat[:, 0])),
                "mean_beta_pi090": float(np.mean(beta_mat[:, 1])),
                "mean_divergence": float(np.mean(divergence)),
                "positive_divergence_rate": float(np.mean(divergence > 0.0)),
            }
            rows.append([
                generator,
                n,
                dataset_reps,
                subsample_reps,
                int(sizes[0]),
                int(sizes[-1]),
                f"{float(np.mean(mean_d_small[:, 0])):.6f}",
                f"{float(np.mean(mean_d_large[:, 0])):.6f}",
                f"{float(np.mean(mean_d_small[:, 1])):.6f}",
                f"{float(np.mean(mean_d_large[:, 1])):.6f}",
                f"{float(np.mean(beta_mat[:, 0])):.6f}",
                f"{float(np.mean(beta_mat[:, 1])):.6f}",
                f"{float(np.mean(divergence)):.6f}",
                f"{float(np.mean(divergence > 0.0)):.4f}",
            ])
    _write_csv(
        out / "diagnostic_probe_sample_size.csv",
        [
            "generator", "n", "dataset_reps", "subsample_reps", "m_min", "m_max",
            "mean_d_pi1_mmin", "mean_d_pi1_mmax", "mean_d_pi090_mmin",
            "mean_d_pi090_mmax", "mean_beta_pi1", "mean_beta_pi090",
            "mean_divergence", "positive_divergence_rate",
        ],
        rows,
    )
    return summary


def write_uniform_draws(path: Path, uniforms: np.ndarray) -> None:
    rows: list[list[object]] = []
    for draw in range(uniforms.shape[0]):
        for case in range(uniforms.shape[1]):
            rows.append([
                draw + 1,
                case + 1,
                f"{float(uniforms[draw, case, 0]):.12f}",
                f"{float(uniforms[draw, case, 1]):.12f}",
            ])
    _write_csv(path, ["draw", "case", "u1", "u2"], rows)


def build_subsample_list(n: int, sizes: list[int], reps: int, seed: int) -> list[np.ndarray]:
    rng = np.random.default_rng(seed)
    subsamples = []
    for size in sizes:
        rows = np.vstack([rng.choice(n, size=size, replace=False) + 1 for _ in range(reps)])
        subsamples.append(rows.astype(int))
    return subsamples


def write_subsamples(path: Path, sizes: list[int], subsamples: list[np.ndarray]) -> None:
    rows: list[list[object]] = []
    for size, matrix in zip(sizes, subsamples):
        for rep in range(matrix.shape[0]):
            for pos in range(matrix.shape[1]):
                rows.append([size, rep + 1, pos + 1, int(matrix[rep, pos])])
    _write_csv(path, ["size", "rep", "pos", "case"], rows)


def parse_sectioned_csv(stdout: str) -> dict[str, list[dict[str, str]]]:
    sections: dict[str, list[str]] = {}
    current: str | None = None
    for raw_line in stdout.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith("[") and line.endswith("]"):
            current = line[1:-1]
            sections[current] = []
        elif current is not None:
            sections[current].append(line)
    parsed: dict[str, list[dict[str, str]]] = {}
    for name, lines in sections.items():
        parsed[name] = list(csv.DictReader(lines))
    return parsed


def run_r_diagnostic_parity(
    data_path: Path,
    uniform_path: Path,
    subsample_path: Path,
    sizes: list[int],
    probe_reps: int,
    n_grid: int,
) -> dict[str, list[dict[str, str]]]:
    code = r'''
args <- commandArgs(trailingOnly = TRUE)
root <- args[[1L]]
data_path <- args[[2L]]
uniform_path <- args[[3L]]
subsample_path <- args[[4L]]
sizes <- as.integer(strsplit(args[[5L]], ",")[[1L]])
probe_reps <- as.integer(args[[6L]])
n_grid <- as.integer(args[[7L]])

source(file.path(root, "r", "R", "qnca.R"))
d <- read.csv(data_path)
u <- read.csv(uniform_path)
M <- max(u$draw)
n <- max(u$case)
uniforms <- array(NA_real_, dim = c(M, n, 2L))
for (i in seq_len(nrow(u))) {
  uniforms[u$draw[i], u$case[i], 1L] <- u$u1[i]
  uniforms[u$draw[i], u$case[i], 2L] <- u$u2[i]
}
s <- read.csv(subsample_path)
subsamples <- vector("list", length(sizes))
for (k in seq_along(sizes)) {
  size_k <- sizes[k]
  part <- s[s$size == size_k, ]
  mat <- matrix(NA_integer_, nrow = probe_reps, ncol = size_k)
  for (i in seq_len(nrow(part))) {
    mat[part$rep[i], part$pos[i]] <- part$case[i]
  }
  subsamples[[k]] <- mat
}
band <- spuriousness_band(
  d$X, d$Y, pi_grid = c(1, 0.95, 0.90), n_grid = n_grid,
  M = M, scope = c(0, 100, 0, 100), uniform_draws = uniforms
)
probe <- consistency_probe(
  d$X, d$Y, pi_pair = c(1, 0.90), sizes = sizes, reps = probe_reps,
  n_grid = n_grid, scope = c(0, 100, 0, 100), subsamples = subsamples
)
cat("[BAND]\n")
write.csv(band$band, row.names = FALSE, quote = FALSE)
cat("[PROBE_SUMMARY]\n")
write.csv(probe$summary, row.names = FALSE, quote = FALSE)
cat("[PROBE_BETA]\n")
cat("pi,beta_drift\n")
for (i in seq_along(probe$pi_values)) {
  cat(sprintf("%.2f,%.15f\n", probe$pi_values[i], probe$beta_drift[i]))
}
cat("[PROBE_DIV]\n")
cat("divergence\n")
cat(sprintf("%.15f\n", probe$divergence))
'''
    proc = subprocess.run(
        [
            "Rscript", "-e", code, str(ROOT), str(data_path), str(uniform_path),
            str(subsample_path), ",".join(str(size) for size in sizes),
            str(probe_reps), str(n_grid),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            "R diagnostic parity run failed\nSTDOUT\n"
            + proc.stdout
            + "\nSTDERR\n"
            + proc.stderr
        )
    return parse_sectioned_csv(proc.stdout)


def run_julia_diagnostic_parity(
    data_path: Path,
    uniform_path: Path,
    subsample_path: Path,
    sizes: list[int],
    probe_reps: int,
    n_grid: int,
) -> dict[str, list[dict[str, str]]]:
    code = r'''
using DelimitedFiles
using Printf
root = ARGS[1]
data_path = ARGS[2]
uniform_path = ARGS[3]
subsample_path = ARGS[4]
sizes = parse.(Int, split(ARGS[5], ","))
probe_reps = parse(Int, ARGS[6])
n_grid = parse(Int, ARGS[7])
include(joinpath(root, "julia", "QuantileNecessity.jl", "src", "QuantileNecessity.jl"))
using .QuantileNecessity

data, _ = readdlm(data_path, ',', header=true)
X = Float64.(data[:, 1])
Y = Float64.(data[:, 2])

uniform_raw, _ = readdlm(uniform_path, ',', header=true)
M = Int(maximum(uniform_raw[:, 1]))
n = Int(maximum(uniform_raw[:, 2]))
uniforms = Array{Float64}(undef, M, n, 2)
for row in eachrow(uniform_raw)
    draw = Int(row[1])
    case = Int(row[2])
    uniforms[draw, case, 1] = Float64(row[3])
    uniforms[draw, case, 2] = Float64(row[4])
end

sub_raw, _ = readdlm(subsample_path, ',', header=true)
subsamples = Vector{Matrix{Int}}(undef, length(sizes))
for (k, size_k) in enumerate(sizes)
    mat = Matrix{Int}(undef, probe_reps, size_k)
    for row in eachrow(sub_raw)
        if Int(row[1]) == size_k
            mat[Int(row[2]), Int(row[3])] = Int(row[4])
        end
    end
    subsamples[k] = mat
end

band = spuriousness_band(
    X, Y; pi_grid=[1.0, 0.95, 0.90], n_grid=n_grid, M=M,
    scope=(0.0, 100.0, 0.0, 100.0), uniform_draws=uniforms
)
probe = consistency_probe(
    X, Y; pi_pair=(1.0, 0.90), sizes=sizes, reps=probe_reps,
    n_grid=n_grid, scope=(0.0, 100.0, 0.0, 100.0), subsamples=subsamples
)

println("[BAND]")
println("pi,d_obs,lower,median,upper,excess,p_null")
for i in eachindex(band.band.pi)
    println(@sprintf(
        "%.2f,%.15f,%.15f,%.15f,%.15f,%.15f,%.15f",
        band.band.pi[i], band.band.d_obs[i], band.band.lower[i],
        band.band.median[i], band.band.upper[i], band.band.excess[i],
        band.band.p_null[i]
    ))
end
println("[PROBE_SUMMARY]")
println("size,pi,mean_d,sd_d")
for k in eachindex(probe.sizes), p in eachindex(probe.pi_values)
    println(@sprintf(
        "%d,%.2f,%.15f,%.15f",
        probe.sizes[k], probe.pi_values[p], probe.summary.mean_d[(p - 1) * length(probe.sizes) + k],
        probe.summary.sd_d[(p - 1) * length(probe.sizes) + k]
    ))
end
println("[PROBE_BETA]")
println("pi,beta_drift")
for i in eachindex(probe.pi_values)
    println(@sprintf("%.2f,%.15f", probe.pi_values[i], probe.beta_drift[i]))
end
println("[PROBE_DIV]")
println("divergence")
println(@sprintf("%.15f", probe.divergence))
'''
    proc = subprocess.run(
        [
            "julia", f"--project={ROOT / 'julia' / 'QuantileNecessity.jl'}", "-e", code, str(ROOT),
            str(data_path), str(uniform_path), str(subsample_path),
            ",".join(str(size) for size in sizes), str(probe_reps), str(n_grid),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            "Julia diagnostic parity run failed\nSTDOUT\n"
            + proc.stdout
            + "\nSTDERR\n"
            + proc.stderr
        )
    return parse_sectioned_csv(proc.stdout)


def max_pairwise_abs(values: list[np.ndarray]) -> float:
    max_diff = 0.0
    for i in range(len(values)):
        for j in range(i + 1, len(values)):
            max_diff = max(max_diff, float(np.max(np.abs(values[i] - values[j]))))
    return max_diff


def probe_summary_matrix(
    rows: list[dict[str, str]],
    sizes: list[int],
    pi_values: tuple[float, float],
    field: str,
) -> np.ndarray:
    out = np.empty((len(sizes), len(pi_values)), dtype=float)
    for row in rows:
        size_idx = sizes.index(int(row["size"]))
        pi_val = round(float(row["pi"]), 2)
        pi_idx = [round(float(pi), 2) for pi in pi_values].index(pi_val)
        out[size_idx, pi_idx] = float(row[field])
    return out


def run_diagnostic_parity(
    out: Path,
    seed: int,
) -> dict[str, float]:
    rows: list[list[object]] = []
    overall = {
        "band_max": 0.0,
        "probe_beta_max": 0.0,
    }
    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        for generator in ("floor", "no_floor"):
            n = 200
            X, Y = make_dataset(generator, n, diagnostic_seed(seed, generator, n, 21))
            data_path = tmpdir / f"{generator}_data.csv"
            _write_csv(
                data_path,
                ["X", "Y"],
                [[f"{float(x):.12f}", f"{float(y):.12f}"] for x, y in zip(X, Y)],
            )
            uniform_seed = diagnostic_seed(seed + 2_000_000, generator, n, 0)
            uniforms = np.random.default_rng(uniform_seed).uniform(0.0, 1.0, size=(19, n, 2))
            uniform_path = tmpdir / f"{generator}_uniforms.csv"
            write_uniform_draws(uniform_path, uniforms)
            sizes = [50, 100, 150, 200]
            subsamples = build_subsample_list(n, sizes, reps=12, seed=diagnostic_seed(seed + 3_000_000, generator, n, 0))
            subsample_path = tmpdir / f"{generator}_subsamples.csv"
            write_subsamples(subsample_path, sizes, subsamples)

            band_py = spuriousness_band(
                X, Y, pi_grid=PI_GRID, n_grid=40, M=uniforms.shape[0], scope=SCOPE,
                uniform_draws=uniforms,
            )
            probe_py = consistency_probe(
                X, Y, pi_pair=PI_PAIR, sizes=sizes, reps=12, n_grid=40, scope=SCOPE,
                subsamples=subsamples,
            )
            rr = run_r_diagnostic_parity(data_path, uniform_path, subsample_path, sizes, 12, 40)
            jl = run_julia_diagnostic_parity(data_path, uniform_path, subsample_path, sizes, 12, 40)

            band_r = rr["BAND"]
            band_j = jl["BAND"]
            band_fields = ("d_obs", "lower", "median", "upper", "excess", "p_null")
            band_diffs: dict[str, float] = {}
            for field in band_fields:
                arrays = [
                    np.asarray(getattr(band_py, field), dtype=float),
                    np.asarray([float(row[field]) for row in band_r], dtype=float),
                    np.asarray([float(row[field]) for row in band_j], dtype=float),
                ]
                band_diffs[field] = max_pairwise_abs(arrays)

            probe_summary_r = rr["PROBE_SUMMARY"]
            probe_summary_j = jl["PROBE_SUMMARY"]
            mean_d_arrays = [
                np.asarray(probe_py.mean_d, dtype=float),
                probe_summary_matrix(probe_summary_r, sizes, PI_PAIR, "mean_d"),
                probe_summary_matrix(probe_summary_j, sizes, PI_PAIR, "mean_d"),
            ]
            sd_d_arrays = [
                np.asarray(probe_py.sd_d, dtype=float),
                probe_summary_matrix(probe_summary_r, sizes, PI_PAIR, "sd_d"),
                probe_summary_matrix(probe_summary_j, sizes, PI_PAIR, "sd_d"),
            ]
            beta_arrays = [
                np.asarray(probe_py.beta_drift, dtype=float),
                np.asarray([float(row["beta_drift"]) for row in rr["PROBE_BETA"]], dtype=float),
                np.asarray([float(row["beta_drift"]) for row in jl["PROBE_BETA"]], dtype=float),
            ]
            divergence_arrays = [
                np.asarray([probe_py.divergence], dtype=float),
                np.asarray([float(rr["PROBE_DIV"][0]["divergence"])], dtype=float),
                np.asarray([float(jl["PROBE_DIV"][0]["divergence"])], dtype=float),
            ]

            rows.append([
                generator,
                n,
                uniforms.shape[0],
                12,
                f"{band_diffs['d_obs']:.2e}",
                f"{band_diffs['lower']:.2e}",
                f"{band_diffs['median']:.2e}",
                f"{band_diffs['upper']:.2e}",
                f"{band_diffs['excess']:.2e}",
                f"{band_diffs['p_null']:.2e}",
                f"{max_pairwise_abs(mean_d_arrays):.2e}",
                f"{max_pairwise_abs(sd_d_arrays):.2e}",
                f"{max_pairwise_abs(beta_arrays):.2e}",
                f"{max_pairwise_abs(divergence_arrays):.2e}",
            ])
            overall["band_max"] = max(overall["band_max"], max(band_diffs.values()))
            overall["probe_beta_max"] = max(overall["probe_beta_max"], max_pairwise_abs(beta_arrays))

    _write_csv(
        out / "diagnostic_cross_language_agreement.csv",
        [
            "generator", "n", "M", "probe_reps", "band_max_abs_d_obs_diff",
            "band_max_abs_lower_diff", "band_max_abs_median_diff",
            "band_max_abs_upper_diff", "band_max_abs_excess_diff",
            "band_max_abs_p_null_diff", "probe_max_abs_mean_d_diff",
            "probe_max_abs_sd_d_diff", "probe_max_abs_beta_diff",
            "probe_max_abs_divergence_diff",
        ],
        rows,
    )
    return overall


def run_calibration_reconciliation(
    out: Path,
    seed: int,
    workers: int,
) -> dict[tuple[str, int, float], float]:
    rows: list[list[object]] = []
    summary: dict[tuple[str, int, float], float] = {}
    with (out / "calibration_null.csv").open(newline="") as f:
        for row in csv.DictReader(f):
            pi = round(float(row["pi"]), 2)
            rate = float(row["fpr"])
            summary[("published_calibration", 150, pi)] = rate
            rows.append([
                "published_calibration",
                "data/calibration_null.csv",
                "independent_uniform",
                150,
                f"{pi:.2f}",
                int(row["reps"]),
                999,
                int(row["rejections_at_0.05"]),
                f"{rate:.4f}",
            ])
    with (out / "power_curve.csv").open(newline="") as f:
        for row in csv.DictReader(f):
            if abs(float(row["ceiling_strength"]) - 0.0) > 1e-12:
                continue
            pi = 0.95
            rate = float(row["power"])
            summary[("published_strength0", 120, pi)] = rate
            rows.append([
                "published_strength0",
                "data/power_curve.csv",
                "strength_zero",
                120,
                f"{pi:.2f}",
                int(row["reps"]),
                int(row["B"]),
                int(row["rejections"]),
                f"{rate:.4f}",
            ])
            break

    scenarios = [
        ("harmonized_strength0", "strength_zero", 150, 300, 999),
        ("harmonized_uniform120", "independent_uniform", 120, 150, 499),
    ]
    for scenario, design, n, reps, B in scenarios:
        tasks = [
            (design, n, pi, B, seed + 4_000_000 + n * 10_000 + rep * 10 + int(round(pi * 100)), 40)
            for pi in PI_GRID
            for rep in range(reps)
        ]
        results = _run_tasks(calibration_task, tasks, workers)
        for pi in PI_GRID:
            rejects = sum(
                1 for design_i, n_i, pi_i, reject in results
                if design_i == design and n_i == n and abs(pi_i - pi) < 1e-12 and reject
            )
            rate = rejects / reps
            summary[(scenario, n, pi)] = rate
            rows.append([scenario, design, n, f"{pi:.2f}", reps, B, rejects, f"{rate:.4f}"])
    _write_csv(
        out / "calibration_reconciliation.csv",
        ["scenario", "source", "design", "n", "pi", "reps", "B", "rejections_at_0.05", "rejection_rate"],
        rows,
    )
    return summary


def build_log(
    illustration: dict[str, dict[float, dict[str, float]]],
    band_summary: dict[tuple[str, int], dict[float, dict[str, float]]],
    probe_summary: dict[tuple[str, int], dict[str, float]],
    parity_summary: dict[str, float],
    calibration_summary: dict[tuple[str, int, float], float],
    seed: int,
    dataset_reps: int,
    sample_M: int,
    illustration_M: int,
    probe_reps: int,
) -> str:
    band_gate = all(
        illustration["floor"][pi]["delta"] > 0.0 and illustration["floor"][pi]["p_null"] <= 0.05
        for pi in PI_GRID
    ) and all(
        illustration["no_floor"][pi]["delta"] <= 0.0 or illustration["no_floor"][pi]["p_null"] > 0.05
        for pi in PI_GRID
    ) and all(
        band_summary[("floor", n)][0.95]["detection_rate"] > band_summary[("no_floor", n)][0.95]["detection_rate"]
        and band_summary[("no_floor", n)][0.95]["detection_rate"] <= 0.20
        for n in SIZES
    ) and all(
        band_summary[("floor", n)][0.95]["detection_rate"] >= 0.80
        for n in (200, 500, 1000)
    )

    probe_gate = all(
        probe_summary[("floor", n)]["mean_divergence"] > 0.0
        and probe_summary[("no_floor", n)]["mean_divergence"] <= 0.0
        for n in SIZES
    )

    lines = [
        "QNCA diagnostic study",
        "=====================",
        f"seed={seed}",
        f"band_illustration_M={illustration_M}",
        f"band_sample_M={sample_M}",
        f"dataset_reps={dataset_reps}",
        f"probe_subsample_reps={probe_reps}",
        "",
        "Spuriousness band illustration",
        "-----------------------------",
    ]
    for generator in ("floor", "no_floor"):
        for pi in PI_GRID:
            item = illustration[generator][pi]
            lines.append(
                f"{generator} pi={pi:.2f} d_obs={item['d_obs']:.6f} upper={item['upper']:.6f} "
                f"delta={item['delta']:.6f} p_null={item['p_null']:.6f}"
            )

    lines.extend([
        "",
        "Band sample-size screen at pi=0.95",
        "---------------------------------",
    ])
    for n in SIZES:
        floor_rate = band_summary[("floor", n)][0.95]["detection_rate"]
        no_floor_rate = band_summary[("no_floor", n)][0.95]["detection_rate"]
        lines.append(
            f"n={n} floor_detection={floor_rate:.4f} no_floor_detection={no_floor_rate:.4f}"
        )

    lines.extend([
        "",
        "Consistency probe screen",
        "------------------------",
    ])
    for n in SIZES:
        floor_item = probe_summary[("floor", n)]
        no_floor_item = probe_summary[("no_floor", n)]
        lines.append(
            f"n={n} floor_divergence={floor_item['mean_divergence']:.6f} "
            f"no_floor_divergence={no_floor_item['mean_divergence']:.6f} "
            f"no_floor_positive_rate={no_floor_item['positive_divergence_rate']:.4f}"
        )

    lines.extend([
        "",
        "Cross-language parity",
        "--------------------",
        f"band_max_abs_diff={parity_summary['band_max']:.2e}",
        f"probe_beta_max_abs_diff={parity_summary['probe_beta_max']:.2e}",
        "",
        "Calibration reconciliation",
        "-------------------------",
    ])
    rates = ", ".join(
        f"pi={pi:.2f}:{calibration_summary[('published_calibration', 150, pi)]:.4f}"
        for pi in PI_GRID
    )
    lines.append(f"published_calibration design=independent_uniform n=150 {rates}")
    lines.append(
        f"published_strength0 design=strength_zero n=120 pi=0.95:"
        f"{calibration_summary[('published_strength0', 120, 0.95)]:.4f}"
    )
    for scenario, design, n in [
        ("harmonized_strength0", "strength_zero", 150),
        ("harmonized_uniform120", "independent_uniform", 120),
    ]:
        rates = ", ".join(
            f"pi={pi:.2f}:{calibration_summary[(scenario, n, pi)]:.4f}"
            for pi in PI_GRID
        )
        lines.append(f"{scenario} design={design} n={n} {rates}")

    lines.extend([
        "",
        "Synthetic verdicts",
        "-----------------",
        f"band_gate={'PASS' if band_gate else 'FAIL'}",
        f"probe_gate={'PASS' if probe_gate else 'FAIL'}",
    ])
    if not probe_gate:
        lines.append(
            "probe_reading=Open problem: no-floor data still show positive divergence, "
            "consistent with top-band, support-bound, and skewness confounds."
        )
    if calibration_summary[("published_calibration", 150, 1.00)] > 0.05 or calibration_summary[("harmonized_strength0", 150, 1.00)] > 0.05:
        lines.append(
            "calibration_reading=pi=1 remains anti-conservative in at least one small-sample null design."
        )
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=20260609)
    parser.add_argument("--illustration-M", type=int, default=999)
    parser.add_argument("--sample-M", type=int, default=199)
    parser.add_argument("--dataset-reps", type=int, default=15)
    parser.add_argument("--probe-reps", type=int, default=100)
    parser.add_argument("--workers", type=int, default=max(1, min(4, os.cpu_count() or 1)))
    parser.add_argument("--out", type=Path, default=ROOT / "data")
    args = parser.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    illustration = run_band_illustration(args.out, args.seed, args.illustration_M, n_grid=50)
    band_summary = run_band_sample_size(
        args.out, args.seed, args.dataset_reps, args.sample_M, n_grid=50, workers=args.workers
    )
    probe_summary = run_probe_sample_size(
        args.out, args.seed, args.dataset_reps, args.probe_reps, n_grid=50, workers=args.workers
    )
    parity_summary = run_diagnostic_parity(args.out, args.seed)
    calibration_summary = run_calibration_reconciliation(args.out, args.seed, workers=args.workers)
    log_text = build_log(
        illustration, band_summary, probe_summary, parity_summary,
        calibration_summary, args.seed, args.dataset_reps, args.sample_M,
        args.illustration_M, args.probe_reps,
    )
    (args.out / "qnca_diagnostics.log").write_text(log_text)
    print(f"wrote diagnostic outputs to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
