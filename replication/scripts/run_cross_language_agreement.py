#!/usr/bin/env python3
"""Run QNCA agreement checks across Python, R and Julia."""

from __future__ import annotations

import argparse
import csv
import subprocess
import tempfile
from pathlib import Path

import numpy as np

from qnca import qnca

ROOT = Path(__file__).resolve().parent.parent.parent
SCOPE = (0.0, 100.0, 0.0, 100.0)
SIZES = (100, 200, 500, 1000)
PI_GRID = (1.00, 0.95, 0.90)


def write_permutations(path: Path, n: int, B: int, seed: int) -> None:
    rng = np.random.default_rng(seed)
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        for _ in range(B):
            writer.writerow((rng.permutation(n) + 1).tolist())


def read_python(data_path: Path, perm_path: Path, n_grid: int) -> dict[float, tuple[float, float]]:
    data = np.loadtxt(data_path, delimiter=",", skiprows=1)
    permutations = np.loadtxt(perm_path, delimiter=",", dtype=int)
    if permutations.ndim == 1:
        permutations = permutations.reshape(1, -1)
    X, Y = data[:, 0], data[:, 1]
    out = {}
    for pi in PI_GRID:
        res = qnca(X, Y, pi=pi, n_grid=n_grid, B=None, scope=SCOPE, permutations=permutations)
        out[pi] = (float(res.d_pi), float(res.p_pi))
    return out


def parse_table(stdout: str) -> dict[float, tuple[float, float]]:
    rows = []
    keep = False
    for line in stdout.splitlines():
        if line.startswith("pi,"):
            keep = True
            rows.append(line)
        elif keep and line[:1].isdigit():
            rows.append(line)
    reader = csv.DictReader(rows)
    return {
        float(row["pi"]): (float(row["d_pi"]), float(row["p_pi"]))
        for row in reader
    }


def run_r(data_path: Path, perm_path: Path, n_grid: int) -> dict[float, tuple[float, float]]:
    code = r'''
args <- commandArgs(trailingOnly = TRUE)
root <- args[[1L]]
data_path <- args[[2L]]
perm_path <- args[[3L]]
n_grid <- as.integer(args[[4L]])
pis <- as.numeric(strsplit(args[[5L]], ",")[[1L]])
source(file.path(root, "r", "R", "datasets.R"))
source(file.path(root, "r", "R", "qnca.R"))
d <- read.csv(data_path)
perm <- as.matrix(read.csv(perm_path, header = FALSE))
storage.mode(perm) <- "integer"
cat("pi,d_pi,p_pi\n")
for (pi in pis) {
  res <- qnca(d$X, d$Y, pi = pi, n_grid = n_grid, B = NA,
              scope = c(0, 100, 0, 100), permutations = perm)
  cat(sprintf("%.2f,%.12f,%.12f\n", pi, res$d_pi, res$p_pi))
}
'''
    cmd = [
        "Rscript", "-e", code, str(ROOT), str(data_path), str(perm_path), str(n_grid),
        ",".join(f"{pi:.2f}" for pi in PI_GRID),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        raise RuntimeError(
            "R agreement run failed\nSTDOUT\n"
            + proc.stdout
            + "\nSTDERR\n"
            + proc.stderr
        )
    return parse_table(proc.stdout)


def run_julia(data_path: Path, perm_path: Path, n_grid: int) -> dict[float, tuple[float, float]]:
    code = r'''
using DelimitedFiles
using Printf
root = ARGS[1]
data_path = ARGS[2]
perm_path = ARGS[3]
n_grid = parse(Int, ARGS[4])
pis = parse.(Float64, split(ARGS[5], ","))
include(joinpath(root, "julia", "QuantileNecessity.jl", "src", "QuantileNecessity.jl"))
using .QuantileNecessity
data, _ = readdlm(data_path, ',', header=true)
X = Float64.(data[:, 1])
Y = Float64.(data[:, 2])
perm_raw = readdlm(perm_path, ',', Int)
perm = ndims(perm_raw) == 1 ? reshape(perm_raw, 1, :) : perm_raw
println("pi,d_pi,p_pi")
for pi in pis
    res = qnca(X, Y; pi=pi, n_grid=n_grid, B=0, scope=(0.0, 100.0, 0.0, 100.0),
               permutations=perm)
    println(@sprintf("%.2f,%.12f,%.12f", pi, res.d_pi, res.p_pi))
end
'''
    cmd = [
        "julia", f"--project={ROOT / 'julia' / 'QuantileNecessity.jl'}", "-e", code, str(ROOT), str(data_path),
        str(perm_path), str(n_grid), ",".join(f"{pi:.2f}" for pi in PI_GRID),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        raise RuntimeError(
            "Julia agreement run failed\nSTDOUT\n"
            + proc.stdout
            + "\nSTDERR\n"
            + proc.stderr
        )
    return parse_table(proc.stdout)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--B", type=int, default=199)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--n-grid", type=int, default=50)
    parser.add_argument("--out", type=Path, default=ROOT / "data" / "cross_language_agreement.csv")
    args = parser.parse_args()

    rows = []
    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        for n in SIZES:
            data_path = ROOT / "data" / f"qnca_dataset_n{n}.csv"
            perm_path = tmpdir / f"permutations_n{n}.csv"
            write_permutations(perm_path, n, args.B, args.seed + n)
            py = read_python(data_path, perm_path, args.n_grid)
            rr = run_r(data_path, perm_path, args.n_grid)
            jl = run_julia(data_path, perm_path, args.n_grid)
            for pi in PI_GRID:
                d_values = (py[pi][0], rr[pi][0], jl[pi][0])
                p_values = (py[pi][1], rr[pi][1], jl[pi][1])
                rows.append([
                    f"n{n}", n, f"{pi:.2f}", args.B, args.n_grid,
                    f"{py[pi][0]:.12f}", f"{rr[pi][0]:.12f}", f"{jl[pi][0]:.12f}",
                    f"{max(d_values) - min(d_values):.2e}",
                    f"{py[pi][1]:.12f}", f"{rr[pi][1]:.12f}", f"{jl[pi][1]:.12f}",
                    str(p_values[0] == p_values[1] == p_values[2]).lower(),
                ])

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "dataset", "n", "pi", "B", "n_grid",
            "d_pi_python", "d_pi_R", "d_pi_julia", "max_abs_d_diff",
            "p_pi_python", "p_pi_R", "p_pi_julia", "p_values_equal",
        ])
        writer.writerows(rows)
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
