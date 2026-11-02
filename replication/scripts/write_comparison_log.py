#!/usr/bin/env python3
"""Write one consolidated QNCA comparison log and bottleneck comparison table."""

from __future__ import annotations

import argparse
import csv
import io
import subprocess
from pathlib import Path

import numpy as np

from qnca import qnca

ROOT = Path(__file__).resolve().parent.parent.parent
SCOPE = (0.0, 100.0, 0.0, 100.0)
SIZES = (100, 200, 500, 1000)
PI_GRID = (1.00, 0.95, 0.90)
BOTTLENECK_LEVELS = tuple(range(0, 101, 10))
PAPER_BOTTLENECK_LEVELS = (0, 20, 50, 70, 80)


def read_text(path: Path) -> str:
    return path.read_text().strip()


def extract_csv(stdout: str, marker: str) -> list[dict[str, str]]:
    lines = stdout.splitlines()
    try:
        start = lines.index(marker) + 1
    except ValueError as exc:
        raise RuntimeError(f"Missing marker {marker!r}") from exc
    reader = csv.DictReader(lines[start:])
    return list(reader)


def python_bottlenecks(n_grid: int) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for n in SIZES:
        data = np.loadtxt(ROOT / "data" / f"qnca_dataset_n{n}.csv", delimiter=",", skiprows=1)
        X, Y = data[:, 0], data[:, 1]
        for pi in PI_GRID:
            res = qnca(X, Y, pi=pi, n_grid=n_grid, B=None, scope=SCOPE)
            by_level = {
                int(round(y)): x
                for y, x in zip(res.y_grid, res.x_required)
                if int(round(y)) in BOTTLENECK_LEVELS and abs(y - round(y)) < 1e-9
            }
            for level in BOTTLENECK_LEVELS:
                rows.append({
                    "dataset": f"n{n}",
                    "n": str(n),
                    "pi": f"{pi:.2f}",
                    "outcome_level": str(level),
                    "python": f"{by_level[level]:.12f}",
                })
    return rows


def r_bottlenecks(n_grid: int) -> list[dict[str, str]]:
    code = r'''
args <- commandArgs(trailingOnly = TRUE)
root <- args[[1L]]
n_grid <- as.integer(args[[2L]])
sizes <- as.integer(strsplit(args[[3L]], ",")[[1L]])
pis <- as.numeric(strsplit(args[[4L]], ",")[[1L]])
source(file.path(root, "r", "R", "qnca.R"))
rows <- list(); idx <- 1L
for (n in sizes) {
  d <- read.csv(file.path(root, "data", paste0("qnca_dataset_n", n, ".csv")))
  for (pi in pis) {
    res <- qnca(d$X, d$Y, pi = pi, n_grid = n_grid, B = NA,
                scope = c(0, 100, 0, 100))
    keep <- res$bottleneck$outcome_level %in% seq(0, 100, 10)
    tab <- res$bottleneck[keep, ]
    for (i in seq_len(nrow(tab))) {
      rows[[idx]] <- data.frame(dataset = paste0("n", n), n = n,
                                pi = sprintf("%.2f", pi),
                                outcome_level = as.integer(tab$outcome_level[i]),
                                R = sprintf("%.12f", tab$X_required[i]))
      idx <- idx + 1L
    }
  }
}
cat("QNCA_R_BOTTLENECKS\n")
write.csv(do.call(rbind, rows), row.names = FALSE)
'''
    proc = subprocess.run(
        [
            "Rscript", "-e", code, str(ROOT), str(n_grid),
            ",".join(str(n) for n in SIZES),
            ",".join(f"{pi:.2f}" for pi in PI_GRID),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError("R bottleneck run failed\n" + proc.stdout + "\n" + proc.stderr)
    return extract_csv(proc.stdout, "QNCA_R_BOTTLENECKS")


def julia_bottlenecks(n_grid: int) -> list[dict[str, str]]:
    code = r'''
using DelimitedFiles
using Printf
root = ARGS[1]
n_grid = parse(Int, ARGS[2])
sizes = parse.(Int, split(ARGS[3], ","))
pis = parse.(Float64, split(ARGS[4], ","))
include(joinpath(root, "julia", "QuantileNecessity.jl", "src", "QuantileNecessity.jl"))
using .QuantileNecessity
println("QNCA_JULIA_BOTTLENECKS")
println("dataset,n,pi,outcome_level,julia")
for n in sizes
    data, _ = readdlm(joinpath(root, "data", "qnca_dataset_n$(n).csv"), ',', header=true)
    X = Float64.(data[:, 1])
    Y = Float64.(data[:, 2])
    for pi in pis
        res = qnca(X, Y; pi=pi, n_grid=n_grid, B=0, scope=(0.0, 100.0, 0.0, 100.0))
        for (y, x) in zip(res.y_grid, res.x_required)
            if isapprox(y, round(y); atol=1e-9) && round(Int, y) % 10 == 0
                println(@sprintf("n%d,%d,%.2f,%d,%.12f", n, n, pi, round(Int, y), x))
            end
        end
    end
end
'''
    proc = subprocess.run(
        [
            "julia", f"--project={ROOT / 'julia' / 'QuantileNecessity.jl'}", "-e", code, str(ROOT), str(n_grid),
            ",".join(str(n) for n in SIZES),
            ",".join(f"{pi:.2f}" for pi in PI_GRID),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError("Julia bottleneck run failed\n" + proc.stdout + "\n" + proc.stderr)
    return extract_csv(proc.stdout, "QNCA_JULIA_BOTTLENECKS")


def dul_bottlenecks() -> list[dict[str, str]]:
    code = r'''
args <- commandArgs(trailingOnly = TRUE)
root <- args[[1L]]
suppressPackageStartupMessages(library(NCA))
sizes <- c(100, 200, 500, 1000)
ceilings_use <- c("ce_fdh", "cr_fdh")
rows <- list(); idx <- 1L
for (n in sizes) {
  d <- read.csv(file.path(root, "data", paste0("qnca_dataset_n", n, ".csv")))
  model <- nca_analysis(
    data = d, x = "X", y = "Y",
    ceilings = ceilings_use,
    scope = c(0, 100, 0, 100),
    bottleneck.x = "actual", bottleneck.y = "actual",
    steps = seq(0, 100, 10),
    test.rep = 0
  )
  for (ceiling in ceilings_use) {
    tab <- model$bottlenecks[[ceiling]]
    for (i in seq_len(nrow(tab))) {
      x_val <- tab$X[i]
      x_txt <- if (is.na(x_val)) {
        "NA"
      } else if (is.numeric(x_val)) {
        sprintf("%.12f", x_val)
      } else {
        as.character(x_val)
      }
      rows[[idx]] <- data.frame(dataset = paste0("n", n), n = n,
                                method = ceiling,
                                outcome_level = as.integer(tab$Y[i]),
                                dul_X_required = x_txt)
      idx <- idx + 1L
    }
  }
}
cat("DUL_BOTTLENECKS\n")
write.csv(do.call(rbind, rows), row.names = FALSE)
'''
    proc = subprocess.run(
        ["Rscript", "-e", code, str(ROOT)],
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError("Dul bottleneck run failed\n" + proc.stdout + "\n" + proc.stderr)
    return extract_csv(proc.stdout, "DUL_BOTTLENECKS")


def merge_qnca_bottlenecks(
    py: list[dict[str, str]],
    rr: list[dict[str, str]],
    jl: list[dict[str, str]],
) -> list[dict[str, str]]:
    key = ("dataset", "n", "pi", "outcome_level")
    rows = {tuple(row[k] for k in key): dict(row) for row in py}
    for source in (rr, jl):
        for row in source:
            ident = tuple(row[k] for k in key)
            rows.setdefault(ident, {k: row[k] for k in key}).update(row)
    merged = []
    for ident in sorted(rows, key=lambda x: (int(x[1]), float(x[2]), int(x[3]))):
        row = rows[ident]
        vals = [float(row[name]) for name in ("python", "R", "julia")]
        row["max_abs_diff"] = f"{max(vals) - min(vals):.2e}"
        merged.append(row)
    return merged


def csv_block(rows: list[dict[str, str]], fieldnames: list[str]) -> str:
    out = io.StringIO()
    writer = csv.DictWriter(out, fieldnames=fieldnames, extrasaction="ignore")
    writer.writeheader()
    writer.writerows(rows)
    return out.getvalue().strip()


def load_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def build_paper_bottleneck_table(
    qnca_bn: list[dict[str, str]],
    dul_bn: list[dict[str, str]],
) -> list[dict[str, str]]:
    qnca_lookup = {
        (row["dataset"], row["pi"], int(row["outcome_level"])): row["R"]
        for row in qnca_bn
    }
    dul_lookup = {
        (row["dataset"], int(row["outcome_level"])): row["dul_X_required"]
        for row in dul_bn
        if row["method"] == "ce_fdh"
    }
    rows: list[dict[str, str]] = []
    dataset = "n1000"
    for level in PAPER_BOTTLENECK_LEVELS:
        rows.append({
            "dataset": dataset,
            "n": "1000",
            "outcome_level": str(level),
            "dul_ce_fdh": dul_lookup[(dataset, level)],
            "qnca_pi1_00": qnca_lookup[(dataset, "1.00", level)],
            "qnca_pi0_95": qnca_lookup[(dataset, "0.95", level)],
            "qnca_pi0_90": qnca_lookup[(dataset, "0.90", level)],
        })
    return rows


def render_paper_bottleneck_tabular(rows: list[dict[str, str]]) -> str:
    body = []
    for row in rows:
        body.append(
            " & ".join([
                row["outcome_level"],
                f"{float(row['dul_ce_fdh']):.3f}",
                f"{float(row['qnca_pi1_00']):.3f}",
                f"{float(row['qnca_pi0_95']):.3f}",
                f"{float(row['qnca_pi0_90']):.3f}",
            ]) + r" \\"
        )
    return "\n".join([
        r"\begin{tabular}{rrrrr}",
        r"\toprule",
        r"Outcome level & Dul CE-FDH & QNCA \(\pi=1.00\) & QNCA \(\pi=0.95\) & QNCA \(\pi=0.90\) \\",
        r"\midrule",
        *body,
        r"\bottomrule",
        r"\end{tabular}",
        "",
    ])


def write_log(path: Path, bottleneck_grid: int) -> None:
    cross = load_csv(ROOT / "data" / "cross_language_agreement.csv")
    dul = load_csv(ROOT / "data" / "dul_nca_results.csv")
    py = python_bottlenecks(bottleneck_grid)
    rr = r_bottlenecks(bottleneck_grid)
    jl = julia_bottlenecks(bottleneck_grid)
    qnca_bn = merge_qnca_bottlenecks(py, rr, jl)
    dul_bn = dul_bottlenecks()
    paper_bn = build_paper_bottleneck_table(qnca_bn, dul_bn)

    max_d_diff = max(float(row["max_abs_d_diff"]) for row in cross)
    p_equal = all(row["p_values_equal"].strip().lower() == "true" for row in cross)
    max_dul_diff = max(
        float(row["abs_diff_native_vs_dul"])
        for row in dul
        if row["method"] == "ce_fdh"
    )
    max_bn_diff = max(float(row["max_abs_diff"]) for row in qnca_bn)

    text = f"""QNCA comparison log

Sources
scripts/run_cross_language_agreement.py
scripts/run_study.py
scripts/write_comparison_log.py
data/qnca_dataset_n100.csv
data/qnca_dataset_n200.csv
data/qnca_dataset_n500.csv
data/qnca_dataset_n1000.csv

Parameters
scope = [0, 100] x [0, 100]
pi_grid = 1.00, 0.95, 0.90
cross_language_B = {cross[0]["B"]}
cross_language_n_grid = {cross[0]["n_grid"]}
bottleneck_n_grid = {bottleneck_grid}
bottleneck_levels = 0, 10, ..., 100
Dul_NCA_test_rep = {dul[0]["test_rep"]}
Dul_NCA_version = {dul[0]["nca_version"]}
R_version = {dul[0]["r_version"]}

Formulation check
The quantile frontier is Q_(1-pi){{X_i : Y_i >= y}} with R type-1 quantiles.
At pi = 1, the quantile is the minimum and the native CE-FDH effect size matches
Dul NCA CE-FDH under the fixed scope. The d formula remains empty-zone area over
scope, but d_pi is pi-relative because lowering pi raises the frontier and
enlarges the empty zone. The p formula remains the finite-sample permutation
tail probability (1 + count)/(B + 1). The whole frontier, isotonic projection,
and d_pi are recomputed inside every permutation.

Gate summary
max_cross_language_d_diff = {max_d_diff:.3e}
cross_language_p_values_equal = {str(p_equal).lower()}
max_cross_language_bottleneck_diff = {max_bn_diff:.3e}
max_native_vs_Dul_CE_FDH_d_diff = {max_dul_diff:.3e}

QNCA d and p comparison
{csv_block(cross, ["dataset", "n", "pi", "B", "n_grid", "d_pi_python", "d_pi_R", "d_pi_julia", "max_abs_d_diff", "p_pi_python", "p_pi_R", "p_pi_julia", "p_values_equal"])}

Dul NCA d and p comparison
{csv_block(dul, ["dataset", "n", "method", "qnca_d_pi1", "native_ce_fdh_d", "dul_d", "dul_p", "abs_diff_qnca_vs_dul", "abs_diff_native_vs_dul", "test_rep", "qnca_n_grid", "r_version", "nca_version"])}

QNCA bottleneck comparison
{csv_block(qnca_bn, ["dataset", "n", "pi", "outcome_level", "python", "R", "julia", "max_abs_diff"])}

Dul NCA bottleneck tables
{csv_block(dul_bn, ["dataset", "n", "method", "outcome_level", "dul_X_required"])}
"""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    with (ROOT / "data" / "bottleneck_comparison_n1000.csv").open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "dataset", "n", "outcome_level",
                "dul_ce_fdh", "qnca_pi1_00", "qnca_pi0_95", "qnca_pi0_90",
            ],
        )
        writer.writeheader()
        writer.writerows(paper_bn)
    tables_dir = ROOT / "replication" / "figures"
    tables_dir.mkdir(parents=True, exist_ok=True)
    (tables_dir / "bottleneck_comparison_n1000_tabular.tex").write_text(
        render_paper_bottleneck_tabular(paper_bn)
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=ROOT / "data" / "qnca_comparison.log")
    parser.add_argument("--bottleneck-grid", type=int, default=101)
    args = parser.parse_args()
    write_log(args.out, args.bottleneck_grid)
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
